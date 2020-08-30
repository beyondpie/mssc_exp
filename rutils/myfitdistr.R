is_outlier <- function(x, up_prob = 0.995) {
  return(x > 5 * quantile(x, up_prob))
}

empirical_prob_zero <- function(x, rmoutliers = T) {
  if (rmoutliers) {
    outliers <- is_outlier(x)
    x <- x[!outliers]
  }
  invisible(1 - length(which(x > 0)) / length(x))
}

fit_poi_with_scalefactors <- function(x, s, rmoutliers = T) {
  if (length(x) != length(s)) {
    stop("unequal length of x and s.")
  }

  if (rmoutliers) {
    outliers <- is_outlier(x)
    x <- x[!outliers]
    s <- s[!outliers]
  }
  sumcnts <- sum(x)
  sumfactors <- sum(s)
  if (sumfactors < 1) {
    stop("all of s are zeros ?")
  }

  estimate <- sumcnts / sumfactors
  sds <- sqrt(estimate / sumfactors)
  names(estimate) <- names(sds) <- "lambda"
  vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list(
    "lambda",
    "lambda"
  ))
  loglik <- sum(dpois(x, lambda = s * estimate, log = TRUE))
  invisible(structure(list(
    estimate = estimate, sd = sds,
    vcov = vc, n = length(x), loglik = loglik
  ), class = "fitdistr"))
}

prob_zero_poi_scale <- function(x, s, rmoutliers = T) {
  thefit <- fit_poi_with_scalefactors(x, s, rmoutliers)
  scaled_mean <- thefit$estimate
  p0 <- exp(-scaled_mean * s)
  invisible(list(
    fitpoiscl = thefit,
    mnp0 = mean(p0),
    mdnp0 = median(p0)
  ))
}

prob_zero_poi <- function(x, rmoutliers = T) {
  if (rmoutliers) {
    outliers <- is_outlier(x)
    x <- x[!outliers]
  }

  thefit <- MASS::fitdistr(x, densfun = "poisson")
  p0 <- exp(-thefit$estimate)
  invisible(list(fitpoi = thefit, p0 = p0))
}

prob_zero_nb <- function(x, rmoutliers = T) {
  if (rmoutliers) {
    outliers <- is_outlier(x)
    x <- x[!outliers]
  }
  nbfit <- MASS::fitdistr(x, densfun = "negative binomial")
  ## r or its reciprocal 1/r is also called dispersion
  ## -- 1/r as dispersion  in the paper
  ##    "Droplet scRNA-Seq is not zero-inflated." Nature Biotech, 2020
  ## -- r as dispersion in stan
  ## Here we choose r as dispersion following stan.
  r <- nbfit$estimate["size"]
  mu <- nbfit$estimate["mu"]
  p <- r / (r + mu)
  v <- mu + mu^2 / r
  p0 <- (r / (r + mu))^r

  invisible(list(
    size = r, r = r, dispersion_stan = r,
    prob = p, variance = v, nbfit = nbfit,
    p0 = p0
  ))
}

## modification sads::dpoilog
dpoilog <- function(x, mu, sig, log = F, verbose = F) {
  #### FIX: poilog::dpoilog throws an error if an invalid parameter is entered
  #### so we have to circumvent the error here:
  if (length(mu) > 1 | length(sig) > 1) {
    stop("Vectorization of parameters not implemented")
  }
  to.NaN <- NULL
  if (!is.finite(mu)) {
    to.NaN <- seq_len(length(x))
  }
  if (!is.finite(sig) | sig <= 0) {
    to.NaN <- seq_len(length(x))
    if (verbose) {
      message("sig is not finite or less than zero: ", sig)
    }
  }
  mu[!is.finite(mu)] <- 1
  sig[!is.finite(sig) | sig <= 0] <- 1
  y <- poilog::dpoilog(x, mu, sig)
  y[to.NaN] <- NaN
  if (any(is.nan(y))) {
    warning("NaNs produced")
  }
  if (log) {
    return(log(y))
  } else {
    return(y)
  }
}

poilog_negsumlld <- function(mu, sig, x, log = T) {
  -sum(dpoilog(x, mu, sig, log))
}

poilog_mle <- function(x) {
  tmp_lambdas <- log((x + 0.01))
  initpars <- list(mu = mean(tmp_lambdas),
                   sig = sd(tmp_lambdas + 0.1))
  tryCatch({
    myfit <- bbmle::mle2(poilog_negsumlld,
                start = initpars,
                optimizer = "optim",
                data = list(x = x),
                lower = list(mu = -Inf, sig = 0.00001),
                method = "L-BFGS-B",
        hessian = F,
        )
    new("fitsad", myfit, sad = "poilog", distr = "discrete",
      trunc = NaN)
    },
    error = function(cond) {
      message(cond)
      return(NULL)
    }
  )
}

prob_zero_poilognm <- function(x) {
  # sads more stable than poilog
  myfit <- poilog_mle(x)
  if (is.null(myfit)) {
    out <- list(p0 = NA)
  } else {
    mle <- coef(myfit) # fit$par
    mu <- mle["mu"]
    sig <- mle["sig"]
    loglik <- as.numeric(logLik(myfit))
    aic <- AIC(myfit)
    p0 <- sads::dpoilog(0, mu = mu, sig = sig, log = F)

    out <- list(
      fit_poilognm = myfit,
      mu = mu,
      sig = sig,
      loglik = loglik,
      aic = aic,
      p0 = p0
    )
  }
  invisible(out)
}

poislog_negsumlld <- function(mu, sig, cnt, tcnt) {
  if (sig <= 0) {
    stop("std should be larger than zero.")
  }
  mymu <- mu + log(tcnt)
  -sum(vapply(seq_len(length(cnt)),
    FUN = function(i) {
      dpoilog(cnt[i], mymu[i], sig, log = TRUE)
    },
    FUN.VALUE = 1.0
  ))
}


## modification sads::fitpoilog
## we add sequence depth per cell as scale factor
## x is the count for a specific gene;
## s is the total count for the cells
prob_zero_poislognm <- function(x, s, rmoutliers = T, method = "L-BFGS-B") {
  if (length(x) != length(s)) {
    stop("x and s are unequal lengths.")
  }
  if (any(s == 0)) {
    stop("s has zero element.")
  }
  if (rmoutliers) {
    outliers <- is_outlier(x)
    x <- x[!outliers]
    s <- s[!outliers]
  }
  tmp_lambdas <- log((x + 0.01) / s)
  initpars <- list(mu = mean(tmp_lambdas), sig = sd(tmp_lambdas + 0.1))

  myfit <- tryCatch({
      bbmle::mle2(poislog_negsumlld,
        start = initpars,
        optimizer = "optim",
        data = list(cnt = x, tcnt = s),
        lower = c(mu = -Inf, sig = 0.00001),
        method = method,
        control = list(),
        hessian = F
      )
    },
    error = function(cond) {
      message(cond)
      return(NULL)
    }
  )
  if (is.null(myfit)) {
    return(invisible(list(mdnp0 = NA, meanp0 = NA)))
  }

  myfitsad <- new("fitsad", myfit,
    sad = "poilog", distr = "discrete",
    trunc = NaN
  )

  mypars <- coef(myfitsad)
  mu <- mypars["mu"]
  sig <- mypars["sig"]
  loglik <- as.numeric(logLik(myfitsad))
  aic <- AIC(myfitsad)

  mdnp0 <- sads::dpoilog(0, mu = mu + log(median(s)), sig = sig, log = F)
  meanp0 <- sads::dpoilog(0, mu = mu + log(mean(s)), sig = sig, log = F)
  invisible(list(
    fitpoislognm = myfitsad,
    mu_init = initpars["mu"],
    sig_init = initpars["sig"],
    mu = mu,
    sig = sig,
    loglik = loglik,
    aic = aic,
    mdnp0 = mdnp0,
    meanp0 = meanp0
  ))
}


## * summary zero ratio estimations.
estimate_zeroratio <- function(x, s, rmoutliers = T) {
  ## x is the cnt for the gene
  ## s is the scale factor (sequence depth) for the cell
  if (rmoutliers) {
    outliers <- is_outlier(x)
    x <- x[!outliers]
    s <- s[!outliers]
  }
  obs_zr <- empirical_prob_zero(x, F)
  poi_zr <- prob_zero_poi(x, F)$p0
  spoi_zr <- prob_zero_poi_scale(x, s, F)$mdnp0
  poilognm_zr <- prob_zero_poilognm(x)$p0
  poislognm_zr <- prob_zero_poislognm(x, s, F)$mdnp0
  nb_zr <- prob_zero_nb(x, F)$p0
  invisible(list(
    obs = obs_zr,
    poi = poi_zr,
    pois = spoi_zr,
    poilognm = poilognm_zr,
    poislognm = poislognm_zr,
    nb = nb_zr
  ))
}
