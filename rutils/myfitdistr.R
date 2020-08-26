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
  ##
  ## -- 1/r as dispersion  in the paper
  ##    "Droplet scRNA-Seq is not zero-inflated." Nature Biotech, 2020
  ## -- r as dispersion in stan
  ## Here we choose r as dispersion following stan.
  r <- nbfit$estimate["size"]
  mu <- nbfit$estimate["mu"]
  p <- r / (r + mu)
  v <- mu + mu^2 / r

  ## when r is extremely large, we can use this as the limitation of r -> inf.
  ## In fact, when r -> inf, NB converges to a Poisson dist, where the lambda
  ## parameter in Poisson is the mean in this NB.
  ## phi <- 1 / r
  ## if (phi == 0.0) {
  ##   p0 <- exp(-mu)
  ## } else {
  ##   p0 <- (r / (r + mu))^r
  ## }

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
  if (length(mu) > 1 | length(sig) > 1) stop("Vectorization of parameters not implemented")
  to.NaN <- NULL
  if (!is.finite(mu)) to.NaN <- 1:length(x)
  if (!is.finite(sig) | sig <= 0) {
    to.NaN <- 1:length(x)
    if (verbose) {
      message("sig is not finite or less than zero: ", sig)
    }
  }
  mu[!is.finite(mu)] <- 1
  sig[!is.finite(sig) | sig <= 0] <- 1
  y <- poilog::dpoilog(x, mu, sig)
  y[to.NaN] <- NaN
  if (any(is.nan(y))) warning("NaNs produced")
  if (log) return(log(y))
  else return(y)
}


## modification sads::fitpoilog
## we add sequence depth per cell as scale factor
## x is the count for a specific gene;
## s is the total count for the cells
prob_zero_poislognm <- function(x, s, rmoutliers = T, method = "BFGS") {
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
  tmp <- poilog::poilogMLE(x, startVals = c(mu = mean(log(x + 0.1)) + log(0.5),
    sig = sd(log(x + 0.1))),
  zTrunc = FALSE, method = "BFGS",
  control = list(maxit = 1000))
  tmp <- tmp$par
  initpars <- list(mu = tmp[1] - log(median(s)), sig = tmp[2])

  negloglikelihood <- function(mu, sig, cnt, tcnt) {
    mymu <- mu + log(tcnt)
    -sum(vapply(seq_len(length(cnt)), FUN = function(i) {
      ## during optimization, sig will sometimes be negative.
      ## so add abs(sig) here, which is different with the sads.
      dpoilog(cnt[i], mymu[i], abs(sig), log = TRUE)},
    FUN.VALUE = 1.0
    ))
  }
  result <- bbmle::mle2(negloglikelihood,
    start = initpars,
    data = list(cnt = x, tcnt = s),
    method = method)

  myfitsad <- new("fitsad", result, sad = "poilog", distr = "discrete",
    trunc = NaN)
  mypars <- coef(myfitsad)
  mu <- mypars["mu"]
  sig <- mypars["sig"]
  loglik <- as.numeric(logLik(myfitsad))
  aic <- AIC(myfitsad)

  mdnp0 <- sads::dpoilog(0, mu = mu + log(median(s)), sig = sig, log = F)
  meanp0 <- sads::dpoilog(0, mu = mu + log(mean(s)), sig = sig, log = F)
  invisible(list(fitpoislognm = myfitsad,
    mu_init = initpars["mu"],
    sig_init = initpars["sig"],
    mu = mu,
    sig = sig,
    loglik = loglik,
    aic = aic,
    mdnp0 = mdnp0,
    meanp0 = meanp0))
}

## *  from quminorm paper with modifications
llcurve_poilog <- function(xmax, lpar, add = TRUE, quadpts = 1000, ...) {
  # Draw the PMF curve on log-log axes for a Poisson-lognormal distribution
  # Curve goes from zero to xmax
  # lpar are the mu,sigma parameters
  # (log scale, default for sads and poilog packages)
  f <- function(t) {
    sads::dpoilog(floor(expm1(t)), mu = lpar[1], sig = lpar[2], log = TRUE)
  }
  curve(f, from = 0, to = log1p(xmax), add = add, ...)
}

llcurve_nb <- function(xmax, lpar, add = TRUE, ...) {
  # Draw the PMF curve on log-log axes for a negative binomial distribution
  # Curve goes from zero to xmax
  # lpar are the size and mu parameters
  f <- function(t) {
    dnbinom(floor(expm1(t)), size = lpar[1], mu = lpar[2], log = TRUE)
  }
  curve(f, from = 0, to = log1p(xmax), add = add, ...)
}

poilog_mle <- function(x, om = "BFGS", ...) {
  # fit<-poilog::poilogMLE(x,startVals=st,zTrunc=FALSE,method=om,...)
  # mle<-fit$par
  # sads more stable than poilog
  fit <- sads::fitpoilog(x, trunc = NULL, method = om, skip.hessian = TRUE, ...)
  mle <- coef(fit) # fit$par
  attr(mle, "loglik") <- as.numeric(logLik(fit)) # fit$logLval
  mle
}

prob_zero_poilognm <- function(x, om = "BFGS", ...) {
  # sads more stable than poilog
  library(sads)
  myfit <- sads::fitpoilog(x, trunc = NULL, method = om,
    skip.hessian = F, ...)
  mle <- coef(myfit) # fit$par
  mu <- mle["mu"]
  sig <- mle["sig"]
  loglik <- as.numeric(logLik(myfit))
  aic <- AIC(myfit)
  p0 <- sads::dpoilog(0, mu = mu, sig = sig, log = F)
  invisible(list(fit_poilognm = myfit,
    mu = mu,
    sig = sig,
    loglik = loglik,
    aic = aic,
    p0 = p0
  ))
}

nb_mle <- function(x, ...) {
  fit <- fitdistrplus::fitdist(x, "nbinom", keepdata = FALSE, ...)
  mle <- coef(fit)
  attr(mle, "loglik") <- logLik(fit)
  mle
}

mle_matrix <- function(m, lik = c("poilog", "nb"), ...) {
  # m a matrix with samples in the columns
  # returns a data frame with nrows=ncols(m)
  # result includes MLEs for each model, log likelihood, and BIC
  lik <- match.arg(lik)
  mle_funcs <- list(poilog = poilog_mle, nb = nb_mle)
  f <- mle_funcs[[lik]]
  mle_func <- function(x) {
    tryCatch({
      mle <- f(x, ...)
      c(mle, loglik = attr(mle, "loglik"))
    },
    error = function(e) {
      rep(NA, 3)
    })
  }
  if (is(m, "sparseMatrix")) {
    apply_func <- function(m) {
      m <- slam::as.simple_triplet_matrix(m)
      res <- slam::colapply_simple_triplet_matrix(m, mle_func)
      as.data.frame(do.call(rbind, res))
    }
  } else {
    apply_func <- function(m) { as.data.frame(t(apply(m, 2, mle_func))) }
  }
  res <- apply_func(m)
  res$bic <- -2 * res$loglik + log(nrow(m)) * 2
  res
}

poilog_mle_matrix <- function(m, ...) {
  # mle_matrix for poisson-lognormal
  # result includes mu,sigma params
  mle_matrix(m, "poilog", ...)
}

nb_mle_matrix <- function(m, ...) {
  mle_matrix(m, "nb", ...)
}

nb_pzero2mu <- function(lpz, size) {
  # Assuming the data follows a negative binomial distribution
  # if we fix the size parameter to a specified value
  # we can infer the mean (mu) parameter from the fraction of zeros
  # lpz is a vector of the log of zero fraction for each cell
  # size can be a scalar or vector
  size * expm1(-lpz / size)
}

poilog_pzero2mu <- function(lpz, sig = 2.5, lims = c(-200, 200)) {
  # Assuming the data follows a Poisson-lognormal
  # if we fix the 'sig' parameter to a specified value
  # we can infer the mu parameter from the fraction of zeros
  # mu != mean of the lognormal, but e^mu is median of lognormal
  # mu can be negative.
  # lpz is a vector of the log of zero fraction for each cell
  # lims are the lower,upper bounds for the mu parameter passed to uniroot
  inner <- function(x, s) {
    # x is an element of lpz
    if (is.na(s)) { return(NA) }
    f <- function(mu) {
      x - sads::dpoilog(0, mu, sig = s, log = TRUE)
    }
    uniroot(f, lims)$root
  }
  if (length(sig) == 1) { # single sig parameter for all cells
    return(vapply(lpz, inner, FUN.VALUE = 1.0, s = sig))
  } else { # different tail parameter for each cell
    return(mapply(inner, lpz, sig))
  }
}
