is_outlier <- function(x, up_prob = 0.995) {
  invisible(x > 5 * quantile(x, up_prob))
}

is_vi_or_opt_success <- function(vi_or_opt) {
  ## https://github.com/stan-dev/cmdstanr/issues/332
  if (vi_or_opt$runset$procs$get_proc(1)$get_exit_status() == 0) {
    return(TRUE)
  }
  return(FALSE)
}

str_glue_vec <- function(nm = "MuInd", a = seq_len(10),
                         lround = "[", rround = "]") {
  invisible(vapply(a, function(i) {
    stringr::str_glue("{nm}{lround}{lround}{i}{rround}")},
    FUN.VALUE = ""
    ))
}

grpoutliers <- function(cntgbc) {
  outliers <- sapply(cntgbc, 1, is_outlier)
  invisible(rowSums(outliers) > 0)
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

  nbfit <- tryCatch(
    {
      MASS::fitdistr(x,
        densfun = "negative binomial",
        lower = c(0.00001, 0.00001)
      )
    },
    error = function(cond) {
      message(cond)
      return(NULL)
    }
  )

  if (is.null(nbfit)) {
    return(invisible(list(p0 = NA)))
  }
  ## r or its reciprocal 1/r is also called dispersion
  ## -- 1/r as dispersion  in the paper
  ##    "Droplet scRNA-Seq is not zero-inflated." Nature Biotech, 2020
  ## -- r as dispersion in stan
  ## Here we choose r as dispersion following stan.
  r <- unname(nbfit$estimate["size"])
  mu <- unname(nbfit$estimate["mu"])
  p <- unname(r / (r + mu))
  v <- unname(mu + mu^2 / r)
  p0 <- unname((r / (r + mu))^r)

  invisible(list(
    size = r, r = r, dispersion_stan = r,
    prob = p, variance = v, nbfit = nbfit,
    p0 = p0
  ))
}


stanfit_scalenb <- function(s, y, scale_nb_model,
                            seed = 355113, numiter = 5000,
                            refresh = 5000, r_default = 10) {
  ## mu in scaled log level, and minus log(s)
  ## use stan to fit
  result <- list(mu = NaN, r = NaN)
  # fit scaled negative binomial using stan
  n <- length(s)
  ## ref: MASS::fitdistr for nb
  m <- mean(y)
  v <- var(y)
  r <- if (v > m) {
    m^2 / (v - m)
  } else {
    r_default
  }
  opt <- scale_nb_model$optimize(
    data = list(n = n, s = s, y = y),
    seed = seed,
    refresh = refresh,
    iter = numiter,
    init = list(list(
      mu = log(m) - log(median(s)),
      r = r
    )),
    algorithm = "lbfgs"
  )
  ## https://github.com/stan-dev/cmdstanr/issues/332
  if (opt$runset$procs$get_proc(1)$get_exit_status()  == 0) {
    t <- opt$mle()
    result$mu <- t["mu"]
    result$r <- t["r"]
  }
  invisible(result)
}

stanfit_scalenb_fixed_r <- function(s, r, y, scale_nb_fixed_r_model,
                                    seed = 355113, numiter = 5000,
                                    refresh = 5000) {
  ## mu in scaled log level, and minus log(s)
  ## use stan to fit
  result <- list(mu = NaN)
  # fit scaled negative binomial using stan
  n <- length(s)
  ## ref: MASS::fitdistr for nb
  m <- mean(y)
  opt <- scale_nb_fixed_r_model$optimize(
    data = list(n = n, s = s, y = y, r = r),
    seed = seed,
    refresh = refresh,
    iter = numiter,
    init = list(list(
      mu = log(m) - log(median(s))
    )),
    algorithm = "lbfgs"
  )
  ## https://github.com/stan-dev/cmdstanr/issues/332
  if (opt$runset$procs$get_proc(1)$get_exit_status() == 0) {
    t <- opt$mle()
    result$mu <- t["mu"]
  }
  invisible(result)
}

stanfit_snb_for_muind <- function(s, y, vec_of_ind, mu, r,
                                  model, nm = "MuInd",
                                  seed = 1L, numiter= 5000,
                                  refresh = 5000) {
  ## fit muind under scale negative binomial scene
  ## vec_of_ind:  [1,1, 2, 2, 3,4] like.

  k <- max(vec_of_ind)
  n <- length(s)
  param_nms <- str_glue_vec(nm = nm, a = seq_len(k))
  mu_ind <- rep(NaN, k)
  names(mu_ind) <- param_nms
  opt <- model$optimize(
                 data = list(n = n, s = s, y= y, k =k,
                             ind = vec_of_ind, mu = mu,
                             r = r),
                 init = list(list(mu = rep(0.0, k))),
                 algorithm = "lbfgs"
               )
  if (is_vi_or_opt_success(opt)) {
    t <- opt$mle()
    mu_ind <- t[param_nms]
  }

  invisible(mu_ind)
}


stanfit_gwsnb_till_cond_level <- function(s, y,
                                          s_control,
                                          y_control, s_case, y_case,
                                          scale_nb_model,
                                          scale_nb_fixed_r_model) {
  ## fit negative binomial, and get the mean and dispersion.
  ## negative binommial scaled (s) and mean in log level (nb v2 as in STAN).
  ## only global and conditional level, no individuals.
  ## no outlier detection.

  ## a simple version for scaled negative binomial estimation.
  ## i.e, fit the nb, then divided by the scale of s.

  result <- list(
    mu0 = NaN, r0 = NaN,
    mu_cond = c(NaN, NaN)
  )

  mur_all <- stanfit_scalenb(s, y, scale_nb_model)
  if (is.nan(mur_all$mu)) {
    message("Cannot fit negative binomial for the overall mean.")
    return(invisible(result))
  }

  result$mu0 <- mur_all$mu
  result$r0 <- mur_all$r

  mur_control <- stanfit_scalenb_fixed_r(
    s_control, mur_all$r, y_control,
    scale_nb_fixed_r_model
  )
  if (!is.nan(mur_control$mu)) {
    result$mu_cond[1] <- mur_control$mu - mur_all$mu
  }

  mur_case <- stanfit_scalenb_fixed_r(
    s_case, mur_all$r, y_case,
    scale_nb_fixed_r_model
  )
  if (!is.nan(mur_case$mu)) {
    result$mu_cond[2] <- mur_case$mu - mur_all$mu
  }
  return(invisible(result))
}


nbfit_mur <- function(x) {
  ## r or its reciprocal 1/r is also called dispersion
  ## -- 1/r as dispersion  in the paper
  ##    "Droplet scRNA-Seq is not zero-inflated." Nature Biotech, 2020
  ## -- r as dispersion in stan
  ## Here we choose r as dispersion following stan.
  fit <- tryCatch(
    {
      MASS::fitdistr(x,
        densfun = "negative binomial",
        lower = c(0.00001, 0.00001)
      )
    },
    error = function(cond) {
      message(cond)
      return(NULL)
    }
  )
  if (is.null(fit)) {
    return(invisible(list(mu = NaN, r = NaN)))
  }
  return(invisible(list(
    mu = fit$estimate["mu"],
    r = fit$estimate["size"]
  )))
}

fit_gwnb_s2_till_cond_level <- function(y, y_control, y_case,
                                        scale_of_s) {
  ## fit negative binomial, and get the mean and dispersion.
  ## negative binommial scaled (s) and mean in log level (nb v2 as in STAN).
  ## only global and conditional level, no individuals.
  ## no outlier detection.

  ## a simple version for scaled negative binomial estimation.
  ## i.e, fit the nb, then divided by the scale of s.

  result <- list(
    mu0 = NaN, r0 = NaN,
    mu_cond = c(NaN, NaN),
    r_cond = c(NaN, NaN)
  )

  mur_all <- nbfit_mur(y)
  if (is.nan(mur_all$mu)) {
    message("Cannot fit negative binomial for the overall mean.")
    return(invisible(result))
  }

  result$mu0 <- log(mur_all$mu / scale_of_s)
  result$r0 <- mur_all$r

  mur_control <- nbfit_mur(y_control)
  if (!is.nan(mur_control$mu)) {
    ## both of mus divided by scale_of_s inner log, then that part
    ## is cancelled.
    result$mu_cond[1] <- log(mur_control$mu) - log(mur_all$mu)
    result$r_cond[1] <- mur_control$r
  }

  mur_case <- nbfit_mur(y_case)
  if (!is.nan(mur_case$mu)) {
    ## both of mus divided by scale_of_s inner log, then that part
    ## is cancelled.
    result$mu_cond[2] <- log(mur_case$mu) - log(mur_all$mu)
    result$r_cond[2] <- mur_case$r
  }
  return(invisible(result))
}

fit_gwnb_s2_ind_mu <- function(y, vec_of_ind_under_cond,
                               mu0, mu_cond, scale_of_s) {
  ## given fitted mu0 and the mu_cond(either case or control)
  ## we estimate the mu_ind by minus mu0 and mu_cond

  nms <- unique(sort(vec_of_ind_under_cond))
  a <- vapply(
    nms, function(nm) {
      mur <- nbfit_mur(y[vec_of_ind_under_cond == nm])
      if (is.nan(mur$mu)) {
        return(c(NaN, NaN))
      }
      return(c(log(nbfit_mur$mu / scale_of_s) - mu0 - mu_cond))
    },
    c(0.0, 0.0)
  )
  colnames(a) <- nms
  rownames(a) <- c("mu", "r")
  invisible(a)
}

dpoilog <- function(x, mu, sig, log = F, verbose = F) {
  ## modification sads::dpoilog
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
  initpars <- list(
    mu = mean(tmp_lambdas),
    sig = sd(tmp_lambdas) + 0.00001
  )
  tryCatch(
    {
      myfit <- bbmle::mle2(poilog_negsumlld,
        start = initpars,
        optimizer = "optim",
        data = list(x = x),
        lower = list(mu = -Inf, sig = 0.00001),
        method = "L-BFGS-B",
        hessian = F,
      )
      new("fitsad", myfit,
        sad = "poilog", distr = "discrete",
        trunc = NaN
      )
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


prob_zero_poislognm <- function(x, s, rmoutliers = T, method = "L-BFGS-B") {
  ## modification sads::fitpoilog
  ## we add sequence depth per cell as scale factor
  ## x is the count for a specific gene;
  ## s is the total count for the cells
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
  initpars <- list(mu = mean(tmp_lambdas), sig = sd(tmp_lambdas) + 0.00001)

  myfit <- tryCatch(
    {
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

estimate_zeroratios <- function(cntgbc, cellmeta_inds,
                                cellmeta_clusters,
                                genes, whichind = NULL, whichcluster = c(2),
                                rmoutliers = T) {
  if (is.null(whichind)) {
    thecells <- cellmeta_clusters %in% whichcluster
  } else {
    thecells <- (cellmeta_clusters %in% whichcluster) &
      (grepl(whichind, cellmeta_inds))
  }

  totcnt <- Matrix::colSums(cntgbc)[thecells]
  zrs <- lapply(genes, FUN = function(g) {
    estimate_zeroratio(
      cntgbc[g, thecells],
      totcnt, rmoutliers
    )
  }) %>% do.call(what = rbind, args = .)
  invisible(zrs)
}



fit_multgenes_nb <- function(cnt, resp, sumcnt,
                             scale_nb_model,
                             seed = 355113,
                             id_control = 0) {
  result <- lapply(seq_len(nrow(cnt)),
    FUN = function(i) {
      unlist(fit_singlegene_nb(i, cnt, resp, sumcnt,
        scale_nb_model, seed, id_control))
    })
  matres <- do.call(rbind, result)
  rownames(matres) <- rownames(cnt)
  invisible(matres)
}
