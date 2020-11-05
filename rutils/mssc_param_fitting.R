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


init_snb_logmu <- function(y, median_of_s) {
  ## init scale negative binomial log mu.
  ## y should not be all zeros.
  return(log(mean(y)) - log(median_of_s))
}
init_snb <- function(s, y, r_default = 10) {
  ## directly use sample mean and variance to
  ## estiamte the parameters for scaled negative binomial
  ## ref: MASS::fitdistr for nb
  m <- mean(y)
  v <- var(y)
  r <- if(v>m) {
         m^2 / (v-m)
       } else {
         r_default
       }
  invisible(list(mu = init_snb_logmu(y, median(s)),
                 r = r))
}


stanfit_scalenb <- function(s, y, scale_nb_model,
                            seed = 355113, numiter = 5000,
                            refresh = 5000, r_default = 10) {
  ## mu in scaled log level, and minus log(s)
  ## use stan to fit
  result <- list(mu = NaN, r = NaN, success = FALSE)
  # fit scaled negative binomial using stan
  init_mur <- init_snb(s, y, r_default)

  opt <- scale_nb_model$optimize(
    data = list(n = length(s), s = s, y = y),
    seed = seed,
    refresh = refresh,
    iter = numiter,
    init = list(init_mur),
    algorithm = "lbfgs"
    )

  if (is_vi_or_opt_success(opt)) {
    t <- opt$mle()
    result$mu <- t["mu"]
    result$r <- t["r"]
    result$success <- TRUE
  } else {
    ## message("STAN cannot fit scaled nb.")
    result$mu <- init_mur$mu
    result$r <- init_mur$r
    result$success <- FALSE
  }
  invisible(result)
}

stanfit_scalenb_fixed_r <- function(s, r, y, model,
                                    seed = 355113, numiter = 5000,
                                    refresh = 5000) {
  ## mu in scaled log level, and minus log(s)
  ## use stan to fit

  result <- list(mu = NaN, success = FALSE)
  # fit scaled negative binomial using stan
  n <- length(s)
  ## ref: MASS::fitdistr for nb
  ## y should be not all zeros.
  init_mu <- init_snb_logmu(y, median(s))
  opt <- model$optimize(
    data = list(n = n, s = s, y = y, r = r),
    seed = seed,
    refresh = refresh,
    iter = numiter,
    init = list(list(
      mu = init_mu
    )),
    algorithm = "lbfgs"
  )
  if (is_vi_or_opt_success(opt)) {
    t <- opt$mle()
    result$mu <- t["mu"]
    result$success <- TRUE
  } else {
    result$mu <- init_mu
  }
  invisible(result)
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
  ## TODO: stanfit logic is changed.
  if (!is.nan(mur_control$mu)) {
    result$mu_cond[1] <- mur_control$mu - mur_all$mu
  }

  mur_case <- stanfit_scalenb_fixed_r(
    s_case, mur_all$r, y_case,
    scale_nb_fixed_r_model
  )
  ## TODO: stanfit logic is changed.
  if (!is.nan(mur_case$mu)) {
    result$mu_cond[2] <- mur_case$mu - mur_all$mu
  }
  return(invisible(result))
}


stanfit_snb_for_muind <- function(s, y, vec_of_ind, mu, r,
                                  model, nm = "MuInd",
                                  seed = 1L, numiter = 5000,
                                  refresh = 5000) {
  ## fit muind under scale negative binomial scene
  ## vec_of_ind:  [1,1, 2, 2, 3,4] like.

  k <- max(vec_of_ind)
  n <- length(s)
  param_nms <- str_glue_vec(nm = nm, a = seq_len(k))
  mu_ind <- rep(NaN, k)
  names(mu_ind) <- param_nms
  opt <- model$optimize(
    data = list(n = n, s = s, y = y, k = k,
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
