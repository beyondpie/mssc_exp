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
    stringr::str_glue("{nm}{lround}{lround}{i}{rround}")
  },
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
  r <- if (v > m) {
    m^2 / (v - m)
  } else {
    r_default
  }
  invisible(list(
    mu = init_snb_logmu(y, median(s)),
    r = r
  ))
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
    result$mu <- init_mur$mu
    result$r <- init_mur$r
    result$success <- FALSE
  }
  invisible(result)
}

stanfit_snb_fr <- function(s, r, y, model,
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

stanfit_gwsnb_to_ind_level <- function(s, y, vec_of_cond,
                                       vec_of_ind, snbm,
                                       snbm_for_mucond,
                                       r_default = 10.0) {
  ## fit mu0, r0; then mu_cond; then mu_ind
  ## index follows stan hbnb, starting from 1.

  ## In most cases, even the fitting is bad,
  ## we'll use init strategy to set up the parameters.
  ## only when y are all zeros, we will return NaN.
  ## but this also means the data is strange.

  k <- max(vec_of_ind) # num of ind
  result <- list(
    mu0 = NaN, r0 = NaN,
    mu_cond = c(NaN, NaN),
    mu_ind = rep(NaN, k)
  )
  if (sum(y) < 1) {
    message("[BUG]: sum of y is zero.")
    return(invisible(result))
  }
  mur_all <- stanfit_scalenb(s, y, snbm)
  result$mu0 <- mur_all$mu
  result$r0 <- mur_all$r

  ## for mu_cond
  for (i in c(1, 2)) {
    index <- (vec_of_cond == i)
    yi <- y[index]
    if (sum(yi) < 1) {
      message("y are zeros.")
      result$mu_cond[i] <- 0.0
    } else {
      mui <- stanfit_snb_fr(
        s = s[index],
        r = result$r0,
        y = yi,
        model = snbm_for_mucond
      )
      result$mu_cond[i] <- mui$mu - result$mu0
    }
  }

  ## for mu_ind
  for (i in seq_len(k)) {
    index <- (vec_of_ind == i)
    yi <- y[index]
    if (sum(yi) < 1) {
      message("y are zeros")
      result$mu_ind[i] <- 0.0
    } else {
      mui <- stanfit_snb_fr(
        s = s[index],
        r = result$r0,
        y = yi,
        model = snbm_for_mucond
      )
      condi <- vec_of_cond[vec_of_ind == i][1]
      result$mu_ind[i] <- mui$mu - result$mu0 - result$mu_cond[condi]
    }
  }
  return(invisible(result))
}


fit_multgenes_nb <- function(cnt, resp, sumcnt,
                             scale_nb_model,
                             seed = 355113,
                             id_control = 0) {
  result <- lapply(seq_len(nrow(cnt)),
    FUN = function(i) {
      unlist(fit_singlegene_nb(
        i, cnt, resp, sumcnt,
        scale_nb_model, seed, id_control
      ))
    }
  )
  matres <- do.call(rbind, result)
  rownames(matres) <- rownames(cnt)
  invisible(matres)
}
