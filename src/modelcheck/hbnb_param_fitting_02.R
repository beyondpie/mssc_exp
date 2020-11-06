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
                            seed = 1, numiter = 5000,
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
                           seed = 1, numiter = 5000,
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

  ## Result must be set even when the fitting is bad,
  ## we'll use init strategy to set up the parameters.

  ## y cannot be all the zeros (not test in the code)

  k <- max(vec_of_ind) # num of ind
  result <- list(
    mu0 = NaN, r0 = NaN,
    mu_cond = c(NaN, NaN),
    mu_ind = rep(NaN, k),
    s1 = FALSE,
    s2 = rep(FALSE, 2),
    s3 = rep(FALSE, k)
  )
  ## if (sum(y) < 1) {
  ##   message("[BUG]: sum of y is zero.")
  ##   return(invisible(result))
  ## }
  mur_all <- stanfit_scalenb(s, y, snbm)
  result$mu0 <- mur_all$mu
  result$r0 <- mur_all$r
  result$s1 <- mur_all$success


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
      result$s2[i] <- mui$success
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
      result$s3[i] <- mui$success
    }
  }
  return(invisible(result))
}

fit_mg_snb <- function(cnt, s, cond, ind,
                       snbm, snbm_for_mucond,
                       murnm = c("mu0", "r0"),
                       mucondnm = str_glue_vec("mu_cond", seq_len(2)),
                       muindnm = str_glue_vec("mu_ind", seq_len(k))) {
  ## fit each gene (row of cnt) for mssc
  ## return each gene in row.

  k <- max(ind)
  t_res <- vapply(seq_len(nrow(cnt)),
    FUN = function(i) {
      r <- stanfit_gwsnb_to_ind_level(
        s = s, y = cnt[i, ], vec_of_cond = cond,
        vec_of_ind = ind,
        snbm = snbm,
        snbm_for_mucond = snbm_for_mucond
      )
      return(invisible(c(r$mu0, r$r0, r$mu_cond, r$mu_ind)))
    },
    FUN.VALUE = rep(0.0, 2 + 2 + k)
  )

  colnames(t_res) <- rownames(cnt)
  rownames(t_res) <- c(murnm, mucondnm, muindnm)
  return(invisible(t(t_res)))
}


est_mu <- function(mu, scale = 1.96^2) {
  m <- mean(mu)
  v <- var(mu) * scale
  return(invisible(c(m, v)))
}

est_r <- function(r) {
  m <- mean(r)
  v <- var(r)
  beta <- m / v
  alpha <- beta * m
  return(invisible(c(alpha, beta)))
}

est_mucond <- function(mucond, scale = 1.96^2) {
  ## mucond: g by 2
  d <- vapply(
    seq_len(nrow(mucond)),
    function(r) {
      max(abs(mucond[r, ]))
    },
    0.0
  )
  v <- max(d * scale)
  return(invisible(v))
}

est_muind <- function(muind, scale = 1.96^2) {
  ## muind: g by k
  v <- vapply(
    seq_len(ncol(muind)),
    function(c) {
      var(muind[, c]) * scale
    },
    0.0
  )
  m <- mean(v)
  vv <- var(v)
  beta <- m / vv
  alpha <- beta * m
  return(invisible(list(
    varofind = v,
    hp_varofind = c(alpha, beta)
  )))
}


init_hbnb_params <- function(est_mg_mat,
                             murnm = c("mu0", "r0"),
                             mucondnm = str_glue_vec("mu_cond", seq_len(2)),
                             muindnm = str_glue_vec("mu_ind", seq_len(k)),
                             scale = 1.96^2) {
  ## generate the initial hbnb params
  ## also set the hp params: mu0

  mu <- est_mg_mat[, murnm[1]]
  r <- est_mg_mat[, murnm[2]]
  mu_cond <- est_mg_mat[, mucondnm]
  mu_ind <- est_mg_mat[, muindnm]

  r1 <- est_mu(mu, scale)
  r2 <- est_r(r)
  r3 <- est_mucond(mu_cond, scale)
  r4 <- est_muind(mu_ind, scale)

  hp_params <- list(mu0 = r1[1])

  init_params <- list(
    hp_r = r2,
    nb_r = r,
    varofmu = r1[2],
    mu = mu,
    varofcond = r3,
    mu_cond = mu_cond,
    hp_varofind = r4$hp_varofind,
    varofind = r4$varofind,
    mu_ind = mu_ind
  )
  return(invisible(
    list(
      hp = hp_params,
      init = init_params
    )
  ))
}
