## use stan for snb fitting
## provide other functions specific for hbnb

## ** default parameters
hpgamma_default <- c(1.0, 1.0)
## default invgamma param
hpinvg_default <- c(1.0, 10.0)
r_default <- 20.0
mu_default <- 0.0
varofmu_default <- 25.0
varofcond_default <- 4.0

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
    stringr::str_glue("{nm}{lround}{i}{rround}")
  },
  FUN.VALUE = ""
  ))
}

str_glue_mat <- function(nr, nc, nm = "MuInd",
                         lround = "[", rround = "]") {
  ## row-wise names
  r <- rep("", nr * nc)
  for (i in 1:nr) {
    for (j in 1:nc) {
      r[nc * (i - 1) + j] <- stringr::str_glue("{nm}{lround}{i},{j}{rround}")
    }
  }
  return(invisible(r))
}


init_snb_logmu <- function(y, median_of_s) {
  ## init scale negative binomial log mu.
  ## y should not be all zeros.
  ## otherwise, will return -20
  if (sum(y) < 1) {
    message("[WARNING]: y are all zeros.")
    return(-20)
  }
  return(log(mean(y)) - log(median_of_s))
}

init_snb <- function(s, y, init_r = r_default) {
  ## directly use sample mean and variance to
  ## estiamte the parameters for scaled negative binomial
  ## ref: MASS::fitdistr for nb
  m <- mean(y)
  v <- var(y)
  r <- if (v > m) {
    m^2 / (v - m)
  } else {
    init_r
  }
  invisible(list(
    mu = init_snb_logmu(y, median(s)),
    r = r
  ))
}


stanfit_scalenb <- function(s, y, scale_nb_model,
                            seed = 1, numiter = 5000,
                            refresh = 5000, init_r = r_default,
                            too_big_r = r_default) {
  ## mu in scaled log level, and minus log(s)
  # fit scaled negative binomial using stan
  result <- list(mu = NaN, r = NaN, success = FALSE)
  init_mur <- init_snb(s, y, init_r)

  capture.output(opt <- scale_nb_model$optimize(
    data = list(n = length(s), s = s, y = y),
    seed = seed,
    refresh = refresh,
    iter = numiter,
    init = list(init_mur),
    algorithm = "lbfgs"
  ))

  if (is_vi_or_opt_success(opt)) {
    t <- opt$mle()
    r <- t["r"]
    if (is.infinite(t["mu"])) {
      ## This should not happen.
      message("[WARNING] Fitting faces mu inifinity, but opt success.")
      result$mu <- init_mur$mu
      result$r <- init_mur$r
      result$success <- FALSE
    } else if (r > too_big_r) {
      message(stringr::str_glue("[WARNING]: Fitted r {r} > {too_big_r}"))
      message(stringr::str_glue("Set r as {r_default}"))
      result$mu <- init_mur$mu
      result$r <- init_mur$r
      result$success <- FALSE
    } else {
      result$mu <- t["mu"]
      result$r <- r
      result$success <- TRUE
    }
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
  ## fitting scaled negative binomial while fixing r

  result <- list(mu = NaN, success = FALSE)
  # fit scaled negative binomial using stan
  n <- length(s)
  ## ref: MASS::fitdistr for nb
  ## y should be not all zeros.
  init_mu <- init_snb_logmu(y, median(s))
  capture.output(opt <- model$optimize(
    data = list(n = n, s = s, y = y, r = r),
    seed = seed,
    refresh = refresh,
    iter = numiter,
    init = list(list(
      mu = init_mu
    )),
    algorithm = "lbfgs"
  ))
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
                                       to_ind = TRUE) {
  ## fit mu0, r0; then mu_cond; then mu_ind
  ## index follows stan hbnb, starting from 1.

  ## Result must be set even when the fitting is bad,
  ## we'll use init strategy to set up the parameters.

  ## y cannot be all the zeros (not test in the code)
  k <- max(vec_of_ind) # num of ind
  j <- max(vec_of_cond) # num of cond, should be 2
  result <- list(
    mu0 = mu_default, r0 = r_default,
    mu_cond = rep(0.0, j),
    mu_ind = rep(0.0, k),
    s1 = FALSE,
    s2 = rep(FALSE, 2),
    s3 = rep(FALSE, k)
  )
  if (sum(y) < 1) {
    message("[WARNING]: [BAD EXPRESSION] sum of y is zero.")
    return(invisible(result))
  }
  mur_all <- stanfit_scalenb(s, y, snbm)
  result$mu0 <- mur_all$mu
  result$r0 <- mur_all$r
  result$s1 <- mur_all$success


  ## for mu_cond
  for (i in seq_len(j)) {
    index <- (vec_of_cond == i)
    yi <- y[index]
    if (sum(yi) < 1) {
      message(stringr::str_glue("[WARNING] Cond [{i}]: y are zeros."))
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
  if (to_ind) {
    for (i in seq_len(k)) {
      index <- (vec_of_ind == i)
      yi <- y[index]
      if (sum(yi) < 1) {
        message(stringr::str_glue("[WARNING] Ind [{i}]: y are zeros."))
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
  }
  return(invisible(result))
}

est_mu <- function(mu, scale = 1.96^2) {
  ## use median instead of mu
  ## - if mu follows normal, median and mean should be almost same
  ## - when mu has Inf element, median is robust.
  m <- median(mu)
  v <- var(mu) * scale
  if (is.nan(v)) {
    ## happens when mu has Inf element.
    ## mu should not have this element
    message(stringr::str_glue("[WARNING]: [Bad Mu] v is NaN; set as {default_v}."))
    v <- varofmu_default
  }
  return(invisible(c(m, v)))
}

est_r <- function(r) {
  ## TODO: check if we need add threshold for alpha, beta
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

  ## ## each individual share the same variance
  ## ## for different genes.
  ## v shape:  k by 1
  ## v <- vapply(
  ##   seq_len(ncol(muind)),
  ##   function(c) {
  ##     var(muind[, c]) * scale
  ##   },
  ##   0.0
  ## )

  ## each gene share the same variance
  ## v shape: g by 1
  v <- vapply(seq_len(nrow(muind)), function(r) {
    var(muind[r, ]) * scale
  }, FUN.VALUE = 0.0)
  m <- mean(v)
  vv <- var(v)
  beta <- m / vv
  alpha <- beta * m

  return(invisible(list(
    varofind = v,
    hp_varofind = c(alpha, beta)
  )))
}


fit_mg_snb <- function(cnt, s, cond, ind,
                       snbm, snbm_for_mucond,
                       murnm = c("mu0", "r0"),
                       mucondnm = str_glue_vec("mu_cond", seq_len(2)),
                       muindnm = str_glue_vec("mu_ind", seq_len(k)),
                       to_ind = TRUE) {
  ## fit each gene (row of cnt) for mssc
  ## return each gene in row.

  k <- max(ind)
  t_res <- vapply(seq_len(nrow(cnt)),
    FUN = function(i) {
      r <- stanfit_gwsnb_to_ind_level(
        s = s, y = cnt[i, ], vec_of_cond = cond,
        vec_of_ind = ind,
        snbm = snbm,
        snbm_for_mucond = snbm_for_mucond,
        to_ind = to_ind
      )
      return(invisible(c(r$mu0, r$r0, r$mu_cond, r$mu_ind)))
    },
    FUN.VALUE = rep(0.0, 2 + 2 + k)
  )

  colnames(t_res) <- rownames(cnt)
  rownames(t_res) <- c(murnm, mucondnm, muindnm)
  return(invisible(t(t_res)))
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
  mu0 <- r1[1]
  varofmu <- r1[2]
  raw_mu <- (mu - mu0) / sqrt(varofmu)

  r2 <- est_r(r)

  varofcond <- est_mucond(mu_cond, scale)
  raw_mu_cond <- mu_cond / sqrt(varofcond)

  ## TODO: set default values if not fitting
  ## to individual level
  r3 <- est_muind(mu_ind, scale)
  ## length equals to gene number
  varofind <- r3$varofind
  ## each row (a gene) divided by the correspond element from varofind
  raw_mu_ind <- mu_ind / sqrt(varofind)
  ## update raw_mu_ind: each gene, share the same varofind
  ## for (i in nrow(mu_ind)) {
  ##   raw_mu_ind[i, ] <- mu_ind[i, ] / sqrt(varofind[i])
  ## }

  hp_params <- list(mu0 = r1[1])

  init_params <- list(
    hp_r = r2,
    nb_r = r,
    varofmu = varofmu,
    mu = mu,
    raw_mu = raw_mu,
    varofcond = varofcond,
    mu_cond = mu_cond,
    raw_mu_cond = raw_mu_cond,
    hp_varofind = r3$hp_varofind,
    varofind = varofind,
    mu_ind = mu_ind,
    raw_mu_ind = raw_mu_ind
  )
  return(invisible(
    list(
      hp = hp_params,
      init = init_params
    )
  ))
}
