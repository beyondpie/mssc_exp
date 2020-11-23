## hbnb model

## param names in model
## 1) names defined in  stan script (by us)
## 2) names returned by stan sampler/opt
##    (stan adds [i,j] or [i] for matrix/vec)
## 3) names from the estimation directly from the data (only part of params)
##    see fig_mg_sng for details.

## hbnb row-major order considers the situation:
## [g1,1] [g1,2], ..., [g1,second_dim], [g2,1], ...
## - this is hbnb param nms order
## - cmdstan will named the varibles column-major order
##   like: [g1,1], [g2,1], [g3, 1], ..., [g1, 2], [g2, 2], ...
## As long as we use variable names to extract the parameters,
## it's fine.

## * load R env
suppressPackageStartupMessages(library(tidyverse))
library(MCMCpack)
library(cmdstanr)
library(grid)
library(gtable)
library(gridExtra)
library(bayesplot)
library(posterior)
library(bbmle)
library(sads)
library(truncnorm)

## warnings/errors traceback settings
options(error = traceback)
options(warn = 1)
options(mc.cores = 3)


options("import.path" = here::here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## * load stan model
nm_params <- c(
  "nb_r", "hp_r", "varofmu", "mu",
  "varofcond", "mu_cond",
  "hp_varofind", "varofind", "mu_ind",
  "raw_mu", "raw_mu_cond", "raw_mu_ind"
)

snbm_for_mur <- cmdstanr::cmdstan_model(
  here::here("stanutils", "scale_nb.stan"),
  compile = T, quiet = F
)

snbm_for_mucond <- cmdstanr::cmdstan_model(
  here::here("stanutils", "scale_nb_fixed_r.stan"),
  compile = T, quiet = F
)

hbnbm <- cmdstanr::cmdstan_model(
  here::here("src", "dirty_stan", "mssc_hbnb_indeff_across_gene.stan"),
  compile = T, quiet = F
)

## * constants
## ** stan training
num_iter <- 50000
## set refresh as zero to reduce the output log from stan
## https://github.com/stan-dev/cmdstanr/issues/341
opt_refresh <- 0
vi_refresh <- 2000
eval_elbo <- 100
## vi_algorithm <- "fullrank"
vi_algorithm <- "meanfield"
output_samples <- 3000
tol_rel_obj <- 0.001
adapt_iter <- 1000
eta <- 0.5

## ** default parameters
hpgamma_default <- c(1.0, 1.0)
## default invgamma param
hpinvg_default <- c(1.0, 1.0)
## how to set r_default appropriately
r_default <- 50.0
mu_default <- 0.0
varofmu_default <- 16.0
varofmu_default_min <- 2.0
varofcond_default <- 4.0
varofcond_default_min <- 1.0
varofind_default <- 1.0
sd_init_mucond <- 0.01
sd_init_muind <- 0.01

## * functions

is_vi_or_opt_success <- function(vi_or_opt) {
  ## https://github.com/stan-dev/cmdstanr/issues/332
  ## https://mc-stan.org/cmdstanr/reference/fit-method-return_codes.html
  ## if (vi_or_opt$runset$procs$get_proc(1)$get_exit_status() == 0) {
  if (vi_or_opt$return_codes() == 0) {
    return(TRUE)
  }
  return(FALSE)
}

str_glue_vec <- function(nm = "MuInd", n = 10,
                         lround = "[", rround = "]") {
  invisible(vapply(seq_len(n), function(i) {
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
  ## otherwise, will return mu_default
  if (sum(y) < 1) {
    warning("[INIT SNB LOGMU]: all the y are zeros.")
    return(mu_default)
  }
  return(log(mean(y)) - log(median_of_s))
}

init_snb_r <- function(y) {
  ## not check if y are all zeros.
  if (sum(y) < 1) {
    warning("[INIT SNB R]: all the y are zeros.")
    return(r_default)
  }
  m <- mean(y)
  v <- var(y)
  r <- if (v > m) {
    m^2 / (v - m)
  } else {
    r_default
  }
  return(invisible(r))
}

init_snb <- function(s, y) {
  ## directly use sample mean and variance to
  ## estiamte the parameters for scaled negative binomial
  ## ref: MASS::fitdistr for nb
  if (sum(y) < 1) {
    warning("[INIT SNB]: all the y are zeros.")
    return(invisible(list(mu = mu_default, r = r_default)))
  }
  invisible(list(
    mu = init_snb_logmu(y, median(s)),
    r = init_snb_r(y)
  ))
}

stanfit_scalenb <- function(s, y, scale_nb_model,
                            seed = 1, numiter = 5000,
                            refresh = 5000,
                            too_big_r = r_default) {
  ## mu in scaled log level, and minus log(s)
  # fit scaled negative binomial using stan
  result <- list(mu = mu_default, r = r_default, success = FALSE)
  if (sum(y) < 1) {
    warning("[STANFIT SCALENB]: all the y are zeros.")
    return(invisible(result))
  }

  init_mur <- init_snb(s, y)

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
    logmu <- t["mu"]
    r <- t["r"]
    if (is.infinite(logmu)) {
      ## This should not happen.
      warning("[STANFIT SNB USING OPT] Mu is inifinity. Using init params.")
      result$mu <- init_mur$mu
      result$r <- init_mur$r
      result$success <- FALSE
    } else if (r > too_big_r) {
      warning(stringr::str_glue(
        "[STANFIT SNB USING OPT]: r {r} > {too_big_r}. ",
        "Using init params"
      ))
      result$mu <- init_mur$mu
      result$r <- init_mur$r
      result$success <- FALSE
    } else {
      result$mu <- logmu
      result$r <- r
      result$success <- TRUE
    }
  } else {
    warning("[STANFIT SNB USING OPT]: Failed. Using init params.")
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

  result <- list(mu = mu_default, success = FALSE)
  if (sum(y) < 1) {
    warning(paste(
      "[STANFIT SNB FIXED R USING OPT]:",
      "all the y are zeros. Using default"
    ))
    return(invisible(result))
  }
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
    warning(paste(
      "[STANFIT SNB FIXED R USING OPT]:",
      "OPT failed. Using init_snb_logmu result."
    ))
    result$mu <- init_mu
  }
  invisible(result)
}

stanfit_gwsnb_to_cond_level <- function(s, y, vec_of_cond,
                                        vec_of_ind, snbm,
                                        snbm_for_mucond) {
  ## fit mu0, r0; then mu_cond
  ## index follows stan hbnb, starting from 1.

  ## Result must be set even when the fitting is bad,
  ## we'll use init strategy to set up the parameters.

  ## [MUIND]: all the muind will be filled with
  ## small values around 0.0

  ## y cannot be all the zeros (not test in the code)
  k <- max(vec_of_ind) # num of ind
  j <- max(vec_of_cond) # num of cond, should be 2
  result <- list(
    mu0 = mu_default, r0 = r_default,
    mu_cond = rnorm(n = j, mean = 0, sd = sd_init_mucond),
    mu_ind = rnorm(n = k, mean = 0, sd = sd_init_muind),
    s1 = FALSE,
    s2 = rep(FALSE, 2)
  )
  if (sum(y) < 1) {
    warning("[STANFIT COND LEVEL]: [FOR A GENE] all obs are zeros.")
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
      warning(stringr::str_glue("[STANFIT COND LEVEL] ",
                                "Cond[{i}]: y are zeros."))
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
    warning("[EST MU0 VAR]: VAR is NaN. Using default")
    v <- varofmu_default
  }
  if (v < varofmu_default_min) {
    warning(stringr::str_glue(
      "[EST Mu0 VAR]: VAR {v} < {varofmu_default_min}",
      " Using the default min."
    ))
    v <- varofmu_default_min
  }
  return(invisible(c(m, v)))
}

est_varofcond <- function(mucond, scale = 1.96^2) {
  ## mucond: g by 2
  d <- vapply(
    seq_len(nrow(mucond)),
    function(r) {
      max(abs(mucond[r, ]))
    },
    0.0
  )
  v <- max(d * scale)
  if (v < varofcond_default_min) {
    warning(stringr::str_glue(
      "[EST MUCOND VAR]: VAR {v} < {varofcond_default_min}",
      "Using the default min."
    ))
    v <- varofcond_default_min
  }
  return(invisible(v))
}

fit_mg_snb <- function(cnt, s, cond, ind,
                       snbm, snbm_for_mucond,
                       murnm,
                       mucondnm,
                       muindnm) {
  ## fit scaled negative binomial dist for each gene (row of cnt) for mssc
  ## return each gene in row.
  ## individual effect will be set as zeros.

  k <- max(ind)
  t_res <- vapply(seq_len(nrow(cnt)),
    FUN = function(i) {
      r <- stanfit_gwsnb_to_cond_level(
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


init_hbnb_params <- function(est_mg_mat,
                             murnm,
                             mucondnm,
                             muindnm,
                             scale = 1.96^2) {
  ## generate the initial hbnb params
  ## also set the hp params: mu0

  mu <- est_mg_mat[, murnm[1]]
  r <- est_mg_mat[, murnm[2]]
  mu_cond <- est_mg_mat[, mucondnm]
  ## all zeros
  mu_ind <- est_mg_mat[, muindnm]

  r1 <- est_mu(mu, scale)
  mu0 <- r1[1]
  varofmu <- r1[2]
  raw_mu <- (mu - mu0) / sqrt(varofmu)

  varofcond <- est_varofcond(mu_cond, scale)
  raw_mu_cond <- mu_cond / sqrt(varofcond)

  ## varofind: k by 1
  varofind <- rep(varofind_default, ncol(mu_ind))
  ## each row (a gene) divided by the correspond element from varofind
  ## here: mu_ind are all zeros by default.
  raw_mu_ind <- mu_ind / sqrt(varofind)

  hp_params <- list(mu0 = r1[1])

  init_params <- list(
    hp_r = hpgamma_default,
    nb_r = r,
    varofmu = varofmu,
    mu = mu,
    raw_mu = raw_mu,
    varofcond = varofcond,
    mu_cond = mu_cond,
    raw_mu_cond = raw_mu_cond,
    hp_varofind = hpinvg_default,
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

get_hbnb_param_nms <- function(k, j, g) {
  hbnb_param_nms <- list(
    hp_r = str_glue_vec(nm = "hp_r", n = 2),
    nb_r = str_glue_vec(nm = "nb_r", n = g),
    varofmu = "varofmu",
    raw_mu = str_glue_vec(nm = "raw_mu", n = g),
    mu = str_glue_vec(nm = "mu", n = g),
    varofcond = "varofcond",
    raw_mu_cond = str_glue_mat(nm = "raw_mu_cond", nr = g, nc = j),
    mu_cond = str_glue_mat(nm = "mu_cond", nr = g, nc = j),
    raw_mu_ind = str_glue_mat(nm = "raw_mu_ind", nr = g, nc = k),
    mu_ind = str_glue_mat(nm = "mu_ind", nr = g, nc = k),
    hp_varofind = str_glue_vec(nm = "hp_varofind", n = 2),
    ## varofind:  in total k
    varofind = str_glue_vec(nm = "varofind", n = k)
  )
  return(invisible(hbnb_param_nms))
}

get_default_hi_params <- function(k, j, g) {
  ## hi: hyper and initial
  ## default hyper params
  dhp <- list(
    mu0 = rep(mu_default, g),
    hp_varofmu = hpinvg_default,
    hp_alpha_r = hpgamma_default,
    hp_beta_r = hpgamma_default,
    hp_alpha_varofind = hpgamma_default,
    hp_beta_varofind = hpgamma_default,
    hp_varofcond = hpinvg_default
  )
  ## default init params
  dip <- list(
    hp_r = hpgamma_default,
    nb_r = rep(r_default, g),
    varofmu = varofmu_default,
    mu = rep(mu_default, g),
    raw_mu = rep(0.0, g),
    mu_cond = array(0.0, dim = c(g, j)),
    raw_mu_cond = array(0.0, dim = c(g, j)),
    varofcond = varofcond_default,
    mu_ind = array(0.0, dim = c(g, k)),
    raw_mu_ind = array(0.0, dim = c(g, k)),
    ## change dim from g to k, since for each individual
    ## all the genes share the same variance.
    varofind = rep(varofind_default, k)
  )
  return(invisible(list(hp = dhp, ip = dip)))
}

set_hi_params <- function(k, j, g,
                          cnt, s, cond, ind,
                          scale = 1.96^2) {
  murnm <- c("mu0", "r0")
  mucondnm <- str_glue_vec("mu_cond", 2)
  muindnm <- str_glue_vec("mu_ind", k)

  dhip <- get_default_hi_params(k, j, g)
  mat <- fit_mg_snb(
    cnt = cnt, s = s, cond = cond, ind = ind,
    snbm = snbm_for_mur,
    snbm_for_mucond = snbm_for_mucond,
    murnm = murnm, mucondnm = mucondnm,
    muindnm = muindnm
  )
  r <- init_hbnb_params(mat,
    murnm = murnm, mucondnm = mucondnm,
    muindnm = muindnm, scale = scale
  )

  ## update the gene scale log mean expression estimation
  dhip$hp$mu0 <- r$hp$mu0

  for (n in nm_params) {
    dhip$ip[[n]] <- r$init[[n]]
  }
  return(invisible(list(hp = dhip$hp, ip = dhip$ip)))
}

to_hbnb_data <- function(cnt, ind, cond, s, hp) {
  ## given the basic data, we translate it into
  ## what hbnb needs.

  return(invisible(c(list(
    n = ncol(cnt),
    k = max(ind),
    j = max(cond),
    g = nrow(cnt),
    s = s,
    cond = cond,
    ind = ind,
    y = t(cnt)
  ), hp)))
}

mu_transform_from_raw <- function(raw_mu, mu0, varofmu) {
  return(invisible(raw_mu * sqrt(varofmu) + mu0))
}

mucond_transfrom_from_raw <- function(raw_mu_cond, varofcond) {
  return(invisible(raw_mu_cond * sqrt(varofcond)))
}

muind_transform_from_raw <- function(raw_mu_ind, varofind) {
  ## raw_mu_ind: g by k
  ## varofind: k by 1
  return(invisible(raw_mu_ind %*% diag(sqrt(varofind))))
}

vi_mu_transform_from_raw <- function(vi_raw_mu, mu0,
                                     vi_varofmu, genenms = NULL) {
  ## get mu: n by g

  ## vi_raw_* are got by vi_fit$draw()

  ## element-wise multiply, when second is column
  ## it will column by column.
  res <- vi_raw_mu * sqrt(as.numeric(vi_varofmu))
  if (!is.null(genenms)) {
    colnames(res) <- genenms
  }
  return(invisible(res))
}

split_matrix_col <- function(mat, second_dim, second_dim_nms = NULL) {
  ## split matrix col into a matrix row-major order

  t <- vapply(seq_len(nrow(mat)), function(i) {
    ## byrow = TRUE, means we order the elements row by row.
    ## so each time we use ncol element to fill a row.
    matrix(mat[i, ], nrow = second_dim, byrow = TRUE)
  }, FUN.VALUE = matrix(mat[1, ], nrow = second_dim, byrow = TRUE))
  r <- aperm(t, c(3, 1, 2))
  if (!is.null(second_dim_nms)) {
    dimnames(r)[[2]] <- second_dim_nms
  }
  return(invisible(r))
}

vi_mucond_transform_from_raw <- function(vi_raw_mu_cond, vi_varofcond,
                                         g, genenms = NULL) {
  ## get mucond: n by g by j (num_of_cond, default is 2)

  ## vi_raw_mu_cond shape: n by g*j
  ## order of the names of vi_raw_mu_cond: by row

  ## vi_varofcond shape: n by 1 if we use cmdstanr draw matrix
  ## then let it be a vector length of n
  t <- vi_raw_mu_cond * sqrt(as.numeric(vi_varofcond))
  r <- split_matrix_col(t, second_dim = g, second_dim_nms = genenms)
  return(invisible(r))
}

vi_muind_transform_from_raw <- function(vi_raw_mu_ind, vi_varofind, g, k,
                                        genenms = NULL) {
  ## vi_raw_mu_ind: n by g * k
  ## each row: [k1, k2, ..., kk] of g1, then of g2, ...

  ## vi_varofind: n by k
  ## repeat_var_per_gene: n by g * k
  ## each row in order:
  ## [k1, k2, ..., kk], this order repeat g times.
  ## [MAJOR] use "times" instead of "each" in rep function.
  repeat_var_per_gene <- t(apply(vi_varofind, 1, rep, times = g))
  ## element-wise multiplification
  t <- vi_raw_mu_ind * sqrt(repeat_var_per_gene)

  ## order the col in row-order
  ## then each row: [k1, k2, ..., kk]
  r <- split_matrix_col(t, second_dim = g, second_dim_nms = genenms)
  return(r)
}

run_hbnb_vi <- function(data, ip,
                        adapt_engaged = TRUE,
                        adapt_iter = adapt_iter,
                        eta = eta,
                        algorithm = vi_algorithm,
                        seed = 1L) {
  ## adapt_iter: 5 (default in cmdstan) * adapt_iter we set
  invisible(
    hbnbm$variational(
      data = data,
      init = list(ip),
      seed = seed,
      refresh = vi_refresh,
      iter = num_iter,
      eval_elbo = eval_elbo,
      adapt_engaged = adapt_engaged,
      adapt_iter = adapt_iter,
      algorithm = algorithm,
      output_samples = output_samples,
      tol_rel_obj = tol_rel_obj,
      eta = eta
    )
  )
}

extract_vifit <- function(vifit,
                          data,
                          param) {
  ## get the draw matrix for the param
  ## if not exist, return NaN

  if (!(param %in% nm_params)) {
    message(stringr::str_glue("{param} is not recognized."))
    return(NaN)
  }
  if (param %in% c("varofmu", "varofcond")) {
    return(invisible(vifit$draws(param)))
  }
  if (param %in% c("hp_r", "hp_varofind")) {
    return(invisible(vifit$draws(str_glue_vec(param, 2))))
  }

  gbc <- t(data$y)
  if (!is.null(rownames(gbc))) {
    genenms <- rownames(gbc)
  } else {
    genenms <- seq_len(data$g)
  }

  if (param == "nb_r") {
    r <- vifit$draws(str_glue_vec(param, data$g))
    colnames(r) <- genenms
    return(invisible(r))
  }

  if (param == "varofind") {
    r <- vifit$draws(str_glue_vec(param, data$k))
    return(invisible(r))
  }


  if (param %in% c("raw_mu", "mu")) {
    r <- vifit$draws(str_glue_vec("raw_mu", data$g))
    if (param == "mu") {
      varofmu <- vifit$draws("varofmu")
      mu0 <- data$mu0
      return(invisible(vi_mu_transform_from_raw(r, mu0, varofmu, genenms)))
    }
    colnames(r) <- genenms
    return(invisible(r))
  }

  if (param %in% c("raw_mu_cond", "mu_cond")) {
    t1 <- vifit$draws(str_glue_mat("raw_mu_cond", nr = data$g, nc = data$j))
    if (param == "mu_cond") {
      t2 <- vifit$draws("varofcond")
      r <- vi_mucond_transform_from_raw(t1, t2, data$g, genenms)
      return(invisible(r))
    }
    return(invisible(split_matrix_col(t1, data$g, genenms)))
  }

  if (param %in% c("raw_mu_ind", "mu_ind")) {
    t1 <- vifit$draws(str_glue_mat("raw_mu_ind", nr = data$g, nc = data$k))
    if (param == "mu_ind") {
      t2 <- vifit$draws(str_glue_vec("varofind", data$k))
      r <- vi_muind_transform_from_raw(t1, t2,
        g = data$g, k = data$k,
        genenms
      )
      return(invisible(r))
    }
    return(invisible(split_matrix_col(t1, data$g, genenms)))
  }
  warning(paste(param, " is missed."))
  return(NaN)
}

get_rank_statistics <- function(mu_cond, c1 = 1, c2 = 2,
                                std_cond = 1) {
  ## mu_cond: n by ngene by ncond
  ## c1, c2 correspond to different conditions

  delta <- as.matrix(mu_cond[, , c1] - mu_cond[, , c2])
  n <- nrow(delta)
  if (!is.null(dimnames(mu_cond)[[2]])) {
    colnames(delta) <- dimnames(mu_cond)[[2]]
  }
  sd_delta <- matrixStats::colSds(delta)
  abs_mean_delta <- abs(colMeans(delta))
  abs_delta <- abs(delta)

  ## z_score (actually t_score)
  z <- abs_mean_delta * sqrt(n) / sd_delta
  m <- abs_mean_delta

  ## probability larger than a given epision
  ## sqrt of var_of_cond
  p <- colSums(abs_delta >= std_cond) / n

  p10 <- colSums(abs_delta >= (std_cond * 1.28)) / n
  p05 <- colSums(abs_delta >= (std_cond * 1.645)) / n
  p025 <- colSums(abs_delta >= (std_cond * 1.96)) / n

  p0 <- colSums(abs_delta > 0.0) / n
  bf <- abs(log(p0 + 1e-06) - log(1 - p0 + 1e-06))

  return(invisible(list(
    z = z, p = p, delta = delta,
    p10 = p10,
    p05 = p05,
    p025 = p025,
    bf = bf,
    m = m
  )))
}

extract_all_params_from_fit <- function(vifit, data) {
  est_params <- lapply(nm_params, function(nm) {
    extract_vifit(vifit, data, nm)
  })
  names(est_params) <- nm_params
  return(invisible(est_params))
}

get_auc <- function(rank_stats, diffg, ndiffg) {
  true_class <- c(rep(TRUE, length(diffg)), rep(FALSE, length(ndiffg)))
  z <- rank_stats$z[c(diffg, ndiffg)]
  m <- rank_stats$m[c(diffg, ndiffg)]
  p <- rank_stats$p[c(diffg, ndiffg)]
  p10 <- rank_stats$p10[c(diffg, ndiffg)]
  p05 <- rank_stats$p05[c(diffg, ndiffg)]
  p025 <- rank_stats$p025[c(diffg, ndiffg)]
  bf <- rank_stats$p[c(diffg, ndiffg)]


  auc_z <- caTools::colAUC(z, true_class)
  auc_m <- caTools::colAUC(m, true_class)
  auc_p <- caTools::colAUC(p, true_class)
  auc_p10 <- caTools::colAUC(p10, true_class)
  auc_p05 <- caTools::colAUC(p05, true_class)
  auc_p025 <- caTools::colAUC(p025, true_class)

  auc_bf <- caTools::colAUC(bf, true_class)
  return(invisible(list(
    auc_z = auc_z,
    auc_m = auc_m,
    auc_p = auc_p,
    auc_p10 = auc_p10,
    auc_p05 = auc_p05,
    auc_p025 = auc_p025,
    auc_bf = auc_bf
  )))
}

get_auc_v2 <- function(rank_stats, diffg, ndiffg) {
  ## simplify the auc calculations
  true_class <- c(rep(TRUE, length(diffg)), rep(FALSE, length(ndiffg)))
  z <- rank_stats$z[c(diffg, ndiffg)]
  m <- rank_stats$m[c(diffg, ndiffg)]
  p <- rank_stats$p[c(diffg, ndiffg)]
  bf <- rank_stats$p[c(diffg, ndiffg)]


  auc_z <- caTools::colAUC(z, true_class)
  auc_m <- caTools::colAUC(m, true_class)
  auc_p <- caTools::colAUC(p, true_class)

  auc_bf <- caTools::colAUC(bf, true_class)
  return(invisible(list(
    auc_z = auc_z,
    auc_m = auc_m,
    auc_p = auc_p,
    auc_bf = auc_bf
  )))
}

## * test
test <- function() {
  pbmc <- readRDS(here::here(
    "src", "modelcheck",
    "snb_pool_ref_pbmc.rds"
  ))
  k <- max(pbmc$ind)
  j <- 2
  g <- nrow(pbmc$y2c)

  hi_params <- set_hi_params(
    k = k, j = j, g = g,
    cnt = pbmc$y2c, s = pbmc$s, cond = pbmc$cond,
    ind = pbmc$ind, scale = 1.96^2
  )

  data <- to_hbnb_data(pbmc$y2c,
    ind = pbmc$ind, cond = pbmc$cond,
    s = pbmc$s, hp = hi_params$hp
  )

  vi_sampler <- run_hbnb_vi(data = data, ip = hi_params$ip)
  est_params <- lapply(nm_params, function(nm) {
    extract_vifit(vi_sampler, data, nm)
  })
  names(est_params) <- nm_params
}
