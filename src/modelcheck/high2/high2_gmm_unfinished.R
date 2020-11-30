## high2 model

## param names in model
## 1) names defined in  stan script (by us)
## 2) names returned by stan sampler/opt
##    (stan adds [i,j] or [i] for matrix/vec)

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
library(R6)

## warnings/errors traceback settings
options(error = traceback)
options(warn = 1)
options(mc.cores = 3)

High2 <- R6Class(
  classname = "High2",
  public = list(
    ## stan script
    stan_snb_path = NULL,
    stan_snb_cond_path = NULL,
    stan_high2_path = NULL,
    ## num of inds
    nind = NULL,
    ## stan model
    snb = NULL,
    snb_cond = NULL,
    high2 = NULL,
    ## vi training parameters
    num_iter = NULL,
    opt_refresh = NULL,
    vi_refresh = NULL,
    algorithm = NULL,
    eval_elbo = NULL,
    output_samples = NULL,
    tol_rel_obj = NULL,
    adapt_iter = NULL,
    ## default hyper parameters
    r = NULL,
    mu = NULL,
    varofmu = NULL,
    min_varofmu = NULL,
    varofcond = NULL,
    min_varofcond = NULL,
    varofind = NULL,
    sd_init_muind = NULL,
    sd_init_mucond = NULL,
    ## other params with default values
    scale_sd_mu = 1.96,
    scale_sd_mucond = 1.96,
    opt_iter = 5000,
    seed = 1L,
    ## high2 parameter names
    murnm = c("mu", "r"),
    mucondnm = "mu_cond",
    ## need initialize
    muindnm = NULL,
    all_params_nms = c(
      "nb_r", "hp_r", "varofmu", "centerofmu", "mu",
      "varofcond", "mu_cond", "weightofcond",
      "hp_varofind", "varofind", "mu_ind",
      "raw_mu", "raw_mu_ind"
    ),
    initialize = function(stan_snb_path,
                          stan_snb_cond_path,
                          stan_high2_path,
                          nind,
                          num_iter = 20000,
                          opt_refresh = 0,
                          vi_refresh = 2000,
                          algorithm = "meanfield",
                          eval_elbo = 100,
                          output_samples = 2000,
                          tol_rel_obj = 0.0001,
                          adapt_iter = 1000,
                          r = 50.0,
                          mu = 0.0,
                          varofmu = 16.0,
                          min_varofmu = 2.0,
                          varofcond = 4.0,
                          min_varofcond = 1.0,
                          varofind = 1.0,
                          sd_init_muind = 0.01,
                          sd_init_mucond = 0.01) {
      self$stan_snb_path <- stan_snb_path
      self$stan_snb_cond_path <- stan_snb_cond_path
      self$stan_high2_path <- stan_high2_path
      self$num_iter <- num_iter
      self$opt_refresh <- opt_refresh
      self$vi_refresh <- vi_refresh
      self$algorithm <- algorithm
      self$eval_elbo <- eval_elbo
      self$output_samples <- output_samples
      self$tol_rel_obj <- tol_rel_obj
      self$adapt_iter <- adapt_iter
      self$r <- r
      self$mu <- mu
      self$varofmu <- varofmu
      self$min_varofmu <- min_varofmu
      self$varofcond <- varofcond
      self$min_varofcond <- min_varofcond
      self$varofind <- varofind
      self$sd_init_muind <- sd_init_muind
      self$sd_init_mucond <- sd_init_mucond
      ## name
      self$nind <- nind
      self$muindnm <- self$str_glue_vec("mu_ind", self$nind)
      ## init stan model
      self$snb <- self$init_stan_model(self$stan_snb_path)
      self$snb_cond <- self$init_stan_model(self$stan_snb_cond_path)
      self$high2 <- self$init_stan_model(self$stan_high2_path)
    },
    init_stan_model = function(model_path) {
      if (file.exists(model_path)) {
        return(invisible(cmdstanr::cmdstan_model(stan_file = model_path,
                                compile = T, quiet = T)))
      } else {
        stop(paste(model_path, "not exist.", sep = " "))
      }
    },
    does_fit_well = function(stanfit) {
      if (stanfit$return_codes() == 0) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    },
    str_glue_vec = function(nm = "MuInd", n = 10,
                            ls = "[", rs = "]") {
      invisible(vapply(seq_len(n), function(i) {
        paste0(nm, ls, i, rs)
      }, FUN.VALUE = "MuInd[1]"))
    },
    str_glue_mat_rowise = function(nm = "MuInd", nr, nc,
                                   ls = "[", rs = "]") {
      r <- rep("", nr * nc)
      for (i in 1:nr) {
        for (j in 1:nc) {
          r[nc * (i-1) + j] <- paste0(nm, ls, i, ",",j, rs)
        }
      }
      return(invisible(r))
    },
    check_y_are_all_zeros = function(y) {
      if (sum(y) < 1) {
        warning("All the y are zeros")
        return(TRUE)
      } else {
        return(FALSE)
      }
    },
    check_s = function(s) {
      if (any(s == 0)) {
        stop("s has zeros. Check the data!")
      }
    },
    ## We might need a matrix version to init all at once
    fallback_init_snb = function(y, s) {
      if (self$check_y_are_all_zeros(y)) {
        warning("Use the dfeault logmu: ", self$mu,
                " and r: ", self$r)
        return(invisibile(list(mu = self$mu, r = self$r)))
      } else {
        m <- mean(y)
        mu <- log(m) - log(median(s))
        v <- var(y)
        r <- ifelse(v>m, m^2/(v-m), self$r)
        return(invisibile(list(mu = mu, r = r)))
      }
    },

    fit_snb = function(y, s) {
      ## mu in scaled log level and furthermore minus log(s)
      self$check_s(s)
      result <- list(mu = self$mu, r = self$r, success = FALSE)
      if (self$check_y_are_all_zeros(y)) {
        return(invisibile(result))
      }
      init_mur <- self$fallback_init_snb(y, s)
      ## get the mle fit
      capture.output(opt <- self$snb$optmize(
        data = list(n = length(y), s = s, y = y),
        seed = self$seed,
        refresh = self$opt_refresh,
        iter = self$opt_iter,
        init = list(init_mur),
        algorithm = "lbfgs"
      ))
      if (!self$dose_fit_well(opt)) {
        warning("Fitting snb failed. Using init mu: ",
          init_mur$mu, " init r: ", init_mur$r)
        result$mu <- init_mur$mu
        result$r <- init_mur$r
        return(invisible(result))
      }
      est_mur <- opt$mle()
      if (est_mur$r > self$big_r) {
        warning("r > ", self$big_r,
                "Will use init mu: ",
                init_mur$mu, " init r: ",
                init_mur$r)
        result$mu <- init_mur$mu
        result$r <- init_mur$r
        return(invisible(result))
      }
      result$mu <- est_mur$mu
      result$r <- est_mur$r
      result$sucess <- TRUE
      return(invisible(result))
    },

    fit_snb_cond = function(y, s, cond, r, mu) {
      self$check_s(s)
      result <- list(mucond = 0.0, success = FALSE)
      if (self$check_y_are_all_zeros(y)) {
        return(invisible(result))
      }
      if (self$check_y_are_all_zeros(y[cond==1]) &
            self$check_y_are_all_zeros(y[cond==2])) {
        warning("In one condition, y are all zeros.")
        return(invisible(result))
      }
      lmu <- log(mean(y[cond==1])) - log(median(s[cond==1]))
      rmu <- log(mean(y[cond==2])) - log(median(s[cond==2]))
      init_mucond <- (lmu - rmu) / 2
      capture.output(opt <- self$snb_cond$optimize(
        data = list(n = length(y), s = s, y = y,
                    cond = cond, r = r, mu = mu),
        seed = self$seed,
        refresh = self$opt_refresh,
        iter = self$opt_iter,
        init = list(list(mucond = init_mucond)),
        algorithm = "lbfgs"
      ))
      if (!self$does_fit_well(opt)) {
        warning("Mucond Fitting failed. Use init: ",
                init_mucond)
        result$mucond <- init_mucond
        return(invisible(result))
      }
      t <- opt$mle()
      result$mucond <- t["mucond"]
      result$success <- TRUE
      return(invisible(result))
    },

    fit_gwsnb_to_cond_level = function(y, s, cond, ind) {
      ## self$check_s(s)
      nind <- max(ind)
      result <- list(
        mu <- self$mu, r = self$r,
        mu_cond <- rnorm(n=1, mean = 0, sd = self$sd_init_mucond),
        mu_ind <- rnorm(n = nind, mean = 0, sd = self$sd_init_muind),
        s1 <- FALSE, s2 <- FALSE
      )
      if (check_y_are_all_zeros(y)) {
        return(invisible(result))
      }
      mur <- self$fit_snb(y=y, s=s)
      result$mu <- mur$mu
      result$r <- mur$r
      result$s1 <- mur$success
      if (!mur$success) {
        return(invisible(result))
      }
      est_mucond <- self$fit_snb_cond(y = y, s = s,
                                  cond = cond, r = reesult$r,
                                  mu = result$mu)
      result$mu_cond <- est_mucond$mucond
      result$s2 <- est_mucond$success
      return(invisible(result))
    },

    est_varofmu = function(mu) {
      v <- var(mu) * (self$scale_sd_mu)^2
      if (is.nan(v)) {
        warning("varofmu is NaN. Set it as default.")
        v <- self$varofmu
      }
      if (v < self$min_varofmu) {
        warning("varofmu: ", v,
                " < ", self$min_varofmu, " and use default.")
        v <- self$min_varofmu
      }
      return(v)
    },

    est_varofcond = function(mucond) {
      ## mucond: g by 1 vector
      v <- max(abs(mucond)) * self$scale_sd_mucond
      v <- v^2
      if (v < self$min_varofcond) {
        warning("varofcond: ", v, " < ", self$min_varofcond,
                " and use the default.")
        v <- self$min_varofcond
      }
      return(invisible(v))
    },

    fit_mg_snb = function(cnt, s, cond, ind) {
      ## return ngene by (2+1+k) estimation
      t_res <- vapply(
        seq_len(nrow(cnt)),
        function(i) {
          r <- self$fit_gwsnb_to_cond_level(
            y = cnt[i,], s = s, cond = cond,
            ind = ind
          )
          return(invisible(c(r$mu, r$r, r$mu_cond, r$mu_ind)))
        }, FUN.VALUE = rep(0.0, 2 + 1 + self$nind))
      ## gene name
      colnames(t_res) <- rownames(cnt)
      ## value name
      rownames(t_res) <- c(self$murnm, self$mucondnm, self$muindnm)
      return(invisible(t(t_res)))
    },

    init_params <- function(est_mg_mat) {
      mu <- est_mg_mat[, self$murnm[1]]
      varofmu <- self$est_varofmu(mu)
      centerofmu <- median(mu)
      raw_mu <- (mu - centerofmu) / sqrt(varofmu)

      r <- est_mg_mat[, self$murnm[2]]

      mu_cond <- est_mg_mat[, self$mucondnm]
      varofcond <- self$est_varofcond(mu_cond)
      raw_mu_cond <- mu_cond / sqrt(varofcond)

      mu_ind <- est_mg_mat[, self$muindnm]
      ## use varofcond to estimate varofind
      varofind <- rep(varofcond, self$nind)
      raw_mu_ind <- mu_ind / sqrt(varofind)

    }
  )
)


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

