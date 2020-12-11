## mssc 2-1 model

## This inherits from high2
## But model the gene-specific variances of conditions,
## which is actually used in the first-version of mssc.


## * load depenences
library(posterior)
library(cmdstanr)
library(loo)
library(R6)

## * settings
options(error = traceback)
options(warn = 1)
options(mc.cores = 2)

## * common functions
init_snb_log_mean <- function(y, s) {
  ## this is similar with log(mean(y/s)) (or better)
  ## y should be: not all zeros
  ## s should no zeros.
  mu <- log(mean(y)) - log(median(s))
  invisible(mu)
}

init_stan_model <- function(model_path) {
  ## load stan model
  if (file.exists(model_path)) {
    return(invisible(cmdstanr::cmdstan_model(
      stan_file = model_path,
      compile = T, quiet = T
    )))
  } else {
    stop(paste(model_path, "not exist.", sep = " "))
  }
}

does_fit_well <- function(stanfit) {
  ## check if the model fit well
  invisible(ifelse(stanfit$return_codes() == 0, TRUE, FALSE))
}

str_glue_vec <- function(nm = "MuInd", n = 10,
                         ls = "[", rs = "]") {
  ## get names of vector for cmdstanr
  invisible(vapply(seq_len(n), function(i) {
    paste0(nm, ls, i, rs)
  }, FUN.VALUE = "MuInd[1]"))
}

str_glue_mat_rowise <- function(nm = "MuInd", nr, nc,
                                ls = "[", rs = "]") {
  ## get names of matrix for cmdstanr
  ## row-wise major
  r <- rep("", nr * nc)
  for (i in 1:nr) {
    for (j in 1:nc) {
      r[nc * (i - 1) + j] <- paste0(nm, ls, i, ",", j, rs)
    }
  }
  return(invisible(r))
}

check_y_are_all_zeros <- function(y) {
  ## return TRUE if y are all zeros
  if (sum(y) < 1) {
    warning("All the y are zeros")
    return(TRUE)
  } else {
    return(FALSE)
  }
}

check_s <- function(s) {
  ## make sure s has no zeros.
  if (any(s == 0)) {
    stop("s has zeros. Check the data!")
  }
}

rep_row <- function(x, n) {
  matrix(rep(x, each = n), nrow = n, byrow = FALSE)
}

rep_col <- function(x, n) {
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

split_matrix_col <- function(mat, second_dim, second_dim_nms = NULL) {
  ## split matrix col into a matrix row-major order
  ## mat: nsample by (second_dim * third_dim)
  ## return: nsample by second_dim by third_dim (i.e., a three-dim array)

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

get_auc <- function(ranking_statistic, c1, c2) {
  ## ranking_statistic: a vector, ngene by 1
  ## c1: index of gene for condition one
  ## c2: index of gene for condition two
  ## return: a vector of AUC value for different columns
  if (is.matrix(ranking_statistic)) {
    t <- ranking_statistic
  } else {
    t <- as.matrix(ranking_statistic, ncol = 1)
  }
  true_class <- c(rep(TRUE, length(c1)), rep(FALSE, length(c2)))
  return(invisible(caTools::colAUC(
    t[c(c1, c2), ],
    true_class
  )))
}

## * Genewisenbfit: initialize the parameters
Genewisenbfit <- R6::R6Class(
  classname = "Genewisenbfit", public = list(
    ## stan models for fitting
    snb = NULL,
    ## gamma alpha and beta for r as hyper prior in stan snb fit
    gamma_alpha = NULL,
    gamma_beta = NULL,
    ## init parameters for snb
    mu = NULL,
    r = NULL,
    big_r = NULL,
    min_varofmu = NULL,
    min_varofcond = NULL,
    min_varofind = NULL,
    min_tau2 = NULL,
    ## parameters for optimization
    seed = 1L,
    opt_iter = 5000,
    opt_refresh = 0,
    initialize = function(stan_snb_path,
                          gamma_alpha = 0.05,
                          gamma_beta = 0.05,
                          mu = 0.0,
                          r = 20,
                          big_r = 500,
                          min_varofmu = 4.0,
                          min_varofcond = 0.25,
                          min_varofind = 0.25,
                          min_tau2 = 0.25) {
      self$snb <- init_stan_model(stan_snb_path)
      self$gamma_alpha <- gamma_alpha
      self$gamma_beta <- gamma_beta
      self$mu <- mu
      self$r <- r
      self$big_r <- big_r
      self$min_varofmu <- min_varofmu
      self$min_varofcond <- min_varofcond
      self$min_varofind <- min_varofind
    },

    init_snb = function(y, s) {
      ## init params for scaled negative binomial:
      ## mean and dispersion
      ## dispersion is defined as the same as in R or stan.
      ## No check if y are all zeros.
      ## Note:
      ## - From observation, mean can be estimate well,
      ##   but if it's below 1, then r tends to be estimated too small.
      mu <- init_snb_log_mean(y = y, s = s)
      ## v equals to m + m^2/r
      v <- var(y)
      m <- mean(y)
      ## dispersion
      r <- ifelse(v > m, m^2 / (v - m), self$r)
      return(invisible(list(mu = mu, r = r)))
    },

    fit_gwsnb = function(y, s, cond, ind) {
      ## fit mu, mucond, muind under scaled negative binomial dist
      ## - mu: scaled log level and furthermore minus log(s)
      ## - r: dispersion
      nind <- max(ind)
      ncond <- max(cond)
      result <- list(
        mu = self$mu,
        r = self$r,
        mucond = rep(0.0, ncond),
        muind = rep(0.0, nind)
      )
      ## state of optimization
      s1 <- FALSE
      if (check_y_are_all_zeros(y)) {
        return(invisible(result))
      }
      ## ** init and opt mu and r
      init_mur <- self$init_snb(y, s)
      capture.output(opt <- self$snb$optimize(
        data = list(
          n = length(y), s = s, y = y,
          hpg = c(self$gamma_alpha, self$gamma_beta)
        ),
        seed = self$seed,
        refresh = self$opt_refresh,
        iter = self$opt_iter,
        init = list(init_mur),
        algorithm = "lbfgs"
      ))
      ## ** update mu and r
      if (does_fit_well(opt)) {
        est_mur <- opt$mle()
        result$mu <- est_mur["mu"]
        ## r might be big due to fitting issue.
        est_r <- est_mur["r"]
        result$r <- ifelse(est_r < self$big_r, est_r, self$r)
        s1 <- ifelse(est_r < self$big_r, TRUE, FALSE)
      } else {
        result$mu <- init_mur$mu
        result$r <- init_mur$r
      }
      ## ** set mucond
      result$mucond <- vapply(
        1:ncond, function(i) {
          ss <- s[cond == i]
          yy <- y[cond == i]
          if (check_y_are_all_zeros(yy)) {
            return(0.0)
          }
          t <- init_snb_log_mean(y = yy, s = ss)
          init_mucond <- t - result$mu
          ## Since fitting is hard when fix r, and limited data
          ## we directly use the init_mucond
          invisible(init_mucond)
        },
        FUN.VALUE = 0.0
      )
      ## ** set muind
      result$muind <- vapply(
        1:nind, function(i) {
          yy <- y[ind == i]
          if (check_y_are_all_zeros(yy)) {
            return(0.0)
          }
          t <- init_snb_log_mean(y = yy, s = s[ind == i])
          return(invisible(t - result$mu - result$mucond[cond[ind == i][1]]))
        },
        FUN.VALUE = 0.0
      )
      return(invisible(result))
    },
    est_varofmu = function(mu) {
      ## mu: ngene by 1
      ## return: mean, var, gamma_alpha, gamma_beta
      m <- median(mu)
      ngene <- length(mu)
      v <- sum((mu - m)^2) / ngene
      v <- max(v, self$min_varofmu)
      ## varofmu prior follows a inv-gamma dist
      ## est hy based on posterior
      alpha <- 1.0 + ngene / 2
      beta <- 1.0 + sum((mu - m)^2) / 2
      return(invisible(c(m, v, alpha, beta)))
    },
    est_varofr = function(r) {
      ## r: ngene by 1
      ## r prior: log normal
      ## return: log level of mean, var, gamma_alpha, gamma_beta
      logr <- log(r)
      invisible(self$est_varofmu(logr))
    },

    ## DIFFERENT with high2
    est_varofcond = function(mucond) {
      ## mucond: ngene by ncond
      ncond <- ncol(mucond)
      ngene <- nrow(mucond)

      ## each gene has its own variance
      ## for different conditions
      varofcond <- vapply(seq_len(ngene), function(i) {
        ## mucond follows a zero mean
        t <- max(abs(mucond[i, ]))
        v <- max(t^2, self$min_varofcond)
      }, FUN.VALUE = 1.0)
      ## varofcond follows inv-gamma(alpha, beta)
      ## est alpha and beta
      ## instead of estimating them from mean and variance of varofcond
      ## here we use the conjugate posterior settings.
      alpha <- 1.0 + ncond / 2
      beta <- 1.0 + quantile(rowSums((mucond - 0.0)^2) / 2, probs = 0.975)
      return(invisible(list(v = varofcond, hp = c(alpha, beta))))
    },

    est_varofind = function(muind) {
      ## muind: ngene by nind
      ## estimate:
      ## - nind by 4: mean of muind, varofmuind, alpha, beta
      ## - tau2, tau2_alpha, tau2_beta
      nind <- ncol(muind)
      ngene <- nrow(muind)
      t_d <- vapply(
        1:nind, function(i) {
          m <- median(muind[, i])
          v <- max(sum((muind[, i] - m)^2) / ngene, self$min_varofind)
          alpha <- 1.0 + ngene / 2
          beta <- 1.0 + sum((muind[, i] - m)^2) / 2
          invisible(c(m, v, alpha, beta))
        },
        FUN.VALUE = rep(1.0, 4)
      )
      ## shape: nind by 4
      r <- t(t_d)
      ## assume muinds follow a N(0.0, tau) (tau is sd)
      tau2 <- max(max(abs(r[, 1]))^2, self$min_tau2)
      ## assume tau2 has a inv-gamma prior
      ## use posterior to set up the hp.
      tau2_alpha <- 1.0 + nind / 2
      tau2_beta <- 1.0 + sum(muind[, 1]^2) / 2
      return(invisible(
        list(
          est_varofind = r,
          est_tau2 = c(tau2, tau2_alpha, tau2_beta)
        )
      ))
    },
    ## this is the function mssc want to use
    fit_mgsnb = function(cnt, s, cond, ind) {
      ## cnt: ngene by ncell
      ## s: ncell by 1; cond: ncell by 1; ind: ncell by 1
      check_s(s)
      ncond <- max(cond)
      nind <- max(ind)
      ngene <- nrow(cnt)
      t_init_mgsnb <- vapply(
        1:ngene, function(i) {
          r <- self$fit_gwsnb(
            y = cnt[i, ], s = s,
            cond = cond, ind = ind
          )
          invisible(unlist(r))
        },
        FUN.VALUE = rep(0.0, 2 + ncond + nind)
      )
      ## shape: ngene by 2 + ncond + nind
      init_mgsnb <- t(t_init_mgsnb)
      init_varofmu <- self$est_varofmu(init_mgsnb[, 1])
      init_varofr <- self$est_varofr(init_mgsnb[, 2])
      ## a list of two elements
      init_varofcond <- self$est_varofcond(init_mgsnb[, 3:(2 + ncond)])
      init_varofind <- self$est_varofind(
        init_mgsnb[, (2 + ncond + 1):ncol(init_mgsnb)]
      )
      return(invisible(list(
        mgsnb = init_mgsnb,
        mu = init_varofmu,
        logr = init_varofr,
        cond = init_varofcond,
        ind = init_varofind
      )))
    }
  ) ## end of publich field
) ## end of Genewisenbfit define

High2 <- R6::R6Class(
  classname = "High2", public = list(
    ## num of inds
    nind = NULL,
    ncond = NULL,
    ## gwsnb model
    gwsnb = NULL,
    ## stan model
    high2 = NULL,
    ### store the fitting result
    high2fit = NULL,
    ## vi training parameters
    num_iter = NULL,
    vi_refresh = NULL,
    algorithm = NULL,
    eval_elbo = NULL,
    output_samples = NULL,
    ## DIFFERENT with high2
    ## used for replacemence sampling
    ndraws = NULL,
    tol_rel_obj = NULL,
    adapt_iter = NULL,
    adapt_engaged = NULL,
    eta = NULL,
    ## parameter names
    ## DIFFERENT with high2
    ## add another variable hp_varofcond
    all_params_nms = c(
      "centerofmu", "varofmu", "mu", "centerofr",
      "varofr", "r", "varofcond", "mucond",
      "tau2", "centerofind", "varofind", "muind",
      "hp_varofcond"
    ),
    initialize = function(## gwsnb parameters
                          stan_snb_path,
                          gamma_alpha = 0.05,
                          gamma_beta = 0.05,
                          r = 20,
                          mu = 0.0,
                          big_r = 500,
                          min_varofmu = 2.0,
                          min_varofcond = 0.25,
                          min_varofind = 0.25,
                          min_tau2 = 0.25,
                          ## high2 related parameters
                          stan_high2_path,
                          nind,
                          ncond = 2,
                          num_iter = 20000,
                          vi_refresh = 2000,
                          algorithm = "meanfield",
                          eval_elbo = 100,
                          output_samples = 2000,
                          tol_rel_obj = 0.0001,
                          adapt_iter = 200,
                          adapt_engaged = TRUE,
                          eta = 0.1) {
      ## initiolize class members
      self$gwsnb <- Genewisenbfit$new(
        stan_snb_path = stan_snb_path,
        mu = mu,
        r = r,
        big_r = big_r,
        min_varofmu = min_varofmu,
        min_varofcond = min_varofcond,
        min_varofind = min_varofind,
        min_tau2 = min_tau2,
        gamma_alpha = gamma_alpha,
        gamma_beta = gamma_beta
      )
      self$high2 <- init_stan_model(stan_high2_path)
      self$num_iter <- num_iter
      self$vi_refresh <- vi_refresh
      self$algorithm <- algorithm
      self$eval_elbo <- eval_elbo
      self$output_samples <- output_samples
      self$ndraws  <- floor(0.75 * self$output_samples)
      self$tol_rel_obj <- tol_rel_obj
      self$adapt_iter <- adapt_iter
      self$adapt_engaged <- adapt_engaged
      self$eta <- eta
      self$nind <- nind
      self$ncond <- ncond
    },

    ## DIFFERENT with high2
    init_params = function(cnt, s, cond, ind) {
      ## init paramters including hyper params.
      init_mgsnb <- self$gwsnb$fit_mgsnb(
        cnt = cnt, s = s,
        cond = cond, ind = ind
      )
      ## * set hyper params
      hp <- list(
        hp_varofmu = init_mgsnb$mu[3:4],
        hp_varofr = init_mgsnb$logr[3:4],
        ## DIFFERENT with high2
        ## hp_varofcond = init_mgsnb$cond[, 2:3],
        hp_varofind = init_mgsnb$ind$est_varofind[, 3:4],
        hp_tau2 = init_mgsnb$ind$est_tau2[2:3]
      )
      ## * set params initiolization
      centerofmu <- init_mgsnb$mu[1]
      varofmu <- init_mgsnb$mu[2]
      mu <- init_mgsnb$mgsnb[, 1]
      raw_mu <- (mu - centerofmu) / sqrt(varofmu)

      ### r in negative binomial
      r <- init_mgsnb$mgsnb[, 2]
      ### log of r level
      centerofr <- init_mgsnb$logr[1]
      varofr <- init_mgsnb$logr[2]
      raw_r <- (log(r) - centerofr) / sqrt(varofr)

      ### ngene by ncond
      mucond <- init_mgsnb$mgsnb[, 3:(2 + self$ncond)]

      ## DIFFERENT with high2
      hp_varofcond <- init_mgsnb$cond$hp
      varofcond <- init_mgsnb$cond$v
      raw_mucond <- diag(1 / sqrt(varofcond)) %*% mucond

      ### scalar
      tau2 <- init_mgsnb$ind$est_tau2[1]
      ### nind by 1
      centerofind <- init_mgsnb$ind$est_varofind[, 1]
      ### nind by 1
      raw_centerofind <- centerofind / sqrt(tau2)
      ### nind by 1
      varofind <- init_mgsnb$ind$est_varofind[, 2]
      ### ngene by nind
      muind <- init_mgsnb$mgsnb[, (2 + self$ncond + 1):
      (2 + self$ncond + self$nind)]
      ngene <- nrow(cnt)
      raw_muind <- (muind - rep_row(centerofind, n = ngene)) %*%
        diag(1 / sqrt(varofind))
      return(invisible(
        list(
          hp = hp,
          ip = list(
            centerofmu = centerofmu,
            varofmu = varofmu,
            raw_mu = raw_mu,
            centerofr = centerofr,
            varofr = varofr,
            raw_r = raw_r,
            varofcond = varofcond,
            hp_varofcond = hp_varofcond,
            raw_mucond = raw_mucond,
            tau2 = tau2,
            raw_centerofind = raw_centerofind,
            varofind = varofind,
            raw_muind = raw_muind
          )
        )
      ))
    },
    to_model_data = function(cnt, ind, cond, s, hp) {
      ## given the basic data, we translate it into
      ## what hbnb needs.
      invisible(c(list(
        ncell = ncol(cnt),
        nind = max(ind),
        ncond = max(cond),
        ngene = nrow(cnt),
        s = s, cond = cond,
        ind = ind, y = t(cnt)
      ), hp))
    },
    run = function(data, list_wrap_ip = NULL) {
      ## set the result of high2 to high2fit
      ## adapt_iter: 5 (default in cmdstan) * adapt_iter we set
      self$high2fit <- self$high2$variational(
        data = data,
        init = list_wrap_ip,
        seed = self$seed,
        refresh = self$vi_refresh,
        iter = self$num_iter,
        eval_elbo = self$eval_elbo,
        adapt_engaged = self$adapt_engaged,
        adapt_iter = self$adapt_iter,
        algorithm = self$algorithm,
        output_samples = self$output_samples,
        tol_rel_obj = self$tol_rel_obj,
        eta = self$eta
      )
    },

    psis = function(takelog = FALSE, donormalize = TRUE) {
      ## Get the Pareto Smoothed Importance Sampling (PSIS) weights.
      ## Weights are normalized and NOT in scale level by default.
      if (is.null(self$high2fit)) {
        warning("High2 has not been run.")
        return(invisible(NA))
      }
      log_ratios <- self$high2fit$lp() -
        self$high2fit$lp_approx()

      ## suppressWarnings for psis
      capture.output(suppressWarnings(r <- loo::psis(
        log_ratios = log_ratios,
        r_eff = NA
      )))
      normweights <- weights(r, log = takelog, normalize = donormalize)
      invisible(list(
        psis = r,
        normweights = normweights
      ))
    },

    extract_draws = function(param,
                             ngene = NULL,
                             genenms = NULL) {
      ## extract draws from model given the param name
      ## after getting the fit
      ## this depends on cmdstanr method, but mofidy the names we need
      if (is.null(self$high2fit)) {
        warning("High2 has not been run.")
        return(invisible(NA))
      }
      check_ngene <- function() {
        if (is.null(ngene)) {
          warning("num of gene is not know, plz set it.")
          return(invisible(NA))
        }
      }
      tryCatch({
          if (param %in% c(
            "centerofmu", "varofmu",
            "centerofr", "varofr", "tau2"
          )) {
            ## when param is scalar
            return(invisible(self$high2fit$draws(param)))
          }
          if (param %in% c("mu", "r", "nb_r")) {
            ## when param is vector len of ngene
            check_ngene()
            t <- self$high2fit$draws(str_glue_vec(param, ngene))
            if (!is.null(genenms)) {
              names(t) <- genenms
            }
            return(invisible(t))
          }

          ## DIFFERENT with high2
          if (param %in% c("hp_varofcond")) {
            return(invisible(self$high2fit$draws(str_glue_vec(param, 2))))
          }
          if (param %in% c("varofcond")) {
            ## when param is vector of len of ncond
            return(invisible(
              ## DIFFERENT with high2
              self$high2fit$draws(str_glue_vec(param, ngene))
            ))
          }

          if (param %in% c("varofind", "centerofind")) {
            ## when param is vector of len of nind
            return(invisible(
              self$high2fit$draws(str_glue_vec(param, self$nind))
            ))
          }

          if (param %in% c("mucond")) {
            ## when param is a matrix of ngene by ncond
            check_ngene()
            t <- self$high2fit$draws(
              str_glue_mat_rowise(param, ngene, self$ncond)
            )
            return(invisible(split_matrix_col(
              mat = t, second_dim = ngene,
              second_dim_nms = genenms
            )))
          }
          if (param %in% c("muind")) {
            ## when param is a matrix of ngene by nind
            check_ngene()
            t <- self$high2fit$draws(
              str_glue_mat_rowise(param, ngene, self$nind)
            )
            return(invisible(split_matrix_col(
              mat = t, second_dim = ngene,
              second_dim_nms = genenms
            )))
          }
        },
        error = function(e) {
          warning(e)
          return(invisible(NA))
        }
      )
    },
    extract_draws_all = function(ngene = NULL,
                                 genenms = NULL) {
      est_params <- lapply(self$all_params_nms, function(nm) {
        self$extract_draws(nm,
          ngene = ngene,
          genenms = genenms
        )
      })
      names(est_params) <- self$all_params_nms
      invisible(est_params)
    },
    get_ranking_statistics = function(mucond, two_hot_vec) {
      ## mucond: nsample by ngene by ncond
      ## two_hot_vec: like (1, -1) or (0, 0, -1, 0, 1, 0)
      ## - i.e., the two conditions we want compare
      ## - one is positive, and the other one is negative
      ## return:
      ## - a matrix: ngene by num_of_statistics with colnames
      if (sum(two_hot_vec) > 1) {
        stop("set 1 and -1 for two conditions.")
      }
      if (length(two_hot_vec) != self$ncond) {
        stop(
          "length of two_hot_vec ",
          length(two_hot_vec), " is not equals to ncond ",
          self$ncond
        )
      }
      ## 3-d array cannot directly mutiply a vector in matrix-multiply way
      ## so we use ours.
      ## r: nsample by ngene
      n <- dim(mucond)[1]
      r <- t(vapply(1:n, function(i) {
        mucond[i, , ] %*% two_hot_vec
      }, FUN.VALUE = rep(0.0, dim(mucond)[2])))

      ## one measure
      abs_colmean <- abs(colMeans(r))

      ## one measure
      ## p0 is the prob that case 1 is largner than case -1
      p0 <- colSums(r > 0.0) / n
      bf <- abs(log(p0 + 1e-06) - log(1 - p0 + 1e-06))

      ## one measure
      ## t statistics
      group1 <- mucond[, , two_hot_vec == 1]
      group2 <- mucond[, , two_hot_vec == -1]
      ngene <- dim(mucond)[2]
      tstat <- vapply(1:ngene, function(i) {
        tryCatch({
          s <- t.test(
            x = group1[, i], y = group2[, i],
            alternative = "two.sided",
            paired = TRUE,
            var.equal = FALSE
          )
          return(invisible(s$statistic))
        },
        error = function(e) {
          warning(e)
          return(invisible(0.0))
        }
        )
      }, FUN.VALUE = 0.0)

      result <- cbind(
        abs_t = abs(tstat),
        bf = bf,
        abs_m = abs_colmean
      )
      if (!is.null(dimnames(mucond)[[2]])) {
        rownames(result) <- dimnames(mucond)[[2]]
      }
      return(invisible(result))
    },

    get_rsis_ranking_statistics = function(ngene, two_hot_vec,
                                           normweights, genenms = NULL,
                                           method = "simple_no_replace"){
      ## method could be "simple_no_replace" or "stratified"
      ## suggest and default: "simple_no_replace"
      ## only for mucond
      t <- self$high2fit$draws(
        str_glue_mat_rowise("mucond", ngene, self$ncond)
      )
      t_rsis <- posterior::resample_draws(
        x = t, weights = normweights,
        method = method, ndraws = self$ndraws
      )

      mucond_rsis <- split_matrix_col(
        mat = t_rsis, second_dim = ngene,
        second_dim_nms = genenms
      )

      invisible(self$get_ranking_statistics(
        mucond = mucond_rsis, two_hot_vec = two_hot_vec
      ))
    },
    get_psis_ranking_statistics = function(mucond, two_hot_vec,
                                           normweights) {
      ## mucond, two_hot_vec: ref get_ranking_statistics.
      ## normweights: getting from psis function, are normalized weights.
      if (sum(two_hot_vec) > 1) {
        stop("Set 1 and -1 for two conditions.")
      }
      if (length(two_hot_vec) != self$ncond) {
        stop(
          "Length of two_hot_vec ",
          length(two_hot_vec), " is not equals to ncond ",
          self$ncond
        )
      }
      len_of_nw <- length(normweights)
      nsample <- dim(mucond)[1]
      if (len_of_nw != nsample) {
        stop("Length of normweights",
             len_of_nm, " != nsample ",
             nsample)
      }
      ngene <- dim(mucond)[2]
      ## r: nsample by ngene
      r <- t(vapply(1:nsample, function(i) {
        mucond[i, , ] %*% two_hot_vec
      }, FUN.VALUE = rep(0.0, 2)))
      ## weighted the samples
      wr <- diag(normweights) %*% r

      ## one measurement
      ## normalized weights, so we directly use sum
      ## to get the weighted sample mean.
      abs_expected_d <- abs(colSums(wr))

      ## one measurement
      ir <- (r > 0.0)
      p0 <- colSums(diag(normweights) %*% ir)
      bf <- abs(log(p0 + 1e-06) - log(1 - p0 + 1e-06))
      ## summary
      result <- cbind(bf = bf, abs_m = abs_expected_d)
      if (!is.null(dimnames(mucond)[[2]])) {
        rownames(result) <- dimnames(mucond)[[2]]
      }
      return(invisible(result))
    }
  ) ## end of public field
) ## end of class high2


test <- function() {
  pbmc <- readRDS(here::here(
    "src", "modelcheck",
    "snb_pool_ref_pbmc.rds"
  ))
  nind <- max(pbmc$ind)
  model <- High2$new(
    stan_snb_path = here::here(
      "src", "mssc_2-1", "stan", "snb.stan"
    ),
    stan_high2_path = here::here(
      "src", "mssc_2-1", "stan", "mssc.stan"
    ),
    nind = nind,
    tol_rel_obj = 0.0001,
    adapt_iter = 200
  )

  init_params <- model$init_params(
    cnt = pbmc$y2c[1:10, ],
    s = pbmc$s,
    cond = pbmc$cond,
    ind = pbmc$ind
  )
  data <- model$to_model_data(
    cnt = pbmc$y2c[1:10, ],
    s = pbmc$s,
    cond = pbmc$cond,
    ind = pbmc$ind,
    hp = init_params$hp
  )
  model$run(data = data, list_wrap_ip = list(init_params$ip))

  est_params <- model$extract_draws_all(
    ngene = 10,
    genenms = rownames(pbmc$y2c[1:10, ])
  )
  str(est_params)

  mucond <- model$extract_draws(
    param = "mucond", ngene = 10,
    genenms = rownames(pbmc$y2c[1:10, ])
  )
  rankings <- model$get_ranking_statistics(
    mucond = mucond,
    two_hot_vec = c(1, -1)
  )
  str(rankings)
  psis <- model$psis()
  print(psis$psis)

  rsis_rankings <- model$get_rsis_ranking_statistics(
    ngene = 10,
    two_hot_vec = c(1, -1),
    genenms = rownames(pbmc$y2c[1:10, ]),
    normweights = psis$normweights
  )

  str(rsis_rankings)
}
test()
