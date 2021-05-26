## mssc v2-0
## or named high2 (hierarchical Bayesian based gene-wise model)
## or named hbnb (hiearchical Bayesian based negative binomial model)

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

## * load dependeces
suppressWarnings(suppressMessages({
  ## library(tidyverse)
  library(cmdstanr)
  library(loo)
  library(posterior)
  library(R6)
}))

## warnings/errors traceback settings
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

init_stan_model <- function(model_path,
                            use_thread = NULL,
                            use_mpi = NULL,
                            use_opencl = NULL) {
  ## load stan model
  ## argument
  ## - use_thread: when stan script uses reduce_sum or map_rect
  ##   set this to be TRUE
  ##   - when running, remember to export STAN_NUM_THREADS=4
  ## - use_mpi: used when stan script uses map_rect 
  ## - use_opencl: enable the OpenCL backend by setting
  ##   STAN_OPENCL=true when compiling
  ##   https://mc-stan.org/docs/2_26/cmdstan-guide/parallelization.html
  ##   The Stan model compiled with STAN_OPENCL can also be supplied the
  ##   OpenCL platform and device IDs of the target device. These IDs
  ##   determine the device on which to run the OpenCL-supported functions
  ##   on. You can list the devices on your system using the clinfo
  ##   program. If the system has one GPU and no OpenCL CPU runtime, the
  ##   platform and device IDs of the GPU are typically 0. In that case you
  ##   can also omit the OpenCL IDs as the default 0 IDs are used in that
  ##   case.
  ##   We supply these IDs when starting the executable as shown below:
  ##   path/to/model sample data file=data.json opencl platform=0 device=1

  ## BUG: set use_xxx = NULL when we don't want to use them.
  ##      setting them FALSE will be equal to setting them TRUE
  
  if (file.exists(model_path)) {
    ## TODO: cmdstanr have output message like
    ## "Model executable is up to date"
    return(invisible(cmdstanr::cmdstan_model(
      stan_file = model_path,
      compile = TRUE,
      quiet = TRUE,
      pedantic = TRUE,
      ## FIXME: when use_thread = FALSE, compling procedures still uses thread
      cpp_options = list(stan_threads = use_thread,
                         stan_map = use_mpi,
                         stan_opencl = use_opencl),
      stanc_options = list(),
      force_recompile = FALSE
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
        ## if y are all zeros, then we use the default ones.
        ## TODO: maybe we can set muind as the mean of muind from our model,
        ## which can be directly estimated by a pure non-mucond effect model.
        ## for example, we can use random effect model?
        return(invisible(result))
      }
      ## ** init and opt mu and r
      init_mur <- self$init_snb(y, s)

      ## NOTE: comment the initialization based on
      ## optimization of scaled negative binomial model,
      ## which usually fails based on the simulation,
      ## and the init_snb function makes sense already.
      
      ## capture.output(opt <- self$snb$optimize(
      ##   data = list(
      ##     n = length(y), s = s, y = y,
      ##     hpg = c(self$gamma_alpha, self$gamma_beta)
      ##   ),
      ##   seed = self$seed,
      ##   refresh = self$opt_refresh,
      ##   iter = self$opt_iter,
      ##   init = list(init_mur),
      ##   algorithm = "lbfgs"
      ## ))
      
      ## ** update mu and r
      ## if (does_fit_well(opt)) {
      ##   est_mur <- opt$mle()
      ##   result$mu <- est_mur["mu"]
      ##   ## r might be big due to fitting issue.
      ##   est_r <- est_mur["r"]
      ##   result$r <- ifelse(est_r < self$big_r, est_r, self$r)
      ##   s1 <- ifelse(est_r < self$big_r, TRUE, FALSE)
      ## } else {
      ##   result$mu <- init_mur$mu
      ##   result$r <- init_mur$r
      ## }
      result$mu <- init_mur$mu
      result$r <- init_mur$r
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
    est_varofcond = function(mucond) {
      ## mucond: ngene by ncond
      ncond <- ncol(mucond)
      ngene <- nrow(mucond)
      ## each condiiton has its own variance.
      ## which follows a inv-gamma prior
      t_d <- vapply(
        1:ncond, function(i) {
          t <- max(abs(mucond[, i]))
          ## set a variance not that small
          v <- max(t^2, self$min_varofcond)
          ## assume the mean of mucond is around 0.0
          ## then use posterior of inv-gamma to set the hyper priors
          alpha <- 1.0 + ngene / 2
          beta <- 1.0 + sum(mucond[, i]^2) / 2
          invisible(c(v, alpha, beta))
        },
        FUN.VALUE = rep(1.0, 3)
      )
      return(invisible(t(t_d)))
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
    ### store the fitting result of vi
    high2fit = NULL,
    ### store the fitting result of opt
    high2optfit = NULL,
    ## vi training parameters
    num_iter = NULL,
    vi_refresh = NULL,
    algorithm = NULL,
    eval_elbo = NULL,
    output_samples = NULL,
    tol_rel_obj = NULL,
    adapt_iter = NULL,
    adapt_engaged = NULL,
    eta = NULL,
    ## high2 parameter names
    all_params_nms = c(
      "centerofmu", "varofmu", "mu", "centerofr",
      "varofr", "r", "varofcond", "mucond",
      "tau2", "centerofind", "varofind", "muind"
    ),
    initialize = function( ## gwsnb parameters
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
                          adapt_engaged = FALSE,
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
      self$tol_rel_obj <- tol_rel_obj
      self$adapt_iter <- adapt_iter
      self$adapt_engaged <- adapt_engaged
      self$eta <- eta
      self$nind <- nind
      self$ncond <- ncond
    },

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
        hp_varofcond = init_mgsnb$cond[, 2:3],
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
      ### NOTE: log of r level
      #### Should change the name of r to logr
      centerofr <- init_mgsnb$logr[1]
      varofr <- init_mgsnb$logr[2]
      raw_r <- (log(r) - centerofr) / sqrt(varofr)

      ### ngene by ncond
      mucond <- init_mgsnb$mgsnb[, 3:(2 + self$ncond)]
      ### ncond by 1
      varofcond <- init_mgsnb$cond[, 1]
      raw_mucond <- mucond %*% diag(1 / sqrt(varofcond))

      ### scalar
      tau2 <- init_mgsnb$ind$est_tau2[1]
      ### nind by 1
      centerofind <- init_mgsnb$ind$est_varofind[, 1]
      ### nind by 1
      raw_centerofind <- centerofind / sqrt(tau2)
      ### nind by 1
      varofind <- init_mgsnb$ind$est_varofind[, 2]
      ### ngene by nind
      ## TODO: muind estimated by the genes from DESEQ2.
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
        eta = self$eta)
    }, ## end of run method

    run_opt = function(data,
                       list_wrap_ip = NULL,
                       threads = NULL,
                       refresh = 10,
                       max_iter = 5000,
                       opt_method = "lbfgs",
                       init_alpha = 0.001,
                       tol_obj = NULL,
                       tol_rel_obj = NULL,
                       tol_grad = NULL,
                       tol_rel_grad = NULL,
                       tol_param = NULL,
                       history_size = 10) {
      ## optimization
      ## - Jacobian adjustment is not an issue: Stan turns off the built-in Jacobian
      ##   adjustments for optimization.
      ## - (L-)BFGS is controlled by a number of tolerance values
      ##   - any one of which being satisfied caused the algorithm to terminate.
      ##   - any of the convergence tests can be disabled by setting its
      ##     corresponding tolerance parameter to zero.
      ## - relative convergence tolerance will be multiplied by an epsilon
      ##   - epsilon around 2 * 10^(-16), is machine precision
      ## - for any given optimization probm it's worthwhile trying the program
      ##   with and without constraint to see which oen is more efficient
      ## arguments
      ## - refresh
      ##   - 0: only error messages will be printed
      ##   - n(>0): number of iterations between printed screen updates
      ## - init_alpha: initial step size, default is 0.001
      ## - tol_obj: convergence tolerance on relative changes in object function
      ##   - NULL: default 10^(-12) used in cmdstan
      ## - tol_rel_obj: convergence tolerance on relative changes in object function
      ##   - NULL: default 10^(4) used in cmdstan
      ## - tol_grad: convergence tolerance on the norm of gradient
      ##   - NULL: default 10^(-8) used in cmdstan
      ## - tol_rel_grad: convergence tolerance on relative norm of the gradient
      ##   - NULL: default 10^(7) used in cmdstan
      ## - tol_param: convergence tolerance on changes in parameter value
      ##   - NULL: default 10^(-8) used in cmdstan
      ## - history_size: used for lbfgs to estimate the pseudo Hessian
      ##   - NULL: default 5 used in cmdstan, should be sufficient around5-10
      ##     must smaller than the dim of parameters
      ##     When lbfgs performs poorly but bfgs performs well, consider to increase
      ##     history_size.
      
      self$high2optfit <- self$high2$optimize(
        data = data,
        init = list_wrap_ip,
        seed = self$seed,
        refresh = refresh,
        save_latent_dynamics = FALSE,
        output_dir = NULL,
        sig_figs = NULL,
        threads = threads,
        algorithm = opt_method,
        init_alpha = init_alpha,
        iter = max_iter,
        tol_obj = tol_obj,
        tol_rel_obj = tol_rel_obj,
        tol_grad = tol_grad,
        tol_rel_grad = tol_rel_grad,
        tol_param = tol_param,
        history_size = history_size)
    }, ## end of run_opt method

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
    }, ## end of psis method

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
        if (param %in% c("varofcond")) {
          ## when param is vector of len of ncond
          return(invisible(
            self$high2fit$draws(str_glue_vec(param, self$ncond))
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
      }, ## end of try block
      error = function(e) {
        warning(e)
        return(invisible(NA))
      }) ## end of try-catch
    }, ## end of extract draws
    
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
    }, ## end of extract_draws_all
    
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

    get_auc = function(ranking_statistic, c1, c2) {
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
    } ## end of get_auc method
  ) ## end of public field
) ## end of class high2


## * test
test <- function() {
  pbmc <- readRDS(here::here(
    "src", "archive", "modelcheck",
    "snb_pool_ref_pbmc.rds"
  ))
  nind <- max(pbmc$ind)
  model <- High2$new(
    stan_snb_path = here::here("src", "mssc", "stan", "snb.stan"),
    stan_high2_path = here::here("src", "mssc", "stan", "mssc_2-0.stan"),
    nind = nind,
    tol_rel_obj = 0.0001,
    adapt_engaged = FALSE
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
  
  ## variational inference
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
  rankings <- model$get_ranking_statistics(
    mucond = mucond,
    two_hot_vec = c(1, -1)
  )
  str(rankings)
  
  ## optimization
  model$run_opt(data = data, list_wrap_ip = list(init_params$ip))
}

test()
