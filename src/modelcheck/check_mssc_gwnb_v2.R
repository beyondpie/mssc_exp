## Check MSSC variational inference model.

## Use simulation based calibration.
## - simulation from prior.
## - [hyper] parameters are setted based on a real dataset:
##   - PBMC dataset

## in the second version, we use stan to support the
## negative binomial fitting.

## * set R environment
import::from(here, here)
suppressPackageStartupMessages(library(tidyverse))
library(MCMCpack)

## develop version of cmdstanr
## devtools::install_github("stan-dev/cmdstanr")
library(cmdstanr)
## use the code below to reload the library
## detach("package:cmdstanr", unload = TRUE)
## library(cmdstanr)

library(grid)
library(gtable)
library(gridExtra)
library(bayesplot)
## for bayesplot plotting
## color_scheme_set("brewer-Spectral")

library(posterior)
library(bbmle)
library(sads)
## install_github("olafmersmann/truncnorm")
library(truncnorm)

## local modules
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
## use the code below to reload the modules
## modules::reload(myfit)
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## warnings/errors traceback settings

options(error = traceback)
options(warn = 0)
options(mc.cores = 3)

## * configs
## for simulation
ncells <- c(50, 100, 150, 200, 300)
ninds <- c(5, 10, 15, 20)
ncond <- 2
sgn <- "NFKB1"
celltype <- "Naive CD4+ T"
scale_factor <- 1e04

## num_of_cell_per_ind <- args$ncell
## num_of_ind_per_cond <- args$nind
## num_of_ind <- num_of_ind_per_cond * num_of_cond

## * load stan models
## if compiling not work, try rebuild_cmdstan()
snbm <- cmdstan_model(
  here::here("stanutils", "scale_nb.stan"),
  compile = T, quiet = T
)

snb_fixr_m <- cmdstan_model(
  here::here("stanutils", "scale_nb_fixed_r.stan"),
  compile = T, quiet = T
)

gwnb_model <- cmdstan_model(here::here(
  "src", "dirty_stan",
  "gwnb_simu_from_prior.stan"
), compile = T, quiet = T)

## * load the dataset
pbmc_seurat <- mypbmc$load_pbmc_seurat() %>%
  mypbmc$extract_from_seurat(pbmc_seurat = .)

subscdata <- mypbmc$get_celltype_specific_scdata(
  pbmc_seurat$cnt,
  pbmc_seurat$resp,
  pbmc_seurat$inds,
  pbmc_seurat$ct,
  ## limit to the cell type
  celltype
)

cnt <- subscdata$cnt
inds <- subscdata$inds
resp <- subscdata$resp
## we can save sumcnt as a pool to sample.
sumcnt <- colSums(cnt) / scale_factor


## functions
fit_singlegene_nb <- function() {
  ## get mu0, mucond, r as in mssc gwnb model
  ## from the real dataset
  ## use stan fit

  d1g_cnt <- cnt[sgn, ]
  outliers <- myfit$is_outlier(d1g_cnt)
  d1g_cnt <- d1g_cnt[!outliers]
  d1g_inds <- inds[!outliers]
  d1g_resp <- resp[!outliers]
  d1g_sumcnt <- sumcnt[!outliers]

  ## fit NB dist
  fit_mur <- myfit$stanfit_gwsnb_till_cond_level(
    s = d1g_sumcnt,
    y = d1g_cnt,
    s_control = d1g_sumcnt[d1g_resp == 0],
    y_control = d1g_cnt[d1g_resp == 0],
    s_case = d1g_sumcnt[d1g_resp == 1],
    y_case = d1g_cnt[d1g_resp == 1],
    scale_nb_model = snbm,
    scale_nb_fixed_r_model = snb_fixr_m
  )
  invisible(list(
    mur = fit_mur,
    cnt = d1g_cnt,
    inds = d1g_inds,
    resp = d1g_resp,
    sumcnt = d1g_sumcnt
  ))
}

set_gwnb_hyper_params <- function(sg_mur, sigmaG0 = 4.0,
                                  default_control_value = 1.0,
                                  default_case_value = -2.0) {
  ## set hyper params based on the fitted sg_mur
  ## from a real dataset

  ## inv-gamma mode estimate.
  ## instead of using var(sg_mur$mu_cond), which might be
  ## small when both of them are in the same sign:
  ## - we assume the mean is around zero

  ## we use the max absolute value, and cover it by 1.5 fold
  mu_control <- sg_mur$mur$mu_cond[1]
  if (is.nan(sg_mur$mur$mu_cond[1])) {
    mu_control <- default_control_value
  }
  mu_case <- sg_mur$mur$mu_cond[2]
  if (is.nan(sg_mur$mur$mu_cond[2])) {
    mu_case <- default_case_value
  }

  alphaTau2G <- 1.0
  betaTau2G <- (alphaTau2G + 1) * max(abs(c(mu_control, mu_case))) * 1.5


  ## when phi2G ~ inv-gamma (alpha, beta)
  ## alphaPhi2G <- 1.0
  ## use mode
  ## betaPhi2G <- min(5.0, (alphaPhi2G + 1) * sg_mur$r0)

  ## when phi2G ~ gamma (alpha, beta)
  betaPhi2G <- 1.0
  ## use mean
  alphaPhi2G <- betaPhi2G * min(10.0, sg_mur$r0)

  invisible(list(
    alphaKappa2G = 1.0,
    betaKappa2G = 1.0,
    alphaTau2G = alphaTau2G,
    betaTau2G = betaTau2G,
    ## control the mean not too small in simulaiton.
    muG0 = max(-7.0, sg_mur$mu0),
    sigmaG0 = sigmaG0,
    alphaPhi2G = alphaPhi2G,
    betaPhi2G = min(5, betaPhi2G)
  ))
}

get_vec_of_repeat_int <- function(to_int, repeatimes, from_int = 1) {
  invisible(c(vapply(
    seq(from_int, to_int),
    FUN = function(i) {
      rep(i, repeatimes)
    }, FUN.VALUE = rep(0L, repeatimes)
  )))
}

simulate_from_gwnb_prior <- function(hp, nind) {
  ## sampling the parameters from the prior
  ## defined by the hyper parameters.

  kappa2g <- min(1.0, rinvgamma(1,
    shape = hp$alphaKappa2G,
    scale = hp$betaKappa2G
  ))
  muind <- rnorm(nind * 2,
    mean = 0.0,
    sd = sqrt(kappa2g)
  )

  ## limit tau size
  tau2g <- min(10.0, rinvgamma(1,
    shape = hp$alphaTau2G,
    scale = hp$betaTau2G
  ))

  ## in reality, mucond should be in different signs.
  mucond <- c(
    max(
      -2.0,
      truncnorm::rtruncnorm(1,
        a = -Inf, b = 0.0, mean = 0.0,
        sd = sqrt(tau2g)
      )
    ),
    min(
      2.0,
      truncnorm::rtruncnorm(1,
        a = 0.0, b = Inf, mean = 0.0,
        sd = sqrt(tau2g)
      )
    )
  )

  invisible(list(
    mug = rnorm(1, hp$muG0, hp$sigmaG0),
    mucond = mucond,
    muind = muind,
    kappa2g = kappa2g,
    tau2g = tau2g,
    phi2g = min(10.0, rinvgamma(1,
      shape = hp$alphaPhi2G,
      scale = hp$betaPhi2G
    ))
  ))
}

generate_gwnc_y <- function(mug, mucond, muind, phi2g,
                            nind, ncell) {
  ## generate data based on the gwnb model
  ## given the model parameters
  ## need muind

  n <- nind * ncell * 2
  sample_sumcnt <- sample(sumcnt, size = n, replace = TRUE)

  ## cond in simulation start from 1 following stan.
  vec_of_cond <- c(rep(1, nind * ncell), rep(2, nind * ncell))

  index_of_ind_per_cond <- get_vec_of_repeat_int(nind, ncell)
  vec_of_ind_under_cond <- c(
    paste0("NR", index_of_ind_per_cond),
    paste0("R", index_of_ind_per_cond)
  )
  vec_of_ind <- get_vec_of_repeat_int(nind * 2, ncell)
  y <- vapply(
    seq_len(n),
    function(i) {
      logmu <- mug + mucond[vec_of_cond[i]] + muind[vec_of_ind[i]]
      invisible(rnbinom(1,
        mu = sample_sumcnt[i] * exp(logmu), size = phi2g
      ))
    }, 0.0
  )
  invisible(list(
    n = n,
    vec_of_cond = vec_of_cond,
    vec_ind_under_cond = vec_of_ind_under_cond,
    vec_of_ind = vec_of_ind,
    s = sample_sumcnt,
    y = y
  ))
}

set_init_params <- function(simu_data, hp, nind, default_control_value = -1.0,
                            default_case_value = 1.0) {
  sy <- simu_data$y
  ss <- simu_data$s
  resp <- simu_data$vec_of_cond

  ## cond in simulation start from 1 following stan.
  index_of_control <- (resp == 1)
  index_of_case <- (resp != 1)
  fit_mur <- myfit$stanfit_gwsnb_till_cond_level(
    s = ss,
    y = sy,
    s_control = ss[index_of_control],
    y_control = sy[index_of_control],
    s_case = ss[index_of_case],
    y_case = sy[index_of_case],
    scale_nb_model = snbm,
    scale_nb_fixed_r_model = snb_fixr_m
  )

  if (is.nan(fit_mur$mu_cond[1])) {
    fit_mur$mucond[1] <- default_control_value
  }

  if (is.nan(fit_mur$mu_cond[2])) {
    fit_mur$mucond[2] <- default_case_value
  }

  ## we use the max absolute value, and cover it by 1.5 fold
  tau2g <- max(abs(fit_mur$mu_cond)) * 1.5
  invisible(list(
    MuG = fit_mur$mu0,
    MuIndRaw = rep(0.0, nind * 2),
    MuCond = fit_mur$mu_cond,
    MuCondRaw = fit_mur$mu_cond / sqrt(tau2g),
    Kappa2G = 1.0,
    Tau2G = tau2g,
    Phi2G = fit_mur$r0
  ))
}

set_gwnb_light_env <- function(nind) {
  ## set hyper params and sampling params from the prior

  ## global variables:
  ## use pbmc data with the specific cell type, and genenm.

  ## use pbmc to estimate and set the part of the parameters.
  fit <- fit_singlegene_nb()
  ## simulation from prior
  hp <- set_gwnb_hyper_params(fit, sigmaG0 = 4.0)
  params <- simulate_from_gwnb_prior(hp, nind)
  invisible(list(
    sg_nbfit = fit,
    hp = hp,
    params_from_prior = params
  ))
}

run_gwnb_model <- function(model, model_env,
                           data, method = "vi",
                           use_init = TRUE) {
  ## run the model using VI or MAP with/without init params
  init_params <- NULL
  if (use_init) {
    t <- set_init_params(data, model_env$hp)
    init_params <- list(t)
  }
  if (method == "vi") {
    fit <- model$variational(
      data = data, init = init_params,
      seed = 1, refresh = 5000,
      iter = 5000, eval_elbo = 100,
      adapt_engaged = TRUE
    )
  } else {
    fit <- model$optimize(
      data = data, init = init_params,
      seed = 1, refresh = 5000,
      iter = 5000, algorithm = "lbfgs"
    )
  }
  invisible(fit)
}


hist_vi_opt <- function(vi, opt, gwnb_env, varnm,
                        lsize = 1.2 * c(1.5, 0.7, 1.2),
                        lalpha = c(0.75, 0.75, 0.75),
                        lcolor = c("red", "gray20", "blue"),
                        ltype = c(1, 2, 2), show_opt = TRUE,
                        show_init = TRUE) {
  ## show lines in order: ground truth; init value; MAP

  if (varnm == "MuG") {
    ground_truth <- gwnb_env$params_from_prior$mug
    init_value <- gwnb_env$init_params[[varnm]]
  }
  if (varnm == "MuCond[1]") {
    ground_truth <- gwnb_env$params_from_prior$mucond[1]
    init_value <- gwnb_env$init_params[["MuCond"]][1]
  }
  if (varnm == "MuCond[2]") {
    ground_truth <- gwnb_env$params_from_prior$mucond[2]
    init_value <- gwnb_env$init_params[["MuCond"]][2]
  }
  if (varnm == "Kappa2G") {
    ground_truth <- gwnb_env$params_from_prior$kappa2g
    init_value <- gwnb_env$init_params[[varnm]]
  }
  if (varnm == "Tau2G") {
    ground_truth <- gwnb_env$params_from_prior$tau2g
    init_value <- gwnb_env$init_params[[varnm]]
  }
  if (varnm == "Phi2G") {
    ground_truth <- gwnb_env$params_from_prior$phi2g
    init_value <- gwnb_env$init_params[[varnm]]
  }

  ## add hist graph with two lines:
  ## simulation truth and designed init params.
  p <- bayesplot::mcmc_hist(vi$draws(varnm)) +
    bayesplot::vline_at(ground_truth,
      size = lsize[1], alpha = lalpha[1],
      color = lcolor[1], linetype = ltype[1]
    )
  if (show_init) {
    p <- p +
      bayesplot::vline_at(init_value,
        size = lsize[2], alpha = lalpha[2],
        color = lcolor[2], linetype = ltype[2]
      )
  }

  if (show_opt) {
    p <- p +
      bayesplot::vline_at(opt$mle(varnm),
        size = lsize[3],
        alpha = lalpha[3], color = lcolor[3], linetype = ltype[3]
      )
  }
  ## p <- p + bayesplot::facet_text(size = 15, hjust = 0.5)
  invisible(p)
}

visualize_vi_init <- function(vi, noi_vi, opt, noi_opt, gwnb_env,
                              varnms = c(
                                "MuG", "MuCond[1]", "MuCond[2]",
                                "Kappa2G", "Tau2G", "Phi2G"
                              )) {
  ## generate a list named with varnms
  ## and each is a list with two element in order:
  ## noivip, vip

  l <- lapply(varnms, function(varnm) {
    p <- list(
      noivip = hist_vi_opt(noi_vi, noi_opt, gwnb_env,
        varnm,
        show_opt = TRUE, show_init = FALSE
      ),
      vip = hist_vi_opt(vi, opt, gwnb_env,
        varnm,
        show_opt = TRUE,
        show_init = TRUE
      )
    )
    invisible(p)
  })
  names(l) <- varnms
  invisible(l)
}


check_model <- function(tag = 1, model = gwnb_model,
                        ninds = ninds, ncells = ncells) {
  ## generat figures: a list of ninds
  ## each is a list of ncells
  ## for each ncell, its also a list of figures for different varnames.
  ## for each varnames, it's a list with two element in order: noivip vip

  l_of_ind <- lapply(ninds, function(nind) {
    ## using simulated data to check vi
    model_env <- set_gwnb_light_env()
    ## ground truth
    gt <- model_env$params_from_prior
    l_of_cell <- lapply(ncells, function(ncell) {
      data <- generate_gwnc_y(
        mug = gt$mug,
        mucond = gt$mucond,
        muind = gt$muind,
        phi2g = gt$phi2g,
        nind = nind,
        ncell = ncell
      )
      vifit <- run_gwnb_model(model,
        model_env,
        data,
        method = "vi",
        use_init = TRUE
      )
      noi_vifit <- run_gwnb_model(model,
        model_env,
        data,
        method = "vi",
        use_init = FALSE
      )
      optfit <- run_gwnb_model(model,
        model_env,
        data,
        method = "opt",
        use_init = TRUE
      )
      noi_optfit <- run_gwnb_model(model,
        model_env,
        data,
        method = "opt",
        use_init = FALSE
      )
      invisible(hist_vi_opt(
        vifit, noi_vifit,
        optfit, noi_optfit, model_env
      ))
    })
    names(l_of_cell) <- ncells
  })
  names(l_of_ind) <- ninds
  invisible(l_of_ind)
}

## for debug
figures <- check_model(1, gwnb_model, ninds, ncells)

## main
## lapply(seq_len(3), function(i) {
##   result <- tryCatch(
##     {
##       figures <- check_model(i, gwnb_model,ninds, ncells)
##       0L
##     },
##     error = function(e) {
##       message(e)
##       1L
##     }
##   )
## })
