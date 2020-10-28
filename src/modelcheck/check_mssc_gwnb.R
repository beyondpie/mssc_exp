## Check MSSC variational inference model.

## Use simulation based calibration.
## - simulation from prior.
## - [hyper] parameters are setted based on a real dataset:
##   - PBMC dataset

## * set R environment

import::from(here, here)
suppressPackageStartupMessages(library(tidyverse))
library(MCMCpack)
library(cmdstanr)

## use the code below to reload the library
## detach("package:cmdstanr", unload = TRUE)
## library(cmdstanr)

library(bayesplot)
## for bayesplot plotting
## color_scheme_set("brewer-Spectral")

library(posterior)
library(bbmle)
library(sads)

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

## * load stan models
mssc_gwnb_muind_model <- cmdstan_model(here::here(
  "src", "dirty_stan",
  "gwnb_simu_from_prior.stan"
), compile = T)

mssc_gwnb_rndeff_model <- cmdstan_model(
  here::here("src", "dirty_stan", "gwnb_rndeffect.stan"),
  compile = T)

## * configs
## for simulation
num_of_cell_per_ind <- 280
num_of_ind <- 10
num_of_ind_per_cond <- 5
num_of_cond <- 2

## for the PBMC dataset
cell_type <- "Naive CD4+ T"
sgn <- "NFKB1"


## functions
fit_singlegene_nb <- function(gn = "NFKB1", cnt = cnt,
                              inds = inds,
                              resp = resp,
                              sumcnt = sumcnt) {
  ## get mu0, mucond, r as in mssc gwnb model
  ## from the real dataset

  d1g_cnt <- cnt[gn, ]
  outliers <- myfit$is_outlier(d1g_cnt)
  d1g_cnt <- d1g_cnt[!outliers]
  d1g_inds <- inds[!outliers]
  d1g_resp <- resp[!outliers]
  d1g_sumcnt <- sumcnt[!outliers]

  ## fit NB dist
  fit_mur <- myfit$fit_gwnb_s2_till_cond_level(
    y = d1g_cnt,
    y_control = d1g_cnt[d1g_resp == 0],
    y_case = d1g_cnt[d1g_resp == 1],
    scale_of_s = median(d1g_sumcnt)
  )
  invisible(list(
    mur = fit_mur,
    cnt = d1g_cnt,
    inds = d1g_inds,
    resp = d1g_resp,
    sumcnt = d1g_sumcnt
  ))
}

set_gwnb_hyper_params <- function(sg_mur, sigmaG0 = 4.0) {
  ## set hyper params based on the fitted sg_mur
  ## from a real dataset

  ## inv-gamma mode estimate.
  ## instead of using var(sg_mur$mu_cond), which might be
  ## small when both of them are in the same sign:
  ## - we assume the mean is around zero

  ## we fix the alpha as 1.0
  alphaTau2G <- 1.0
  ## we use the max absolute value, and cover it by 1.5 fold
  betaTau2G <- (alphaTau2G + 1) * max(abs(sg_mur$mu_cond)) * 1.5
  alphaPhi2G <- 1.0
  betaPhi2G <- min(5.0, (alphaPhi2G + 1) * sg_mur$r0)

  invisible(list(alphaKappa2G = 1.0,
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
  invisible(c(vapply(seq(from_int, to_int),
    function(i) {
      rep(i, repeatimes)
    }, rep(0L, repeatimes))
  ))
}

generate_gwnc_y <- function(mug, mucond, kappa2g, vec_sumcnt) {
  ## geneerate data based on the gwnb model
  ## given the parameters.

  half_n <- num_of_ind_per_cond * num_of_cell_per_ind
  n <- 2 * half_n
  vec_of_cond <- c(rep(1, half_n), rep(2, half_n))
  index_of_ind_per_cond <- get_vec_of_repeat_int(num_of_ind_per_cond,
    num_of_cell_per_ind)
  vec_of_ind_under_cond <- c(paste0("NR", index_of_ind_per_cond),
    paste0("R", index_of_indi_per_cond))
  vec_of_ind <- get_vec_of_repeat_int(num_of_ind, num_of_cell_per_ind)
  y <- vapply(seq_len(n),
    function(i) {
      logmu <- rnorm(1,
        mean = mug + mucond[vec_of_cond[i]],
        sd = sqrt(kappa2g)
      )
      rnbinom(
        1, mu = vec_sumcnt[i] * exp(logmu), size = phi2g
      )
    }, 0L)

  invisible(n = n,
    vec_of_cond = vec_of_cond,
    vec_ind_under_cond = vec_of_ind_under_cond,
    vec_of_ind = vec_of_ind,
    s = vec_sumcnt,
    y = y)
}

simulate_from_gwnb_prior <- function(hp, num_of_cond = 2) {
  ## sampling the parameters from the prior
  ## defined by the hyper parameters.
  tau2g <- rinvgamma(1, shape = hp$alphaTau2G,
    scale = hp$betaTau2G)
  mucond <- rnorm(num_of_cond, 0.0, sqrt(tau2g))

  invisible(list(mug = max(-7.0,
    rnorm(1, hp$muG0, hp$sigmaG0)),
  kappa2g = rinvgamma(1, shape = hp$alphaKappa2G,
    scale = hp$betaKappa2G),
  tau2g = taug2g,
  phi2g = rinvgamma(1, shape = hp$alphaPhi2G,
    scale = hp$betaPhi2G),
  mucond = mucond
  ))
}

set_init_params <- function(simu_data, hp) {
  sy <- simu_data$y
  kappa2g <- 1.0
  fit_mur <- myfit$fit_gwnb_s2_till_cond_level(
    y = sy,
    y_control = sy[simu_data$vec_of_cond == 1],
    y_case = sy[simu_data$vec_of_cond != 1],
    scale_of_s = median(simu_data$s)
  )
  ## we use the max absolute value, and cover it by 1.5 fold
  tau2g <- max(abs(fit_mur$mu_cond)) * 1.5
  invisible(list(
    MuG = fit_mur$mu0,
    MuIndRaw = rnom(simu_data$n, 0.0, sqrt(kappa2g)),
    MuCondRaw = fit_mur$mu_cond / sqrt(tau2g),
    Kappa2G = 1.0,
    Tau2G = tau2g,
    Phi2G = fit_mur$r0
  ))
}

## ** load the dataset
pbmc_seurat <- mypbmc$load_pbmc_seurat() %>%
  mypbmc$extract_from_seurat(pbmc_seurat = .)

subscdata <- mypbmc$get_celltype_specific_scdata(
  pbmc_seurat$cnt,
  pbmc_seurat$resp,
  pbmc_seurat$inds,
  pbmc_seurat$ct,
  ## limit to the cell type
  cell_type
)

## sample cell numbers
sample_cells <- mypbmc$sample_cells_per_ind(
  subscdata$inds,
  num_of_cell_per_ind
)
## note individual order is changed according to sample_cells
cnt <- subscdata$cnt[, sample_cells]
inds <- subscdata$inds[sample_cells]
resp <- subscdata$resp[sample_cells]
sumcnt <- colSums(cnt)

## ** use pbmc to estimate and set the parameters.
pbmc_sg_fit <- fit_singlegene_nb(gn = sgn)

## summary of the fit
str(pbmc_sg_fit)

par(mfrow = c(1, 2))
hist(pbmc_sg_fit$cnt[pbmc_sg_fit$resp == 0])
hist(pbmc_sg_fit$cnt[pbmc_sg_fit$resp == 1])

## * simulation from prior
gwnb_hyperparmas <- set_gwnb_hyper_params(pbmc_sg_fit, sigmaG0 = 4.0)
params_from_prior <- simulate_from_gwnb_prior(gwnb_hyperparmas)
gwnb_genrt_data <- generate_gwnc_y(params_from_prior$mug,
  params_from_prior$mucond,
  params_from_prior$kappa2g,
  vec_sumcnt = sumcnt)

## set data for mssc.
input_data <- c(list(
    N = gwnb_genrt_data$n,
    K = num_of_ind,
    J = num_of_cond,
    S = gwnb_genrt_data$s,
    Ind = gwnb_genrt_data$vec_of_ind,
    Cond = gwnb_genrt_data$vec_of_cond,
    y = gwnb_genrt_data$y),
  gwnb_hyperparmas)

## * set initial params
init_params <- set_init_params(gwnb_genrt_data,
  gwnb_hyperparmas)

## * run models
## ** hmc
## ** optmize
## ** vi

