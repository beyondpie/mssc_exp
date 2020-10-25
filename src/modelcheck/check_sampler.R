# Following Jun's suggestions, here we check
# the sampler firstly by generating from the prior.

# Set R environment

## Following Jun's suggestions
## Eval the sampler firstly by generating from the prior.

import::from(here, here)
suppressPackageStartupMessages(library(tidyverse))
library(MCMCpack)
library(cmdstanr)
## detach("package:cmdstanr", unload = TRUE)
## library(cmdstanr)
library(bayesplot)
library(posterior)
library(bbmle)
library(sads)

## local modules
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
## modules::reload(myfit)
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## warnings/errors traceback settings
options(error = traceback)
options(warn = 0)
options(mc.cores = 3)
set.seed(355113)

# Setting functions
get_fitted_mu <- function(y, vec_of_cond,
                          vec_of_ind_under_cond,
                          vec_of_s,
                          num_of_ind_per_cond = 5,
                          label_of_control = 0,
                          nm_of_control = "NR",
                          nm_of_case = "R") {
  ## fit NB distribution, then
  ## extract the different means for MSSC model.

  ## mean of individuals will be 5 control; 5 case in order.

  y_fit <- myfit$prob_zero_nb(y, F)
  mu0 <- log(y_fit$nbfit$estimate["mu"] / median(vec_of_s))

  index_of_cond <- vec_of_cond == label_of_control
  y_control_fit <- myfit$prob_zero_nb(y[index_of_cond], F)
  mu_control <- log(y_control_fit$nbfit$estimate["mu"] /
    median(vec_of_s[index_of_cond])) - mu0
  mu_control_ind <- vapply(
    paste0(nm_of_control, seq_len(num_of_ind_per_cond)),
    FUN = function(x) {
      index <- vec_of_ind_under_cond == x
      fit <- myfit$prob_zero_nb(y[index], F)
      invisible(log(fit$nbfit$estimate["mu"] /
        median(vec_of_s[index])) - mu_control - mu0)
    },
    FUN.VALUE = 0.0
  )

  index_of_case <- vec_of_cond != label_of_control
  y_case_fit <- myfit$prob_zero_nb(y[index_of_case], F)
  mu_case <- log(y_case_fit$nbfit$estimate["mu"] /
    median(vec_of_s[index_of_case])) - mu0
  mu_case_ind <- vapply(
    paste0(nm_of_case, seq_len(num_of_ind_per_cond)),
    FUN = function(x) {
      index <- vec_of_ind_under_cond == x
      fit <- myfit$prob_zero_nb(y[index], F)
      invisible(log(fit$nbfit$estimate["mu"] /
        median(vec_of_s[index])) - mu_case - mu0)
    },
    FUN.VALUE = 0.0
  )
  invisible(list(mu0 = mu0,
    mu_cond = c(mu_control, mu_case),
    mu_ind = c(mu_control_ind, mu_case_ind)
  ))
}

## * load pbmc for parameter estimate and setting.
## Classical genes as DE
##   - SNHG16, OASL, NAMPT, NFKB1, BCL2L11, TRAF4, ICAM1, XCL2, XCL1, CCL3L3, CCL3L1
## Strong individual effect genes
##   - HBA1, HBA2, HBD
## possible non-DE
##   - TOX, YIPF5, CCL3, KDM6A, HDDC2

cell_type <- "Naive CD4+ T"
num_of_cell_per_ind <- 280
num_of_ind <- 10
num_of_ind_per_cond <- 5
num_of_cond <- 2
# a DE gene: DE one gene
d1g <- "NFKB1"

## the whole dataset
pbmc_seurat <- mypbmc$load_pbmc_seurat() %>%
  mypbmc$extract_from_seurat(pbmc_seurat = .)
## limit to the cell type
subscdata <- mypbmc$get_celltype_specific_scdata(
  pbmc_seurat$cnt,
  pbmc_seurat$resp,
  pbmc_seurat$inds,
  pbmc_seurat$ct,
  cell_type
)

## ** sample cell numbers
sample_cells <- mypbmc$sample_cells_per_ind(
  subscdata$inds,
  num_of_cell_per_ind
)
## note individual order is changed according to sample_cells
cnt <- subscdata$cnt[, sample_cells]
inds <- subscdata$inds[sample_cells]
resp <- subscdata$resp[sample_cells]
colsumcnt <- colSums(cnt)

## ** limit to one DE gene
d1g_cnt <- cnt[d1g, ]
outliers <- myfit$is_outlier(d1g_cnt)
d1g_cnt <- d1g_cnt[!outliers]

par(mfrow = c(1, 2))
hist(d1g_cnt[resp == 1])
hist(d1g_cnt[resp == 0])

## remove outliers
d1g_inds <- inds[!outliers]
d1g_resp <- resp[!outliers]
d1g_sumcnt <- colsumcnt[!outliers]

## fit NB dist
nb_fit_all <- myfit$prob_zero_nb(d1g_cnt, F)
nb_fit_control <- myfit$prob_zero_nb(d1g_cnt[d1g_resp == 0], F)
nb_fit_case <- myfit$prob_zero_nb(d1g_cnt[d1g_resp == 1], F)

nb_fit_case_ind_mu <- vapply(
  paste0("R", seq_len(5)),
  FUN = function(x) {
    fit <- myfit$prob_zero_nb(d1g_cnt[d1g_inds == x])
    invisible(log(fit$nbfit$estimate["mu"] / median(d1g_sumcnt[d1g_inds == x])))
  },
  FUN.VALUE = 0.0
)

nb_fit_control_ind_mu <- vapply(
  paste0("NR", seq_len(5)),
  FUN = function(x) {
    fit <- myfit$prob_zero_nb(d1g_cnt[d1g_inds == x])
    invisible(log(fit$nbfit$estimate["mu"] / median(d1g_sumcnt[d1g_inds == x])))
  }, FUN.VALUE = 0.0
)

## -------
## fit NB for a list of genes, add see the distribution
## of the dispersion.
## -------

## estimate the mu in the model.
mu0 <- log(nb_fit_all$nbfit$estimate["mu"] / median(d1g_sumcnt))
mu_control <- log(nb_fit_control$nbfit$estimate["mu"] /
  median(d1g_sumcnt[d1g_resp == 0])) - mu0
mu_control_ind <- nb_fit_control_ind_mu -
  log(nb_fit_control$nbfit$estimate["mu"] / median(d1g_sumcnt[d1g_resp == 0]))
mu_case <- log(nb_fit_case$nbfit$estimate["mu"] /
  median(d1g_sumcnt[d1g_resp == 1])) - mu0
mu_case_ind <- nb_fit_case_ind_mu -
  log(nb_fit_case$nbfit$estimate["mu"] / median(d1g_sumcnt[d1g_resp == 1]))


# set data for check model from prior
Ind <- c(vapply(seq_len(num_of_ind),
  FUN = function(x) {
    rep(x, num_of_cell_per_ind)
  },
  FUN.VALUE = c(rep(0.0, num_of_cell_per_ind))
))

## * load stan model
sbc_gwnb_model <- cmdstan_model(here(
  "src", "dirty_stan",
  "gwnb_simu_from_prior.stan"
), compile = T, force_recompile = T)

## * set gwnb hyper prior
## for muG
muG0 <- mu0
sigmaG0 <- 4.0

## for MuInd
alphaKappa2G <- 1.0
betaKappa2G <- 1.0

## for MuCond
alphaTau2G <- 1.0
betaTau2G <- 1.0

## for NB
alphaPhi2G <- 1.0
betaPhi2G <- 1.0

## * sampling parameters from the prior
mug <- rnorm(1, mean = muG0, sd = sigmaG0)
kappa2g <- rinvgamma(1,
  shape = alphaKappa2G,
  scale = betaKappa2G
)
tau2g <- rinvgamma(1,
  shape = alphaTau2G,
  scale = betaTau2G
)
phi2g <- rinvgamma(1,
  shape = alphaPhi2G,
  scale = betaPhi2G
)

muind <- rnorm(num_of_ind,
  mean = 0.0,
  sd = sqrt(kappa2g)
)
mucond <- rnorm(num_of_cond,
  mean = 0.0,
  sd = sqrt(tau2g)
)
y <- vapply(seq_len(num_of_cell_per_ind * num_of_ind),
  FUN = function(i) {
    rnbinom(
      n = 1,
      mu = d1g_sumcnt[i] * exp(mug + muind[Ind[i]] + mucond[Cond[i]]),
      size = phi2g
    )
  },
  FUN.VALUE = 0.0
)

## set the initial values
vec_of_cond <- c(rep(1, num_of_cell_per_ind * num_of_ind_per_cond),
  rep(2, num_of_cell_per_ind * num_of_ind_per_cond))
vec_of_ind_under_cond <- c(paste0("NR",
                                  rep(seq_len(num_of_ind_per_cond),
                                      num_of_cell_per_ind)),
  paste0("R", rep(seq_len(num_of_ind_per_cond), num_of_cell_per_ind)))

sbc_mu_list <- get_fitted_mu(y, vec_of_cond, vec_of_ind_under_cond,
  vec_of_s = d1g_sumcnt,
  label_of_control = 1)

kappa2g <- 1.0
tau2g <- 1.0
phi2g <- 1.0

init_params <- list(
  MuG = sbc_mu_list$mu0,
  MuIndRaw = sbc_mu_list$mu_ind / sqrt(kapp2g),
  MuCondRaw = sbc_mu_list$mu_cond / sqrt(tau2g),
  Kappa2G = kapp2g,
  Tau2G = tau2g,
  Phi2G = phi2g
)

## set the data
vec_of_ind <- c(vapply(seq_len(num_of_ind),
  FUN = function(x) {
    rep(x, num_of_cell_per_ind)
  },
  FUN.VALUE = c(rep(0.0, num_of_cell_per_ind))
))

mydata <- list(
  N = num_of_cell_per_ind * num_of_ind,
  K = num_of_ind,
  J = num_of_cond,
  S = d1g_sumcnt,
  Ind = vec_of_ind,
  Cond = vec_of_cond,
  kappa2g = kappa2g,
  tau2g = tau2g,
  phi2g = phi2g,
  mug = mug,
  muind = muind,
  mucond = mucond,
  y = y,
  muG0 = muG0,
  sigmaG0 = sigmaG0,
  alphaKappa2G = alphaKappa2G,
  betaKappa2G = betaKappa2G,
  alphaTau2G = alphaTau2G,
  betaTau2G = betaTau2G,
  alphaPhi2G = alphaPhi2G,
  betaPhi2G = betaPhi2G
)

## ** analyze the variational inference
sbc_gwnb_vi_sampler <- sbc_gwnb_model$variational(
  data = mydata,
  seed = 1,
  init = list(init_params)
)
sbc_gwnb_vi_draws <- sbc_gwnb_vi_sampler$draws()
