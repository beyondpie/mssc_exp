## Following Jun's suggestions
## Eval the sampler firstly by generating from the prior.

import::from(here, here)
suppressPackageStartupMessages(library(tidyverse))
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

## ## set cmdstan path
set_cmdstan_path(path = paste(Sys.getenv("HOME"),
  "softwares",
  "cmdstan-2.23.0",
  sep = "/"
))

## * load pbmc for parameter estimate and setting.
## classical genes as DE
## SNHG16, OASL, NAMPT, NFKB1, BCL2L11, TRAF4, ICAM1,
## XCL2, XCL1, CCL3L3, CCL3L1

## strong individual effect genes
## HBA1, HBA2, HBD

## possible non-DE
## TOX, YIPF5, CCL3, KDM6A, HDDC2

cell_type <- "Naive CD4+ T"
sample_cellnum <- 280
num_of_ind <- 10
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
  sample_cellnum
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


## * set data for eval model from prior
Ind <- c(vapply(seq_len(num_of_ind),
  FUN = function(x) {
    rep(x, sample_cellnum)
  },
  FUN.VALUE = c(rep(0.0, sample_cellnum))
))
Cond <- c(vapply(seq_len(num_of_cond),
  FUN = function(x) {
    rep(x, sample_cellnum * num_of_ind / 2)
  },
  FUN.VALUE = c(rep(0.0, sample_cellnum * num_of_ind / 2))
))

mydata <- list(
  N = sample_cellnum * num_of_ind,
  K = num_of_ind,
  J = num_of_cond,
  S = d1g_sumcnt,
  Ind = Ind,
  Cond = Cond
)
## * analyze simulation-based calibration
sbc_gwnb_model <- cmdstan_model(here(
    "src", "dirty_stan",
    "gwnb_simu_from_prior.stan"
),
  compile = T,
  force_recompile = T
)


















