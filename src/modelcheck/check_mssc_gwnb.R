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


## * load the dataset
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

## * use pbmc to estimate and set the parameters.
pbmc_sg_fit <- fit_singlegene_nb(gn = sgn)

## summary of the fit
str(pbmc_sg_fit)

par(mfrow = c(1, 2))
hist(pbmc_sg_fit$cnt[pbmc_sg_fit$resp == 0])
hist(pbmc_sg_fit$cnt[pbmc_sg_fit$resp == 1])




