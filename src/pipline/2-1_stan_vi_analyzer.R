options(error = traceback)
options(warn = -1)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(bayesplot))
suppressPackageStartupMessages(library(data.table))

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(import::from(here, here))
import::from(optparse, make_option, OptionParser, parse_args)
import::from(stringr, str_glue)

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myroc <- modules::import("roc")

## * configs
## ** data related configs
data_dir <- "data"
subdir <- "UM"

de_outfnm <- "tcga_diffexp_genes.rds"
fpde_outfnm <- "tcga_fp_diffexp_genes.rds"
tnde_outfnm <- "tcga_tn_diffexp_genes.rds"
scdata_sumfnm <- "sampled_scRNAseq_summary.rds"

## ** stan result configs
exp_dir <- "exps"
exp_sub_dir <- "UM"
stan_dir <- "stan"
mc_dir <- "mc"
vi_dir <- "vi"

## * load genes information
deg <- readRDS(here(data_dir, subdir, de_outfnm))
fpdeg <- readRDS(here(data_dir, subdir, fpde_outfnm))
tndeg <- readRDS(here(data_dir, subdir, tnde_outfnm))
sc_data_list <- readRDS(here(data_dir, subdir, scdata_sumfnm))
sc_genes <- rownames(sc_data_list$cnt)

## * retrieve de/nonde-related genes
deg <- myt$stat_geneset(sc_genes, deg$genesymbol)
fpdeg <- myt$stat_geneset(sc_genes, fpdeg$genesymbol)
tndeg <- myt$stat_geneset(sc_genes, tndeg$genesymbol)


## * utils for gettign rstan results
load_stan_vi <- function(path) {
  rstan::read_stan_csv(path)
}

load_stan_mc <- function(dirpath, modelnm) {
  csvfiles <- dir(
    path = dirpath, pattern = paste0(modelnm, "[0-9].csv"),
    full.names = T
  )
  rstan::read_stan_csv(csvfiles)
}

## TODO: add explanation about which is case and control
## ctrlmnscase: control minus case
get_ctrlmnscase_par <- function(mystanfit, par = "MuCond") {
  mus <- rstan::extract(mystanfit, pars = par)[[par]]
  delta <- as.data.frame(mus[, , 1] - mus[, , 2])
  colnames(delta) <- sc_genes
  return(delta)
}

## point relative to regions
mypntrela2rgn <- function(myintvals, myprobs = c(0.025, 0.975), pnt = 0.0) {
  q <- as.data.frame(apply(myintvals, 2, quantile, probs = myprobs))
  apply(q, 2, myt$pntrela2rgn, pnt = pnt)
}

## Bayesin posterial quatile evaluation
bayesqtlevl <- function(modelnm = "v1-1", method = "vi", par = "MuCond",
                        myprobs = c(0.025, 0.975), pnt = 0.0) {
  if (method == "vi") {
    vifnm <- paste0(modelnm, ".csv")
    mystanfit <- load_stan_vi(here(
      exp_dir, exp_sub_dir,
      stan_dir, vi_dir, vifnm
    ))
  }
  if (method == "mc") {
    mcprefix <- paste0(modelnm, "_chain_")
    mystanfit <- load_stan_mc(
      dirpath = here(
        exp_dir,
        exp_sub_dir, stan_dir, mc_dir
      ),
      modelnm = mcprefix
    )
  }
  ## dmucond: delta mucond
  dmucond <- get_ctrlmnscase_par(mystanfit = mystanfit, par = par)
  zerorela2dmucond <- mypntrela2rgn(dmucond, myprobs = myprobs, pnt = pnt)
  pred_deg <- sc_genes[zerorela2dmucond != 0]
  pred_ndeg <- setdiff(sc_genes, pred_deg)

  pred_tpg <- intersect(pred_deg, deg)
  pred_tng <- intersect(pred_ndeg, tndeg)
  pred_fpg <- setdiff(pred_deg, deg)
  pred_fng <- setdiff(pred_ndeg, tndeg)

  tp <- length(pred_tpg)
  tn <- length(pred_tng)
  fp <- length(pred_fpg)
  fn <- length(pred_fng)

  tpr <- myroc$tpr(tp, fp, tn, fn)
  fpr <- myroc$fpr(tp, fp, tn, fn)
  fdr <- myroc$fdr(tp, fp, tn, fn)
  f1 <- myroc$f1(tp, fp, tn, fn)
  message(str_glue("model {modelnm} with method {method}"))
  message(str_glue("parameter: {par}"))
  message(str_glue("Bayesian quantile: ({myprobs[1]}, {myprobs[2]})"))
  message(str_glue("TPR({tpr}), FPR({fpr}), FDR({fdr}). F1({f1})"))
}

## * main
## ** performance analyze
bayesqtlevl("v1-1", "vi", myprobs = c(0.01, 0.99))
bayesqtlevl("v1-1", "vi", myprobs = c(0.025, 0.975))
bayesqtlevl("v1-1", "vi", myprobs = c(0.05, 0.95))
