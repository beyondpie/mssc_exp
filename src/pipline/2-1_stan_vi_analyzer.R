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

## * load genes information
deg <- readRDS(here(data_dir, subdir, de_outfnm))
fpdeg <- readRDS(here(data_dir, subdir, fpde_outfnm))
tndeg <- readRDS(here(data_dir, subdir, tnde_outfnm))
sc_data_list <- readRDS(here(data_dir, subdir, scdata_sumfnm))
sc_genes <- rownames(sc_data_list$cnt)

## * retrieve de/nonde-related genes
degnms <- myt$stat_geneset(sc_genes, deg$genesymbol)
fpdegnms <- myt$stat_geneset(sc_genes, fpdeg$genesymbol)
tndegnms <- myt$stat_geneset(sc_genes, tndeg$genesymbol)

## TODO: add explanation about which is case and control
## ctrlmnscase: control minus case
get_ctrlmnscase_par <- function(mystanfit, par = "MuCond") {
  mus <- rstan::extract(mystanfit, pars = par)[[par]]
  delta <- as.data.frame(mus[, , 1] - mus[, , 2])
  colnames(delta) <- sc_genes
  return(delta)
}

## simple t statistics
calt <- function(delta, fn=matrixStats::colMedians) {
  fnhat <- fn(delta)
  std_hat <- matrixStats::colSds(delta + 1e-10)
  return(fnhat / (sqrt(nrow(delta) * std_hat)))
}

## AUC analysis
calauc <- function(scores, backend) {
  caTools::colAUC(scores, backend)
}

## eval scores based on posterior samples using AUC
evalstat <- function(modelnm = "v1-1", method = "vi", par = "MuCond",
                     fn = matrixStats::colMeans,
                     degnms=degnms, ndegnms=tndegnms) {
  mystanfit <- myt$load_stan(
    here(exp_dir, exp_sub_dir, stan_dir),
    modelnm, method
    )
  dmucond <- get_ctrlmnscase_par(mystanfit = mystanfit, par = par)
  dmut <- calt(dmucond, fn)
  mybackend <- c(rep(TRUE, length(degnms)), rep(FALSE, length(ndegnms)))
  bgnms <- c(degnms, ndegnms)
  names(mybackend) <- bgnms
  myauc <- calauc(dmut[names], mybackend)
  message(str_glue("model {modelnm} with method {method}"))
  message(str_glue("parameter: {par}"))
  message(str_glue("AUC: {myauc}"))
  return(myauc)
}

## point relative to regions
mypntrela2rgn <- function(myintvals, myprobs = c(0.025, 0.975), pnt = 0.0) {
  q <- as.data.frame(apply(myintvals, 2, quantile, probs = myprobs))
  apply(q, 2, myt$pntrela2rgn, pnt = pnt)
}

## Bayesin posterial quatile evaluation
evalqtl <- function(modelnm = "v1-1", method = "vi", par = "MuCond",
                    myprobs = c(0.025, 0.975), pnt = 0.0,
                    deg = degnms, ndeg = tndegnms) {
  mystanfit <- myt$load_stan(here(exp_dir, exp_sub_dir, stan_dir),
                             modelnm, method)
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

  tpr <- myt$fmtflt(myroc$tpr(tp, fp, tn, fn))
  fpr <- myt$fmtflt(myroc$fpr(tp, fp, tn, fn))
  fdr <- myt$fmtflt(myroc$fdr(tp, fp, tn, fn))
  f1 <- myt$fmtflt(myroc$f1(tp, fp, tn, fn))

  message(str_glue("model {modelnm} with method {method}"))
  message(str_glue("parameter: {par}"))
  message(str_glue("Bayesian quantile: ({myprobs[1]}, {myprobs[2]})"))
  message(str_glue("TPR({tpr}), FPR({fpr}), FDR({fdr}). F1({f1})"))
  message(str_glue("TP({tp}), FP({fp}), TN({tn}), FN({fn})"))
}

## * main
## ** performance analyze
## *** quantile based
evalqtl("v1-1", "vi", myprobs = c(0.01, 0.99))
evalqtl("v1-1", "vi", myprobs = c(0.025, 0.975))
evalqtl("v1-1", "vi", myprobs = c(0.05, 0.95))

## *** statisitc rank based, i.e., AUC
message("Using mean divided by std with true negative genes")
evalstat("v1-1", "vi", par = "MuCond", fn = colMeans,
         degnms = degnms, ndegnms = tndegnms)
