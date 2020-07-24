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
exp_dir <- "data"
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
calt <- function(delta, fn = matrixStats::colMedians) {
  fnhat <- fn(as.matrix(delta))
  std_hat <- matrixStats::colSds(as.matrix(delta) + 1e-10)
  sts <- fnhat / (sqrt(nrow(delta)) * std_hat)
  names(sts) <- colnames(delta)
  return(sts)
}

## AUC analysis
calauc <- function(scores, backend) {
  caTools::colAUC(scores, backend)
}

getndegnms <- function(myndegnm = "extreme") {
  if (myndegnm == "extreme") {
    ndegnms <- tndegnms
  }
  if (myndegnm == "nearpositive") {
    ndegnms <- fpdegnms
  }
  if (myndegnm == "all") {
    ndegnms <- c(tndegnms, fpdegnms)
  }
  message(str_glue("using {myndegnm} negatives"))
  return(ndegnms)
}

## eval scores based on posterior samples using AUC
evalstat <- function(modelnm = "v1-1", method = "vi", par = "MuCond",
                     fnm = "mean", myndegnm = "extreme",
                     mydegnms = degnms) {
  mystanfit <- myt$load_stan(
    here(exp_dir, exp_sub_dir, stan_dir),
    modelnm, method
  )
  dmucond <- get_ctrlmnscase_par(mystanfit = mystanfit, par = par)
  if (fnm == "mean") {
    fn <- colMeans
  }
  if (fnm == "median") {
    fn <- matrixStats::colMedians
  }
  dmut <- calt(dmucond, fn)

  ndegnms <- getndegnms(myndegnm)
  mybackend <- c(rep(TRUE, length(mydegnms)), rep(FALSE, length(ndegnms)))
  bgnms <- c(mydegnms, ndegnms)
  names(mybackend) <- bgnms
  myauc <- myt$fmtflt(calauc(dmut[bgnms], mybackend))
  message(str_glue("model {modelnm} with method {method}"))
  message(str_glue(
    "parameter: {par} with stats {fnm}"
  ))
  message(str_glue("AUC: {myauc}"))
  return(list(auc = myauc, sts = dmut))
}

## point relative to regions
mypntrela2rgn <- function(myintvals, myprobs = c(0.025, 0.975), pnt = 0.0) {
  q <- as.data.frame(apply(myintvals, 2, quantile, probs = myprobs))
  apply(q, 2, myt$pntrela2rgn, pnt = pnt)
}

## Bayesin posterial quatile evaluation
evalqtl <- function(modelnm = "v1-1", method = "mc", par = "MuCond",
                    myprobs = c(0.025, 0.975), pnt = 0.0,
                    deg = degnms, myndegnm = "extreme") {
  mystanfit <- myt$load_stan(
    here(exp_dir, exp_sub_dir, stan_dir),
    modelnm, method
  )
  ## dmucond: delta mucond
  dmucond <- get_ctrlmnscase_par(mystanfit = mystanfit, par = par)
  zerorela2dmucond <- mypntrela2rgn(dmucond, myprobs = myprobs, pnt = pnt)

  ndegnms <- getndegnms(myndegnm)
  pred_deg <- sc_genes[zerorela2dmucond != 0]
  pred_ndeg <- setdiff(sc_genes, pred_deg)

  pred_tpg <- intersect(pred_deg, deg)
  pred_tng <- intersect(pred_ndeg, ndegnms)
  pred_fpg <- intersect(pred_deg, ndegnms)
  pred_fng <- intersect(pred_ndeg, deg)

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
  message(str_glue("TPR({tpr}), FPR({fpr}), FDR({fdr}), F1({f1})"))
  message(str_glue("TP({tp}), FP({fp}), TN({tn}), FN({fn})"))
}

## * main
## ** performance analyze
## *** quantile based
## region to use (0.01, 0.99), (0.025, 0.975), (0.05, 0.95)
for (myndeg in c("extreme", "nearpositive", "all")) {
  for (mymethod in c("mc")) {
    tmp <- evalqtl("v1-1", mymethod,
      myprobs = c(0.25, 0.75), myndegnm = myndeg
    )
  }
}

## *** statisitc rank based, i.e., AUC
for (myndeg in c("extreme", "nearpositive", "all")) {
  for (mymethod in c("vi", "mc")) {
    tmp <- evalstat("v1-1", mymethod,
      par = "MuCond",
      fnm = "mean", myndegnm = myndeg
    )
  }
}
