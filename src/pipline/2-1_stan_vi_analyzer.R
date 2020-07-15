options(error = traceback)
options(warn = -1)

library(Seurat)
library(rstan)
library(ggplot2)
library(bayesplot)
library(data.table)
library(pROC)

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
mcprefix <- "v1-1_chain_"
vi_dir <- "vi"
vifnm <- "v1-1.csv"

## * load genes information
deg <- readRDS(here(data_dir, subdir, de_outfnm))
fpdeg <- readRDS(here(data_dir, subdir, fpde_outfnm))
tndeg <- readRDS(here(data_dir, subdir, tnde_outfnm))
sc_data_list <- readRDS(here(data_dir, subdir, scdata_sumfnm))
sc_genes <- rownames(sc_data_list$cnt)

## * retrieve de/nonde-related genes
deg <- myt$stat_geneset(sc_genes, deg)
fpdeg <- myt$stat_geneset(sc_genes, fpdeg)
tndeg <- myt$stat_geneset(sc_genes, tndeg)



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


## * main
## ** vifit
mystanvifit <- load_stan_vi(here(
  exp_dir, exp_sub_dir,
  stan_dir, vi_dir, vifnm
))

## ** mcfit
mystanmcfit <- load_stan_mc(here(
  exp_dir, exp_sub_dir,
  stan_dir, mc_dir
), mcprefix)

## ** performance analyze
mystanfit <- mystanvifit

## dmucond: delta mucond
dmucond <- get_ctrlmnscase_par(mystanfit, par = "MuCond")
zerorela2dmucond <- mypntrela2rgn(dmucond, myprobs = c(0.025, 0.975), pnt = 0.0)
pred_deg <- sc_genes[zerorela2dmucond != 0]
pred_ndeg <- setdiff(sc_genes, pred_deg)

pred_tpg <- intersect(pred_deg, tpdeg)
pred_tng <- intersect(pred_ndeg, tndeg)
pred_fpg <- setdiff(pred_deg, tpdeg)
pred_fng <- setdiff(pred_ndeg, tndep)

tp <- length(pred_tpg)
tn <- length(pred_tng)
fp <- length(pred_fpg)
fn <- length(pred_fng)

myevals <- c(
  myroc$tpr(tp, fp, tn, fn),
  myroc$fpr(tp, fp, tn, fn),
  myroc$fdr(tp, fp, tn, fn),
  myroc$f1(tp, fp, tn, fn)
)

names(myevals) <- c("TPR", "FPR", "FDR", "F1")
message(str(myevals))
