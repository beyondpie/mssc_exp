options(error = traceback)
options(warn = -1)

library(Seurat)
library(rstan)
library(ggplot2)
library(bayesplot)
library(data.table)

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(import::from(here,here))
import::from(optparse, make_option, OptionParser, parse_args)
import::from(stringr, str_glue)

options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * load genes information
data_dir <- "data"
subdir <- "UM"

de_outfnm <- "tcga_diffexp_genes.rds"
fpde_outfnm <- "tcga_fp_diffexp_genes.rds"
tnde_outfnm <- "tcga_tn_diffexp_genes.rds"

deg <- readRDS(here(data_dir, subdir, de_outfnm))
fpdeg <- readRDS(here(data_dir, subdir, fpde_outfnm))
tndeg <- readRDS(here(data_dir, subdir, tnde_outfnm))

scRNAseq_sumfnm <- "sampled_scRNAseq_summary.rds"
sc_data_list <- readRDS(here(data_dir, subdir, scRNAseq_sumfnm))
sc_genes <- rownames(sc_data_list$cnt)

## ** retrieve de/nonde-related genes
deg <- myt$stat_geneset(sc_genes, deg)
fpdeg <- 


## * utils for gettign rstan results
load_stan_vi <- function(path) {
  rstan::read_stan_csv(path)
}

load_stan_mc <- function(dirpath, modelnm) {
  csvfiles <- dir(path = dirpath, pattern = paste0(modelnm, "[0-9].csv"),
                  full.names = T)
  rstan::read_stan_csv(csvfiles)
}

## TODO: add explanation about which is case and control
get_ctrlmnscase_par <- function(mystanfit, par="MuCond") {
  mus <- rstan::extract(mystanfit, pars = par)[[par]]
  ## TODO: add colnames of delta: gsymbols
  delta <- as.data.frame(mus[, , 1] - mus[, , 2])
  return(delta)
}

## point relative to regions
mypntrela2rgn <- function(myintvals, myprobs = c(0.025, 0.975), pnt=0.0) {
  q <- as.data.frame(apply(myintvals, 2, quantile, probs = myprobs))
  apply(q, 2, myt$pntrela2rgn, pnt = pnt)
}






