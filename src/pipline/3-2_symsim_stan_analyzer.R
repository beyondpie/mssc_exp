options(error = traceback)
options(warn = -1)
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(bayesplot))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(import::from(here, here))
import::from(optparse, make_option, OptionParser, parse_args)
import::from(stringr, str_glue)

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myroc <- modules::import("roc")
mysymsim <- modules::import("mysymsim")

## * options
option_list <- list(
  make_option(c("--myseed"),
    action = "store",
    type = "integer",
    default = 1
  ),
  make_option(c("--version"),
    action = "store",
    type = "character",
    default = "v1-1"),
  make_option(c("--method"),
    action = "store",
    type = "character",
    default = "vi")
)

args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()

## * configs
myseed <- args$myseed
mvrsn <- args$version
method <- args$method
symsim_data_dir <- here("data", "symsim", "twostage_be_symsim", "data")
symsim_exp_dir <- here("exps", "symsim", "stan", myseed)
symsim_exp_vi_dir <- paste(symsim_exp_dir, "vi")
symsim_exp_mc_dir <- paste(symsim_exp_dir, "mc")

## * get truth
symsimtrue <- readRDS(file =
  paste(symsim_data_dir,
    str_glue("symsim_true_{myseed}.rds"),
    sep = "/"))
symsim_dea <- mysymsim$symsim_de_analysis(symsimtrue,
  popA_idx = which(symsimtrue$cell_meta$pop == 1),
  popB_idx = which(symsimtrue$cell_meta$pop == 2))

numofgene <- nrow(symsimtrue$counts)

symsim_degenes <- mysymsim$get_symsim_degenes(symsim_dea,
  nDiffEVF = 1,
  logFC = 0.6) %>% which(. == T)

symsim_ndegs <- setdiff(seq_len(numofgene), symsim_degenes)

symsim_strict_ndegenes <- mysymsim$get_symsim_strict_ndegenes(symsim_dea,
  nDiffEVF = 0,
  logFC = 0.5) %>% which(. == T)

symsim_zerodiffevf_genes <- mysymsim$get_symsim_ndiffevf_genes(symsim_dea) %>%
  which(. == T)

symsim_sampled_ndegs_1 <- sample(symsim_ndegs, size = length(symsim_degenes),
  replace = F)

symsim_sampled_ndegs_2 <- sample(symsim_ndegs, size = length(symsim_degenes),
  replace = F)

symsim_sampled_ndegs_3 <- sample(symsim_ndegs, size = length(symsim_degenes),
  replace = F)

## * load simulated counts
symsimumi <- readRDS(
  file = paste(symsim_data_dir,
    str_glue("symsim_umi_{myseed}.rds"),
    sep = "/"))
symsimbe <- readRDS(
  file = paste(symsim_data_dir,
    str_glue("symsim_be_{myseed}.rds"),
    sep = "/"))
symsim2be <- readRDS(
  file = paste(symsim_data_dir,
    str_glue("symsim_2be_{myseed}.rds"),
    sep = "/"))

## * get results
getauc <- function(resdir, modelvsn = "v1-1", method = "vi", par = "MuCond",
                   degs, ndegs, seed = myseed) {
  stanfit <- myt$load_stan(dirnm = resdir, modelnm = modelvsn, method = method,
    vi_dir = "vi", mc_dir = "mc")
  dcond <- myt$get_ctrlmnscase_par(stanfit, par = par)
  dt <- abs(myt$calt(dcond, fn = colMeans))
  dt <- dt[c(degs, ndegs)]
  labels <- c(rep(1L, length(degs)), rep(0L, length(ndegs)))
  auc <- myt$fmtflt(caTools::colAUC(dt, labels))
  result <- list(modelvsn = modelvsn,
    method = method,
    par = par,
    numofdeg = length(degs),
    numofndeg = length(ndegs),
    auc = auc,
    symsim_seed = seed
  )
  return(result)
}

## TODO: merge multiple pairs of degs and ndegs
symsimvi_auc <- function() {
  par <- "MuCond"
  symsim_vi_auc_strict <- getauc(resdir = symsim_exp_dir,
    modelvsn = mvrsn, method = "vi", par = "MuCond",
    degs = symsim_degenes, ndegs = symsim_strict_ndegenes,
    seed = myseed)
  message(
    str_glue("symsim {mvrsn} {method} on {par}: strict ndegs"))
  str(symsim_vi_auc_strict)

  symsim_vi_auc_all <- getauc(resdir = symsim_exp_dir,
    modelvsn = mvrsn, method = "vi", par = "MuCond",
    degs = symsim_degenes, ndegs = symsim_ndegs,
    seed = myseed)
  message(
    str_glue("symsim {mvrsn} {method} on {par}: all ndegs"))
  str(symsim_vi_auc_all)

  symsim_vi_auc_zerodiffevf <- getauc(resdir = symsim_exp_dir,
    modelvsn = mvrsn, method = "vi", par = "MuCond",
    degs = symsim_degenes, ndegs = symsim_zerodiffevf_genes,
    seed = myseed)
  message(
    str_glue("symsim {mvrsn} {method} on {par}: zero diffevf ndegs"))
  str(symsim_vi_auc_zerodiffevf)

    symsim_vi_auc_sample <- getauc(resdir = symsim_exp_dir,
      modelvsn = mvrsn, method = "vi", par = "MuCond",
      degs = symsim_degenes,
      ndegs = symsim_sampled_ndegs_1,
      seed = myseed)
    message(
      str_glue("symsim {mvrsn} {method} on {par}: sampled ndegs"))
    str(symsim_vi_auc_sample)


    symsim_vi_auc_sample <- getauc(resdir = symsim_exp_dir,
      modelvsn = mvrsn, method = "vi", par = "MuCond",
      degs = symsim_degenes,
      ndegs = symsim_sampled_ndegs_2,
      seed = myseed)
    message(
      str_glue("symsim {mvrsn} {method} on {par}: sampled ndegs"))
    str(symsim_vi_auc_sample)

    symsim_vi_auc_sample <- getauc(resdir = symsim_exp_dir,
      modelvsn = mvrsn, method = "vi", par = "MuCond",
      degs = symsim_degenes,
      ndegs = symsim_sampled_ndegs_3,
      seed = myseed)
    message(
      str_glue("symsim {mvrsn} {method} on {par}: sampled ndegs"))
    str(symsim_vi_auc_sample)

}

symsimvi_auc()
