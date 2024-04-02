options(error = traceback)
options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(import::from(here, here))
suppressPackageStartupMessages(library(DESeq2))
import::from(optparse, make_option, OptionParser, parse_args)
import::from(stringr, str_glue)

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myroc <- modules::import("roc")
mysymsim <- modules::import("mysymsim")
mypseudo <- modules::import("pseudobulk")

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
    default = "v1-1")
)

args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()

## * configs
myseed <- args$myseed
mvrsn <- args$version
symsim_data_dir <- here("data", "symsim", "twostage_be_symsim", "data")

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

## * pseudo bulk analysis
mycnt <- symsim2be$counts
mybatches <- symsim2be$batch_meta$batch
myconds <- symsim2be$cell_meta$pop
pseudo_deseq2_res <- mypseudo$pseudobulk_deseq2(mycnt, mybatches, myconds)

tmp <- mypseudo$calc_auc(pseudo_deseq2_res, symsim_degenes,
                         symsim_strict_ndegenes)

tmp <- mypseudo$calc_auc(pseudo_deseq2_res, symsim_degenes,
                         symsim_ndegs)

tmp <- mypseudo$calc_auc(pseudo_deseq2_res, symsim_degenes,
                        symsim_zerodiffevf_genes)

tmp <- mypseudo$calc_auc(pseudo_deseq2_res, symsim_degenes,
  symsim_sampled_ndegs_1)

tmp <- mypseudo$calc_auc(pseudo_deseq2_res, symsim_degenes,
  symsim_sampled_ndegs_2)

tmp <- mypseudo$calc_auc(pseudo_deseq2_res, symsim_degenes,
  symsim_sampled_ndegs_3)
