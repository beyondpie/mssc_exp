options(error = traceback)
options(warn = -1)
library(optparse)
suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)
import::from(stringr, str_glue)
options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * options from outside
option_list <- list(
  make_option(c("--myseed"),
    action = "store",
    type = "integer",
    default = 1
  )
)

args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()

## * configs
symsim_data_dir <- here("data", "symsim",
  "twostage_be_symsim", "data")
fprefix <- "symsim_2be"
myseed <- args$myseed

## * loading data
symsimumi2be <- readRDS(paste(symsim_data_dir,
  str_glue("{fprefix}_{myseed}.rds"),
  sep = "/"
))

## * transform to rstan
cnt <- symsimumi2be$counts
batches <- symsimumi2be$batch_meta$batch
conds <- symsimumi2be$cell_meta$pop
totcntpcell <- colSums(cnt)

myt$to_bagwiff(
  cnt, batches, conds, totcntpcell,
  outf = paste(symsim_data_dir,
               str_glue("{fprefix}_{myseed}.rdump"),
               sep = "/"),
  rdump = T
)
