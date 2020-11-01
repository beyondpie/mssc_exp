options(error = traceback)
options(warn = -1)
library(optparse)
suppressPackageStartupMessages(library(SymSim))
suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)
import::from(stringr, str_glue)
suppressPackageStartupMessages(library(ggpubr))

options("import.path" = here("rutils"))
myt <- modules::import("transform")
mysymsim <- modules::import("mysymsim")

## * configs
todatadir <- here("data", "symsim", "nobatch")
ncell <- 2000
ngene <- 300
nbatch <- 10
nevf <- 10
n_de_evf <- 6
sigma <- 0.2
vary <- "all"
myseed <- 0L

## * generate data
symsimtrue <- mysymsim$sim_symsim_true(myseed = myseed,
  ncell = ncell,
  ngene = ngene,
  hasgenemodule = F,
  npop = 2,
  nevf = nevf,
  n_de_evf = n_de_evf,
  sigma = sigma, vary = vary)

symsimumi <- mysymsim$sim_symsim_obs("UMI", symsimtrue)
batchids <- mysymsim$assign_batches_2pop(symsimumi, nbatch)
symsimumi$batch_meta <- list(batch = batchids)
cnt <- symsimumi$counts
totcntpcell <- colSums(cnt)

## * output the train data
saveRDS(object = symsimumi,
  file = paste(todatadir, str_glue("{myseed}.rds"), sep = "/"))
myt$to_bagwiff(cnt, batchids, symsimumi$cell_meta$pop,
  totcntpcell,
  outf = paste(todatadir, str_glue("{myseed}.rdump"), sep = "/"),
  rdump = T)
