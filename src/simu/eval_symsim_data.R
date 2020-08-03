library(SymSim)
suppressPackageStartupMessages(library(tidyverse))
library(ggplot2)
library(RColorBrewer)
import::from(here, here)
options("import.path" = here("rutils"))
myt <- modules::import("transform")
## use modules::reload(myt) to update myt

## * Configs
myseed <- 0L

## * load saved data
symsim_true <- readRDS("symsim_true.rds")
symsim_umi <- readRDS("symsim_umi.rds")
## put full-length support for later
## symsim_fullen <- readRDS("symsim_fullen.rds")

## * cell related
## ** SymSim_true
logcpm_true <- myt$getlogcpm(myt$getcpm(symsim_true$counts,
                                        myt$umi_scale_factor))
tsne_true <- myt$gettsne(logcpm_true, myseed)
## ptsne_true <- myt$plottsne(tsne_true$Y, symsim_true$cell_meta$pop,
    ## title = "SymSim_True TSNE", fnm = "symsim_true_tsne.pdf"
## )
ptsne_true <- myt$plottsne(tsne_true$Y, symsim_true$cell_meta$pop,
    title = "SymSim_True TSNE"
)

## *** cell library size
phist_true <- myt$plothist(symsim_true$counts, mytitle = "SymSim_True Cell Library Size",
                           mybinwidth = 1000, max_xlim = 2e+05,
                           fnm = "symsim_true_cnt_hist.pdf")

## *** zeros in cells
true_zeros <- symsim_true$counts
true_zeros[symsim_true$counts == 0] <- 1
true_zeros[symsim_true$counts > 0] <- 0
true_zeros <- true_zeros / nrow(true_zeros)
phist_zeros_true <- myt$plothist(true_zeros,
    mytitle = "SymSim_TRUE Zeros in Cells",
    mybinwidth = 0.05, max_xlim = 1,
    fnm = "symsim_true_zeros_hist.pdf"
)

## ** SymSim_true
logcpm_true <- myt$getlogcpm(myt$getcpm(
    symsim_umi$counts,
    myt$umi_scale_factor
))
tsne_umi <- myt$gettsne(logcpm_umi, myseed)
ptsne_umi <- myt$plottsne(tsne_umi$Y, symsim_umi$cell_meta$cell_meta.pop,
                          title = "SymSim_UMI TSNE",
                          fnm = "symsim_umi_tsne.pdf"
)
ptsne_batch_umi <- myt$plottsne(tsne_umi$Y, symsim_umi$cell_meta$batch,
                                title="SymSim_UMI Batch effect TSNE",
                                fnm="symsim_umi_batch_effect.pdf")
## *** cell library size
phist_umi <- myt$plothist(symsim_umi$counts,
    mytitle = "SymSim_UMI Cell Library Size",
    max_xlim = 1e+04,
    mybinwidth = 300,
    fnm = "symsim_umi_cnt_hist.pdf"
)

## *** zeros in cells
umi_zeros <- symsim_umi$counts
umi_zeros[symsim_umi$counts == 0]  <- 1
umi_zeros[symsim_umi$counts > 0 ]  <- 0
umi_zeros <- umi_zeros / nrow(umi_zeros)
phist_zeros_umi <- myt$plothist(umi_zeros, mytitle = "SymSim_UMI Zeros in Cells",
                                mybinwidth = 0.05, max_xlim = 1,
                                fnm = "symsim_umi_zeros_hist.pdf")
## * gene related


## ** gene module and PCA
