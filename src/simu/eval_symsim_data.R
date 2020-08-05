suppressPackageStartupMessages(library(SymSim))
suppressPackageStartupMessages(library(tidyverse))
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
import::from(here, here)
import::from(stringr, str_glue)
options("import.path" = here("rutils"))
myt <- modules::import("transform")
mysymsim <- modules::import("mysymsim")
## use modules::reload(myt) to update myt

## * configs
myseed <- 0L

## * funcitons
symsim2ptsne <- function(symsimdata, issimtrue = T, showbatch = F,
                         seed = myseed, scale_factor = myt$umi_scale_factor,
                         title = "SymSim",
                         fnm = NULL) {
    if (issimtrue) {
        celltypes <- symsimdata$cell_meta$pop
    } else {
        celltypes <- symsimdata$cell_meta$cell_meta.pop
    }
    if (showbatch) {
        celltypes <- symsimdata$cell_meta$batch
    }

    cnt <- symsimdata$counts
    logcpm <- myt$getlogcpm(cnt, scale_factor)
    ttsne <- myt$gettsne(logcpm, seed)
    myt$plottsne(ttsne$Y, celltypes, title, fnm)
}

load_symsim <- function(seed = 1, issimtrue = T,
                        ncell = 5000, ngene = 270,
                        npop = 5, nminmodule = 50) {
    prefix <- ifelse(issimtrue, "true", "umi")
    suffix <- str_glue("{ncell}_{ngene}_{npop}_{nminmodule}_{seed}.rds")
    ## hard code dirnm
    fnm <- here("src", "simu", "symsimdata", str_glue("{prefix}_{suffix}"))
    readRDS(fnm)
}

eval_symsim_true <- function() {
    ## ** SymSim_true
    p_symsimtrue_list <- lapply(1:10, FUN = function(i) {
        symsim_true <- load_symsim(i)
        symsim2ptsne(symsim_true, seed = myseed)
    })

    p_allsymsim <- ggarrange(plotlist = p_symsimtrue_list, nrow = 3, ncol = 4)
    ggsave(here("src", "simu", "symsim_data_quality_pdf", "symsim_true_tsnes.pdf"),
        plot = p_allsymsim, width = 20, height = 10
    )
}
eval_symsim_umi <- function() {
    ## ** SymSim_umi
    ## view tsne of umi for cell types
    p_symsimumi_list <- lapply(1:10, FUN = function(i) {
        symsim_true <- load_symsim(i, issimtrue = F)
        symsim2ptsne(symsim_true, seed = myseed, issimtrue = F)
    })

    p_allsymsim_umis <- ggarrange(plotlist = p_symsimumi_list, nrow = 3, ncol = 4)
    ggsave(here("src", "simu", "symsim_data_quality_pdf", "symsim_umi_tsnes.pdf"),
        plot = p_allsymsim_umis, width = 20, height = 10
    )

    ## view tsne of umi for batch effects
    p_symsimumi_batch_list <- lapply(1:10, FUN = function(i) {
        symsim_true <- load_symsim(i, issimtrue = F)
        symsim2ptsne(symsim_true, seed = myseed, issimtrue = F, showbatch = T)
    })

    p_allsymsim_umis_batch <- ggarrange(plotlist = p_symsimumi_batch_list, nrow = 3, ncol = 4)
    ggsave(here("src", "simu", "symsim_data_quality_pdf", "symsim_umi_batch_tsnes.pdf"),
        plot = p_allsymsim_umis_batch, width = 20, height = 10
    )

    ## *** cell library size
    ## we use the second simualtion since it has good isolation of cell population.
    symsim_true <- load_symsim(2)
    phist_true <- myt$plothist(symsim_true$counts,
        mytitle = "SymSim_True Cell Library Size",
        mybinwidth = 1000, max_xlim = 2e+05
    )

    symsim_umi <- load_symsim(2, issimtrue = F)
    phist_umi <- myt$plothist(symsim_umi$counts,
        mytitle = "SymSim_UMI Cell Library Size",
        max_xlim = 2e+04,
        mybinwidth = 300
    )

    phists <- ggarrange(plotlist = list(phist_true, phist_umi), nrow = 1, ncol = 2)
    ggsave(here("src", "simu", "symsim_data_quality_pdf", "symsim_hist_compare.pdf"),
        plot = phists, width = 15, height = 10
    )

    ## *** zeros in cells
    true_zeros <- symsim_true$counts
    true_zeros[symsim_true$counts == 0] <- 1
    true_zeros[symsim_true$counts > 0] <- 0
    true_zeros <- true_zeros / nrow(true_zeros)
    phist_zeros_true <- myt$plothist(true_zeros,
        mytitle = "SymSim_TRUE Zeros in Cells",
        mybinwidth = 0.05, max_xlim = 1
    )

    umi_zeros <- symsim_umi$counts
    umi_zeros[symsim_umi$counts == 0] <- 1
    umi_zeros[symsim_umi$counts > 0] <- 0
    umi_zeros <- umi_zeros / nrow(umi_zeros)
    phist_zeros_umi <- myt$plothist(umi_zeros,
        mytitle = "SymSim_UMI Zeros in Cells",
        mybinwidth = 0.05, max_xlim = 1
    )
    phist_zeros <- ggarrange(
        plotlist = list(phist_zeros_true, phist_zeros_umi),
        nrow = 1, ncol = 2
    )
    ggsave(here("src", "simu", "symsim_data_quality_pdf", "symsim_zeros_compare.pdf"),
        plot = phist_zeros, width = 15, height = 10
    )
}
## * gene related
## ** test DE genes on any groups of cells
test_de_anygroup <- function() {
    tmpsymsim_true <- load_symsim(2)
    degenes <- mysymsim$getDEgenes(tmpsymsim_true, 1:100, 200:300)
}
## ** TODO gene module and PCA

## * simulation design
## ** case and control experiments (true positive)
## ** add batch effect (false positive)
## * evaluate simulation

























