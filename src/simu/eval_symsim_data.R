suppressPackageStartupMessages(library(SymSim))
suppressPackageStartupMessages(library(tidyverse))
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
import::from(here, here)
import::from(stringr, str_glue)
options("import.path" = here("rutils"))
myt <- modules::import("transform")
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
getDEgenes <- function(true_counts_res, popA_idx, popB_idx) {
    meta_cell <- true_counts_res$cell_meta
    meta_gene <- true_counts_res$gene_effects
    ## popA_idx <- which(meta_cell$pop == popA)
    ## popB_idx <- which(meta_cell$pop == popB)
    ngenes <- dim(true_counts_res$gene_effects[[1]])[1]

    DEstr <- sapply(strsplit(colnames(meta_cell)[which(grepl("evf", colnames(meta_cell)))], "_"), "[[", 2)
    param_str <- sapply(strsplit(colnames(meta_cell)[which(grepl("evf", colnames(meta_cell)))], "_"), "[[", 1)
    n_useDEevf <- sapply(1:ngenes, function(igene) {
        return(sum(abs(meta_gene[[1]][igene, DEstr[which(param_str == "kon")] == "DE"]) - 0.001 > 0) +
            sum(abs(meta_gene[[2]][igene, DEstr[which(param_str == "koff")] == "DE"]) - 0.001 > 0) +
            sum(abs(meta_gene[[3]][igene, DEstr[which(param_str == "s")] == "DE"]) - 0.001 > 0))
    })

    kon_mat <- true_counts_res$kinetic_params[[1]]
    koff_mat <- true_counts_res$kinetic_params[[2]]
    s_mat <- true_counts_res$kinetic_params[[3]]

    logFC_theoretical <- sapply(1:ngenes, function(igene) {
        return(log2(mean(s_mat[igene, popA_idx] * kon_mat[igene, popA_idx] / (kon_mat[igene, popA_idx] + koff_mat[igene, popA_idx])) /
            mean(s_mat[igene, popB_idx] * kon_mat[igene, popB_idx] / (kon_mat[igene, popB_idx] + koff_mat[igene, popB_idx]))))
    })

    true_counts <- true_counts_res$counts
    true_counts_norm <- t(t(true_counts) / colSums(true_counts)) * 10^6

    wil.p_true_counts <- sapply(1:ngenes, function(igene) {
        return(wilcox.test(true_counts_norm[igene, popA_idx], true_counts_norm[igene, popB_idx])$p.value)
    })

    wil.adjp_true_counts <- p.adjust(wil.p_true_counts, method = "fdr")

    return(list(nDiffEVF = n_useDEevf, logFC_theoretical = logFC_theoretical, wil.p_true_counts = wil.adjp_true_counts))
}

test_de_anygroup <- function() {
    tmpsymsim_true <- load_symsim(2)
    degenes <- getDEgenes(tmpsymsim_true, 1:100, 200:300)
}
## ** TODO gene module and PCA

## * batch effect
DivideBatches <- function(observed_counts_res, nbatch,
                          batch_effect_size = 1) {
    ## add batch effects to observed counts
    # use different mean and same sd to generate
    # the multiplicative factor for different gene in different batch
    observed_counts <- observed_counts_res[["counts"]]
    meta_cell <- observed_counts_res[["cell_meta"]]
    ncells <- dim(observed_counts)[2]
    ngenes <- dim(observed_counts)[1]
    batchIDs <- sample(1:nbatch, ncells, replace = TRUE)
    meta_cell2 <- data.frame(batch = batchIDs, stringsAsFactors = F)
    meta_cell <- cbind(meta_cell, meta_cell2)

    mean_matrix <- matrix(0, ngenes, nbatch)
    gene_mean <- rnorm(ngenes, 0, 1)
    temp <- lapply(1:ngenes, function(igene) {
        return(runif(nbatch, min = gene_mean[igene] - batch_effect_size, max = gene_mean[igene] + batch_effect_size))
    })
    mean_matrix <- do.call(rbind, temp)

    batch_factor <- matrix(0, ngenes, ncells)
    for (igene in 1:ngenes) {
        for (icell in 1:ncells) {
            batch_factor[igene, icell] <- rnorm(n = 1, mean = mean_matrix[igene, batchIDs[icell]], sd = 0.01)
        }
    }
    observed_counts <- round(2^(log2(observed_counts) + batch_factor))
    return(list(counts = observed_counts, cell_meta = meta_cell))
}

## * case / control design
Phyla3 <- function(plotting = F) {
    # par(mfrow=c(2,2))
    phyla <- rtree(2)
    phyla <- compute.brlen(phyla, 1)
    tip <- rtree(2)
    tip <- compute.brlen(phyla, 1)
    phyla <- bind.tree(phyla, tip, 1)
    phyla <- compute.brlen(phyla, c(1, 1, 1, 2))
    edges <- cbind(phyla$edge, phyla$edge.length)
    edges <- cbind(c(1:length(edges[, 1])), edges)
    connections <- table(c(edges[, 2], edges[, 3]))
    root <- as.numeric(names(connections)[connections == 2])
    tips <- as.numeric(names(connections)[connections == 1])
    phyla$tip.label <- as.character(tips)

    if (plotting == T) {
        plot(phyla, show.tip.label = F, lwd = 2)
        tiplabels(cex = 2)
        nodelabels(cex = 2)
    }
    return(phyla)
}
