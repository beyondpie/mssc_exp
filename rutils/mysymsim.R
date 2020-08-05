phyla2 <- function() {
    tmp <- ape::rtree(2, tip.label = c(1, 2))
    ape::compute.brlen(tmp, 1)
}

plotphylo <- function(mytree) {
    ape::plot.phylo(mytree, show.tip.label = F, lwd = 2)
    ape::tiplabels(cex = 2)
    ape::nodelabels(cex = 2)
}

## consider case and control in setting for cells
## case and control as two different cell populations.
assign_batch_for_cells <- function(symsimobs, nbatch) {
    left_max_batch <- floor(nbatch / 2)
    ncell <- ncol(symsimobs$counts)
    cellpops <- symsimobs$cell_meta$cell_meta.pop
    control_batches <- seq_len(left_max_batch)
    control_cells <- which(cellpops == 1)
    case_batches <- seq(left_max_batch + 1, nbatch)
    case_cells <- which(cellpops == 2)
    batchids <- rep(1, ncell)
    batchids[control_cells] <- sample(control_batches,
        length(control_cells),
        replace = T
    )
    batchids[case_cells] <- sample(case_batches,
        length(case_cells),
        replace = T
    )
    return(batchids)
}

## by setting oncells and batchids, we can add different batch effects
## to different group of cells.
add_batch_effect <- function(symsimobs, nbatch,
                             ongenes = seq_len(nrow(symsimobs$counts)),
                             oncells = seq_len(ncol(symsimobs$counts)),
                             batchids = assign_batch_for_cells(
                                 symsimobs, nbatch
                             ),
                             batch_effect_size = 1,
                             sd = 0.01) {
    ## add batch effects to observed counts
    # use different mean and same sd to generate the
    # multiplicative factor for different gene in different batch
    observed_counts <- symsimobs[["counts"]]
    ncells <- ncol(symsimobs$counts)
    ngenes <- nrow(symsimobs$counts)

    ## set gene-batch effect matrix
    ## only the ongenes and oncells will non-zeros have batch effects
    gene_mean <- rnorm(ngenes, 0, 1)
    mean_matrix <- matrix(0, ngenes, nbatch)
    for (i in ongenes) {
        for (j in oncells) {
            mean_matrix[i, j] <- runif(1,
                min = gene_mean[i] - batch_effect_size,
                max = gene_mean[i] + batch_effect_size
            )
        }
    }

    ## add batch effects the entire observed count matrix
    batch_factor <- matrix(0, ngenes, ncells)
    for (igene in seq_len(ngenes)) {
        for (icell in seq_len(ncells)) {
            batch_factor[igene, icell] <- rnorm(
                n = 1,
                mean = mean_matrix[igene, batchids[icell]],
                sd = sd
            )
        }
    }

    observed_counts <- round(2^(log2(observed_counts) + batch_factor))
    return(list(
        counts = observed_counts,
        gene_mean = gene_mean,
        mean_matrix = mean_matrix,
        batch_factor = batch_factor,
        batchids = batchids
    ))
}

## set getDEgenes on any two groups of cells
getDEgenes <- function(true_counts_res, popA_idx, popB_idx) {
    meta_cell <- true_counts_res$cell_meta
    meta_gene <- true_counts_res$gene_effects
    ## popA_idx <- which(meta_cell$pop == popA)
    ## popB_idx <- which(meta_cell$pop == popB)
    ngenes <- dim(true_counts_res$gene_effects[[1]])[1]

    DEstr <- sapply(
        strsplit(colnames(meta_cell)[which(
            grepl("evf", colnames(meta_cell))
        )], "_"), "[[", 2
    )
    param_str <- sapply(
        strsplit(colnames(meta_cell)[which(
            grepl("evf", colnames(meta_cell))
        )], "_"), "[[", 1
    )
    n_useDEevf <- sapply(1:ngenes, function(igene) {
        return(sum(abs(meta_gene[[1]][
            igene,
            DEstr[which(param_str == "kon")] == "DE"
        ]) - 0.001 > 0) +
            sum(abs(meta_gene[[2]][
                igene,
                DEstr[which(param_str == "koff")] == "DE"
            ]) - 0.001 > 0) +
            sum(abs(meta_gene[[3]][
                igene,
                DEstr[which(param_str == "s")] == "DE"
            ]) - 0.001 > 0))
    })
    kon_mat <- true_counts_res$kinetic_params[[1]]
    koff_mat <- true_counts_res$kinetic_params[[2]]
    s_mat <- true_counts_res$kinetic_params[[3]]
    logFC_theoretical <- sapply(1:ngenes, function(igene) {
        return(log2(mean(s_mat[igene, popA_idx] * kon_mat[igene, popA_idx] /
            (kon_mat[igene, popA_idx] + koff_mat[igene, popA_idx])) /
            mean(s_mat[igene, popB_idx] * kon_mat[igene, popB_idx]
                / (kon_mat[igene, popB_idx] + koff_mat[igene, popB_idx]))))
    })
    true_counts <- true_counts_res$counts
    true_counts_norm <- t(t(true_counts) / colSums(true_counts)) * 10^6
    wil.p_true_counts <- sapply(
        1:ngenes,
        function(igene) {
            return(wilcox.test(
                true_counts_norm[igene, popA_idx],
                true_counts_norm[igene, popB_idx]
            )$p.value)
        }
    )
    wil.adjp_true_counts <- p.adjust(wil.p_true_counts, method = "fdr")
    return(list(
        nDiffEVF = n_useDEevf,
        logFC_theoretical = logFC_theoretical,
        wil.p_true_counts = wil.adjp_true_counts
    ))
}
