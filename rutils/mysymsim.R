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
