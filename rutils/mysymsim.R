
myggtitle <- theme(plot.title = element_text(size = 15, hjust = 0.5))

phyla2 <- function() {
  tmp <- ape::rtree(2, tip.label = c(1, 2))
  ape::compute.brlen(tmp, 1)
}

plotphylo <- function(mytree) {
  ape::plot.phylo(mytree, show.tip.label = F, lwd = 2)
  ape::tiplabels(cex = 2)
  ape::nodelabels(cex = 2)
}


sim_symsim_obs <- function(protocol, symsimtrue) {
  data(gene_len_pool, package = "SymSim")
  ngene <- nrow(symsimtrue$counts)
  gene_len <- sample(gene_len_pool, ngene, replace = FALSE)

  ## UMI settings
  depth_mean_umi <- 45000
  depth_sd_umi <- 4500
  alpha_mean_umi <- 0.1
  ## fullength settings
  depth_mean_fullength <- 1e+05
  depth_sd_fullength <- 10000
  alpha_mean_fullength <- 0.4

  depth_mean <- ifelse(protocol == "UMI",
    depth_mean_umi, depth_mean_fullength
  )
  depth_sd <- ifelse(protocol == "UMI", depth_sd_umi, depth_sd_fullength)
  alpha_mean <- ifelse(protocol == "UMI",
    alpha_mean_umi, alpha_mean_fullength
  )
  SymSim::True2ObservedCounts(
    true_counts = symsimtrue$counts,
    meta_cell = symsimtrue$cell_meta,
    protocol = protocol, alpha_mean = alpha_mean,
    alpha_sd = 0.02, gene_len = gene_len, depth_mean = depth_mean,
    depth_sd = depth_sd, nPCR1 = 14
  )
}


## consider case and control in setting for cells
## case and control as two different cell populations.
assign_batch_for_cells <- function(symsimobs, nbatch) {
  left_max_batch <- floor(nbatch / 2)
  ncell <- ncol(symsimobs$counts)
  cellpops <- symsimobs$cell_meta$pop
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
                             onbatches = c(1, 2, 5),
                             batchids = assign_batch_for_cells(
                               symsimobs, nbatch
                             ),
                             batch_effect_size = rep(1.0, nbatch),
                             sd = 0.01,
                             gene_mean_std = 0.2,
                             isg2brandom = F) {
  ## add batch effects to observed counts
  # use different mean and same sd to generate the
  # multiplicative factor for different gene in different batch
  observed_counts <- symsimobs[["counts"]]
  ncells <- ncol(symsimobs$counts)
  ngenes <- nrow(symsimobs$counts)

  ## set gene-batch effect matrix
  gene_mean <- rnorm(ngenes, 0, gene_mean_std)
  mean_matrix <- matrix(0, ngenes, nbatch)
  for (i in ongenes) {
    for (j in seq_len(nbatch)) {
      if (isg2brandom) {
        mean_matrix[i, j] <- runif(1,
          min = gene_mean[i] - batch_effect_size[j],
          max = gene_mean[i] + batch_effect_size[j]
        )
      } else {
        mean_matrix[i, j] <- gene_mean[i] + batch_effect_size[j]
      }
    }
  }

  ## add batch effects the entire observed count matrix
  batch_factor <- matrix(0, ngenes, ncells)
  for (i in seq_len(ngenes)) {
    for (j in seq_len(ncells)) {
      batch_factor[i, j] <- rnorm(
        n = 1,
        mean = mean_matrix[i, batchids[j]],
        sd = sd
      )
    }
  }

  observed_counts <- round(2^(log2(observed_counts) + batch_factor))

  ## double check observed_counts counts data as integer
  intc <- apply(observed_counts, c(1, 2), function(x) {
    ifelse(x > 0, as.integer(x + 1), 0L)
  })

  ## convert to symsim related structure
  result <- symsimobs
  result$counts <- intc
  result$batch_meta <- list(
    batch = batchids,
    g2b = mean_matrix,
    gmean = gene_mean,
    g2c = batch_factor,
    sd = sd
  )
  invisible(result)
}

## set getDEgenes on any two groups of cells
symsim_de_analysis <- function(true_counts_res, popA_idx, popB_idx) {
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

get_symsim_degenes <- function(symsim_dea, nDiffEVF = 1, logFC = 0.6) {
  invisible((symsim_dea$nDiffEVF >= nDiffEVF) &
    (symsim_dea$logFC_theoretical >= logFC))
}

get_symsim_strict_ndegnes <- function(symsim_dea, nDiffEVF = 0, logFC = 0.1) {
  invisible((symsim_dea$nDiffEV <= nDiffEVF) &
    (symsim_dea$logFC_theoretical <= logFC))
}

plotviolin <- function(symsimdata, genes) {
  library(tidyverse)
  ## library(ggpubr)
  n <- length(genes)
  gnm <- paste0("gene", genes)
  select_cnt <- symsimdata$counts[genes, ]
  cell2batch <- symsimdata$batch_meta$batch
  nbatch <- max(cell2batch)
  plotdata <- as.data.frame(t(select_cnt))
  colnames(plotdata) <- gnm

  ## use factor to order the plot
  ## in the order ind1, ind2, ..., not ind1. ind10, ind2, ...
  plotdata$pop <- factor(paste0("ind", cell2batch),
    levels = paste0("ind", seq_len(nbatch)),
    ordered = TRUE
  )

  plist <- lapply(seq_len(n), FUN = function(i) {
    tmp <- plotdata[c(gnm[i], "pop")]
    ggplot(tmp, aes_string(x = "pop", y = gnm[i], fill = "pop")) +
      geom_violin() +
      myggtitle +
      theme(
        legend.position = "none",
        axis.title.x = element_blank()
      )
  })
  return(plist)
}
