## simulate individual effect (bf) v2 using SymSim

## TODO:
## - model gene modules

## * set R environment
suppressWarnings(suppressMessages({
  library(SymSim)
  library(tidyverse)
  library(grid)
  library(gtable)
  library(gridExtra)
  library(ggpubr)
  library(ape)
  library(Rtsne)
  library(cmdstanr)
  library(bayesplot)
  library(posterior)
}))

## warnings/errors traceback settings
options(error = traceback)
options(warn = -1)
options(mc.cores = 3)

## load help functions
options("import.path" = list(here::here("rutils"),
  here::here("src", "mssc")))
mypseudo <- modules::import("pseudobulk")

## * modify symsim functions
## load from a local file
library(phytools)
library(MASS)

mysymsim <- modules::import("simulation_functions")
## modules::reload(mysymsim)

## * help functions
set_population_struct_using_phylo <- function(w = rep(0.1, 6)) {
  ## two populations (two conditiosn), while each has some subpopulations
  ## used when individual effects as subpopulations in phylogenetic tree
  ## w: vector of branch length in phylo tree, such as (0.1,0.2,0.3, 0.2,0.1,0.2)
  ##   - one minus the branch length describes the correlation
  ##     between the parent and the tip, so should be less than 1
  ##   - 1 to length(w)/2 describes the subpopulation structure in one case,
  ##     rest the other

  n <- ceiling(length(w) / 2)
  sub1 <- 1:n
  sub2 <- (n + 1):length(w)
  ## generate newick tree format
  newick_tree <- paste("((", paste(sub1, w[sub1], sep = ":", collapse = ","), "):1, (",
    paste(sub2, w[sub2], sep = ":", collapse = ","), "):1);")
  p <- ape::read.tree(text = newick_tree)
  p$pop <- list(sub1 = sub1, sub2 = sub2)
  p$w <- w
  return(invisible(p))
}

get_symsim_umi <- function(phyla, bimod = 0, capt_alpha = 0.1,
                           gene_module_prop = 0.0,
                           ncell_per_subpop = 100,
                           ngene = 100,
                           seed = 1L, nevf = 10, n_de_evf = 7,
                           sigma = 0.6, vary = "s") {
  ## - bimod: bimodility in expressions
  ##   - In SymSim, be default, half of the genes will use this feature
  ##   - Though we can use continuous value between [0, 1], in SymSim examples,
  ##     they only use 0 or 1.
  ## - capt_alpha: capture efficiency mean
  ##   - used to generate the observed counts, default is 0.1
  ##   - we could vary this from (0.0, 0.2] by 0.05 step size.
  ## - sigma: std of evf within the same population
  ##   - It determines gene expression variations whithn one population.
  ##   - default as 0.2 or 0.6

  npop <- length(phyla$tip.label)
  ncell_total <- npop * ncell_per_subpop
  ## simulate the true counts per cell
  symsim_true <- SymSim::SimulateTrueCounts(
    ncells_total = ncell_total,
    ngene = ngene,
    phyla = phyla,
    bimod = bimod,
    ## scale kinetic parameter:s
    scale_s = 1,
    param_realdata = "zeisel.imputed",
    geffect_mean = 0,
    gene_effects_sd = 1,
    ## prob gene_effect not zero
    gene_effect_prob = 0.3,
    evf_center = 1,
    Sigma = sigma,
    evf_type = "discrete",
    nevf = nevf,
    n_de_evf = n_de_evf,
    vary = vary,
    min_popsize = floor(ncell_total / npop),
    prop_hge = 0.0,
    gene_module_prop = gene_module_prop,
    randseed = seed
  )
  ## simulate reads per cell under UMI-based scRNA sequencing
  data(gene_len_pool, package = "SymSim")
  gene_len <- sample(gene_len_pool, ngene, replace = FALSE)
  symsim_umi <- SymSim::True2ObservedCounts(
    true_counts = symsim_true$counts,
    meta_cell = symsim_true$cell_meta,
    protocol = "UMI",
    alpha_mean = capt_alpha,
    alpha_sd = 0.002,
    gene_len = gene_len,
    depth_mean = 45000,
    depth_sd = 4500,
    lenslope = 0.02,
    nbins = 20,
    amp_bias_limit = c(-0.2, 0.2),
    nPCR1 = 14,
    nPCR2 = 10,
    LinearAmp = F,
    LinearAmp_coef = 2000)
  return(invisible(list(true = symsim_true, umi = symsim_umi)))
}

getDEgenes <- function(true_counts_res, popA, popB) {
  meta_cell <- true_counts_res$cell_meta
  meta_gene <- true_counts_res$gene_effects
  ## when popA or popB has multiple subpopulations
  popA_idx <- which(meta_cell$pop %in% popA)
  popB_idx <- which(meta_cell$pop %in% popB)
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

  logFC_theoretical <- sapply(1:ngenes, function(igene)
    return(log2(mean(s_mat[igene, popA_idx] * kon_mat[igene, popA_idx] / (kon_mat[igene, popA_idx] + koff_mat[igene, popA_idx])) /
      mean(s_mat[igene, popB_idx] * kon_mat[igene, popB_idx] / (kon_mat[igene, popB_idx] + koff_mat[igene, popB_idx])))))

  true_counts <- true_counts_res$counts
  true_counts_norm <- t(t(true_counts) / colSums(true_counts)) * 10^6

  ## wil.p_true_counts <- sapply(1:ngenes, function(igene)
  ##   return(wilcox.test(true_counts_norm[igene, popA_idx], true_counts_norm[igene, popB_idx])$p.value))

  ## wil.adjp_true_counts <- p.adjust(wil.p_true_counts, method = 'fdr')

  return(list(nDiffEVF = n_useDEevf, logFC_theoretical = logFC_theoretical))
}

set_mssc_meta <- function(symsim, phyla, logfc_threshold = 0.6) {
  ## dea: DE gene analysis results
  dea <- getDEgenes(symsim$true, popA = phyla$pop$sub1, popB = phyla$pop$sub2)

  ## TODO: set several criterias
  diffg <- which((dea$logFC_theoretical >= logfc_threshold) & (dea$nDiffEVF > 0))
  nondiffg <- setdiff(seq_len(length(dea$logFC_theoretical)), diffg)

  ind <- symsim$umi$cell_meta$pop
  
  ncell <- length(ind)
  cond <- floor(rep(0, ncell))
  cond[ind %in% phyla$pop$sub1] <- 1
  cond[ind %in% phyla$pop$sub2] <- 2

  npop <- length(phyla$tip.label)
  cond_of_ind <- floor(rep(1, npop))
  cond_of_ind[phyla$pop$sub1] <- 1
  cond_of_ind[phyla$pop$sub2] <- 2
  return(invisible(list(dea = dea,
                        diffg = diffg,
                        nondiffg = nondiffg,
                        logfc_threshold = logfc_threshold,
                        ind = ind,
                        cond = cond,
                        cond_of_ind = cond_of_ind)))
}


## plot utilities
plot_tree <- function(p) {
  ape::plot.phylo(p, show.tip.label = F, lwd = 2)
  ape::nodelabels(cex = 1)
  ape::tiplabels(cex = 2)
  ape::edgelabels(text = sprintf("%0.2f", p$edge.length),
    col = "black", bg = "lightgreen", font = 1, adj = c(0.5, 1.5))
}

plot_tsne <- function(symsim_umi, color_values) {
  ## cnt should be ngene by ncell
  cnt <- symsim_umi$counts
  tpm <- scale(cnt, center = F, scale = colSums(cnt)) * 10000
  logtpm <- log(tpm + 1)
  types <- symsim_umi$cell_meta$pop
  tsne_logtpm <- Rtsne::Rtsne(t(logtpm), pca_scale = T, pca_center = T, initial_dims = 50,
    pca = T, check_duplicates = F)
  dat <- data.frame(cbind(tsne_logtpm$Y), types)
  colnames(dat) <- c("TSNE1", "TSNE2", "Type")
  p <- ggplot(dat, aes(x = TSNE1, y = TSNE2, color = factor(types))) +
    geom_point() +
    scale_colour_manual(values = color_values) +
    theme_bw()
  return(invisible(p))
}

plotviolin <- function(cnt, ind, genes, logtpm = F,
                       myggtitle = theme(
                         plot.title = element_text(
                           size = 15, hjust = 0.5))) {
  n <- length(genes)
  gnm <- paste0("gene", genes)
  select_cnt <- cnt[genes, ]
  if (logtpm) {
    select_cnt <- log(scale(select_cnt, center = F, scale = colSums(cnt)) * 10000 + 1)
  }
  total_ind <- max(ind)
  plotdata <- as.data.frame(t(select_cnt))
  colnames(plotdata) <- gnm

  ## use factor to order the plot
  ## in the order ind1, ind2, ..., not ind1. ind10, ind2, ...
  plotdata$pop <- factor(paste0("ind", ind),
    levels = paste0("ind", seq_len(total_ind)),
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

plot_genes_after_batcheffect <- function(symsim_umi,
                                         ind,
                                         diffg,
                                         nondiffg,
                                         nde = 40, nnde = 40,
                                         pnrow = 5,
                                         save_path = NULL,
                                         save_figure = TRUE,
                                         fnm_prefix = "symsim_violin",
                                         width = 20, height = 10) {
  ## show batch effect using violin plot for different genes.
  ## nde: num of differential genes to show
  ## nnde: num of non-differential genes to show

  ## violinplot for diffg
  pviolin_symsimdeg <- plotviolin(symsim_umi$counts, ind, diffg)
  n <- ifelse(nde > length(diffg), length(diffg), nde)
  sampled_pvd <- ggarrange(
    plotlist = pviolin_symsimdeg[sample(length(diffg), n, replace = F)],
    nrow = pnrow,
    ncol = ceiling(nde / pnrow)
  )

  ## violinplot for nondiffg
  pviolin_symsim_sndeg <- plotviolin(symsim_umi$counts, ind, nondiffg)
  m <- ifelse(nnde > length(nondiffg), length(nondiffg), nnde)
  sampled_pvsnd <- ggarrange(
    plotlist = pviolin_symsim_sndeg[sample(length(nondiffg), m, replace = F)],
    nrow = pnrow,
    ncol = ceiling(nnde / pnrow)
  )

  if (save_figure) {
    ggsave(file.path(save_path,
      stringr::str_glue(fnm_prefix, "dg.pdf", .sep = "_")),
    plot = sampled_pvd, width = width, height =  height)
    ggsave(file.path(save_path,
      stringr::str_glue(fnm_prefix, "nondg.pdf", .sep = "_")),
    plot = sampled_pvsnd, width = width, height = height)
  }

  invisible(list(
    pvln_alldegs = pviolin_symsimdeg,
    pvln_allndegs = pviolin_symsim_sndeg,
    spvln_degs = sampled_pvd,
    spvln_ndegs = sampled_pvsnd
  ))
}

## * main
phyla3 <- set_population_struct_using_phylo(w = c(0.2, 0.5, 0.4, 0.4, 0.5, 0.3))
plot_tree(phyla3)
symsim <- get_symsim_umi(phyla = phyla3, bimod = 1, capt_alpha = 0.1,
  ngene = 100, ncell_per_subpop = 60)
color_values <- c("brown", "brown2", "brown3",
  "chartreuse", "chartreuse3", "darkolivegreen1")
p_tsne <- plot_tsne(symsim_umi = symsim$umi, color_values)
p_tsne

## de analysis
mssc_meta <- set_mssc_meta(symsim = symsim, phyla = phyla3)
p <- plot_genes_after_batcheffect(symsim_umi = symsim$umi,
                                  ind = mssc_meta$ind,
                                  diffg = mssc_meta$diffg,
                                  nondiffg = mssc_meta$nondiffg,
                                  nde = 20, nnde = 20,
                                  save_figure = F)
p$spvln_degs
p$spvln_ndegs

## violin
genes <- sample(seq_len(100), size = 50, replace = F)
p_violins <- plotviolin(cnt = symsim$umi$counts, ind = symsim$umi$cell_meta$pop,
  genes = genes)
p <- ggarrange(plotlist = p_violins, nrow = 10, ncol = 5)
p



## logtpm
## p_violins <- plotviolin(cnt = symsim$umi$counts, ind = symsim$umi$cell_meta$pop,
##                         logtpm = T,
##            genes = genes)
## p <- ggarrange(plotlist = p_violins, nrow = 10, ncol = 5)
## p
