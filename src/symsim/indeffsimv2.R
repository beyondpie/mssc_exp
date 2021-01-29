## simulate individual effect (bf) v2 using SymSim

## * set R environment
suppressWarnings(suppressMessages({
  library(SymSim)
  library(tidyverse)
  library(grid)
  library(gtable)
  library(gridExtra)
  library(ggpubr)
  library(ape)
  library(ggtree)
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

## * meta
color3 <- c("brown", "brown2", "brown3",
            "chartreuse", "chartreuse3", "darkolivegreen1")

color5 <- c("brown", "brown2", "brown3", "brown4", "deeppink4",
            "chartreuse", "chartreuse3", "darkolivegreen1", "darkolivegreen3", "aquamarine3")

color10 <- c("chocolate", paste0("chocolate", 1:4),
             "coral", paste0("coral", 1:4),
             "deepskyblue", paste0("deepskyblue", 1:4),
             "darkslategray", paste0("darkslategray", 1:4))

## * help functions
logtpm <- function(cnt, scale = 10000) {
  tpm <- scale(cnt, center = F, scale = colSums(cnt)) * scale
  return(invisible(log(tpm + 1)))
}
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

get_symsim_de_analysis <- function(true_counts_res, popA, popB) {
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

  return(list(nDiffEVF = n_useDEevf, logFC_theoretical = logFC_theoretical))
}

set_mssc_meta <- function(symsim_umi, phyla) {
  ind <- symsim_umi$cell_meta$pop
  sub1 <- phyla$pop$sub1
  sub2 <- phyla$pop$sub2

  ncell <- length(ind)
  cond <- floor(rep(0, ncell))
  cond[ind %in% sub1] <- 1
  cond[ind %in% sub2] <- 2

  npop <- length(sub1) + length(sub2)
  cond_of_ind <- floor(rep(1, npop))
  cond_of_ind[sub1] <- 1
  cond_of_ind[sub2] <- 2
  return(invisible(list(ind = ind,
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

plot_tree_v2 <- function(phyla) {
  ## use ggtree package to plot tree in ggplot way
  p <- ggtree(phyla) + theme_tree2() +
    geom_nodepoint(color = "green", size = 5) +
    geom_rootpoint(color = "red", size = 6) +
    geom_tippoint(color = "blue", size = 4) +
    geom_tiplab(align = TRUE, linetype = 'dashed', linesize = .3, hjust =  -1) +
    geom_treescale(fontsize=5, linesize=2) +
    theme(axis.text = element_text(size = 20))
  return(invisible(p))
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

violin_plot <- function(cnt, ind, genes, logtpm = F,
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

plot_de_violin <- function(symsim_umi,
                           ind,
                           diffg,
                           nondiffg,
                           logfc,
                           nde = 20, nnde = 20,
                           pnrow = 5) {
  ## show batch effect using violin plot for different genes.
  ## nde: num of differential genes to show
  ## nnde: num of non-differential genes to show

  ## violinplot for diffg
  pviolin_symsimdeg <- violin_plot(symsim_umi$counts, ind, diffg)
  n <- ifelse(nde > length(diffg), length(diffg), nde)
  pvsd <- arrangeGrob(
    grobs = pviolin_symsimdeg[sample(length(diffg), n, replace = F)],
    nrow = pnrow,
    ncol = ceiling(nde / pnrow),
    top = paste0("DE under logFC >= ", logfc)
  )

  ## violinplot for nondiffg
  pviolin_symsim_sndeg <- violin_plot(symsim_umi$counts, ind, nondiffg)
  m <- ifelse(nnde > length(nondiffg), length(nondiffg), nnde)
  pvsnd <- arrangeGrob(
    grobs = pviolin_symsim_sndeg[sample(length(nondiffg), m, replace = F)],
    nrow = pnrow,
    ncol = ceiling(nnde / pnrow),
    top = paste0("Non-DE under logFC < ", logfc)
  )
  return(invisible(list(pvsd, pvsnd)))
}

simu <- function(ncell_per_ind = 30, nind_per_cond = 3, brn_len = 0.5,
                 bimod = 0, sigma = 0.6, capt_alpha = 0.2, ngene = 200,
                 seed = 1L) {
  ## simulate individual effect using symsim
  ncond <- 2
  nind <- nind_per_cond * ncond
  ncell <- nind * ncell_per_ind
  w <- rep(brn_len, nind)
  phyla <- set_population_struct_using_phylo(w = brn)
  symsim <- get_symsim_umi(phyla = phyla, bimod = bimod,
    capt_alpha = capt_alpha, gene_module_prop = 0.0,
    ncell_per_subpop = ncell_per_ind,
    sigma = sigma,
    ngene = ngene,
    seed = seed,
    nevf = 10,
    n_de_evf = 7,
    vary = "s")
  dea <- get_symsim_de_analysis(true_counts_res = symsim$true,
    popA = phyla$pop$sub1,
    popB = phyla$pop$sub2)
  return(invisible(list(phyla = phyla,
    symsim = symsim,
    dea = dea)))
}

load_mssc <- function(nind = 10, mssc_version = "mssc_2-0",
                      tol_rel_obj = 1e-06,
                      num_iter = 20000, output_samples = 3000) {
  ## nind: total number individuals
  ## tol_rel_obj, num_iter, output_samples are used
  ## to tune the model
  ## return:
  ## - a mssc model
  mssc_path <- here::here("src", "mssc")
  model <- modules::import_(mssc_version)
  invisible(model$High2$new(
    stan_snb_path = file.path(mssc_path, "stan", "snb.stan"),
    stan_high2_path = file.path(mssc_path, "stan",
                                paste0(mssc_version, ".stan")),
    nind = nind,
    tol_rel_obj = tol_rel_obj,
    algorithm = "meanfield",
    adapt_engaged = FALSE,
    eta = 0.1,
    ## or iter more steps
    num_iter = num_iter,
    output_samples = output_samples
  ))
}

set_result_array <- function(rpt = 5, ncells = c(20, 40, 80),
                             methods = c("mssc_2-0", "pseudo", "wilcox", "t")) {
  ## setup the result array:
  ## - num_measurement by length(ncell) by repeat_num
  invisible(array(data = NA, dim = c(length(methods), length(ncells), rpt),
                  dimnames = list(methods, ncells, seq_len(rpt))))
}

## * main
## de analysis
mssc_meta <- set_mssc_meta(symsim = symsim, phyla = phyla3)
p <- plot_de_violin(symsim_umi = symsim$umi,
  ind = mssc_meta$ind,
  diffg = mssc_meta$diffg,
  nondiffg = mssc_meta$nondiffg,
  nde = 20, nnde = 20,
  save_figure = F)
p$spvln_degs
p$spvln_ndegs

main <- function(nind_per_cond,
                 brn_len,
                 bimod,
                 sigma,
                 ncells,
                 capt_alpha,
                 ngene = 100,
                 rpt = 5,
                 save_figure = TRUE,
                 fig_width = 20,
                 fig_height = 10,
                 nde_plt = 20,
                 nnde_plt = 20,
                 pnrow = 5,
                 save_symsim_data = FALSE,
                 save_mssc_model = FALSE) {
  ## argument:
  ## - ncells: vector of ncell_per_ind
  
  ## set local path to save results
  result_dir <- here::here("src", "symsim",
                           paste0("symsim_", format(Sys.time(), format = "%Y%m%d%H%M")))
  simu_fig_dir <- file.path(result_dir, "simu_figs")
  simu_data_dir <- file.path(result_dir, "simu_data")
  dea_dir <- file.path(result_dir, "dea")
  mssc_v2_dir <- file.path(result_dir, "mssc_v2")

  for ( d in c(result_dir, simu_fig_dir, simu_data_dir, dea_dir, mssc_v2_dir)) {
    if (!dir.exists(d)) {
      dir.create(path = d)
    }
  }
  
  ## load mssc model
  nind_all <- nind_per_cond * 2
  mssc_20 <- load_mssc(nind = nind_all, mssc_version = "mssc_2-0")

  ## declare the result
  auc06 <- set_result_array(rpt = rpt, ncells = ncells,
                        methods = c("mssc_2-0", "pseudo", "wilcox", "t"))
  auc08 <- set_result_array(rpt = rpt, ncells = ncells,
                        methods = c("mssc_2-0", "pseudo", "wilcox", "t"))
  auc10 <- set_result_array(rpt = rpt, ncells = ncells,
                        methods = c("mssc_2-0", "pseudo", "wilcox", "t"))

  for (i in seq_len(rpt)) {
    auc06_i <- matrix(NA, nrow = dim(r)[1], ncol = dim(r)[2])
    auc08_i <- matrix(NA, nrow = dim(r)[1], ncol = dim(r)[2])
    auc10_i <- matrix(NA, nrow = dim(r)[1], ncol = dim(r)[2])
    rownames(auc06_i) <- rownames(auc06)
    rownames(auc08_i) <- rownames(auc08)
    rownames(auc10_i) <- rownames(auc10)
    for (j in seq_along(ncells)) {
      ncell <- ncells[j]
      ## simulate data
      symsim_prefix <- str_glue("{ngene}gene", "{nind_all}ind",
                                  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
                                  "{sigma}sigma", "{capt_alpha}alpha", "{i}seed",
                                "{rpt}rpt", .sep = "_")
      message(str_glue("SymSim experiment: {symsim_prefix}."))
      mysimu <- simu(ncell_per_ind = ncell,
                     nind_per_cond = nind_per_cond,
                     brn_len = brn_len,
                     bimod = bimod,
                     sigma = sigma,
                     capt_alpha = capt_alpha,
                     ngene = ngene,
                     seed = i)
      ## get differentially expressed genes based different fold change levels
      diffg_06 <- which((mysimu$dea$logFC_theoretical >= 0.6) & (mysimu$dea$nDiffEVF >  0))
      nondiffg_06 <- setdiff(seq_along(mysimu$dea$logFC_theoretical), diffg_06)
      message(str_glue("SymSim DE analysis under logFC 0.6: ",
                       "{length(diffg_06)} diffg and {length(nondiffg_06)} nondiffg",
                       .sep = "\n"))

      diffg_08 <- which((mysimu$dea$logFC_theoretical >= 0.8) & (mysimu$dea$nDiffEVF >  0))
      nondiffg_08 <- setdiff(seq_along(mysimu$dea$logFC_theoretical), diffg_08)
      message(str_glue("SymSim DE analysis under logFC 0.8: ",
                       "{length(diffg_08)} diffg and {length(nondiffg_08)} nondiffg",
                       .sep = "\n"))

      diffg_10 <- which((mysimu$dea$logFC_theoretical >= 1.0) & (mysimu$dea$nDiffEVF >  0))
      nondiffg_10 <- setdiff(seq_along(mysimu$dea$logFC_theoretical), diffg_10)
      message(str_glue("SymSim DE analysis under logFC 1.0: ",
                       "{length(diffg_10)} diffg and {length(nondiffg_10)} nondiffg",
                       .sep = "\n"))

      ## prepare mssc meta
      mssc_meta <- set_mssc_meta(symsimumi = mysimu$symsim$umi, phyla = mysimu$phyla)
    
      ## save data
      if (save_figure) {
        pdf(file = file.path(simu_fig_dir, paste0(symsim_prefix, ".pdf")),
            width = fig_width, height = fig_height)
        ## phylo tree
        p_phylo <- plot_tree_v2(phyla = mysimu$phyla)
        colors <- ifelse(test = (nind_per_cond == 3), yes = color3,
                         no = ifelse(test = (nind_per_cond == 5),
                                    yes = color5,
                                    no = color10))
        ## tsne of cells under population
        p_tsne <- plot_tsne(symsim_umi = mysimu$symsim$umi, color_values = colors)
        grid.arrange(grobs = list(p_phylo, p_tsne), nrow  = 1, ncol = 2,
                     top = "Population structure and t-SNE in SymSim")
        ## violin of (non-)differentially expression genes.
        pv_06 <- plot_de_violin(symsim_umi = mysimu$symsim$umi,
                                ind = mssc_meta$ind,
                                diffg = diffg_06, nondiffg = nondiffg_06,
                                nde = nde_plt, nnde = nnde_plt, pnrow = pnrow,
                                logFC = 0.6)
        print(pv_06[[1]])
        print(pv_06[[2]])
        pv_08 <- plot_de_violin(symsim_umi = mysimu$symsim$umi,
                                ind = mssc_meta$ind,
                                diffg = diffg_08, nondiffg = nondiffg_08,
                                nde = nde_plt, nnde = nnde_plt, pnrow = pnrow,
                                logFC = 0.8)
        print(pv_08[[1]])
        print(pv_08[[2]])
        pv_10 <- plot_de_violin(symsim_umi = mysimu$symsim$umi,
                                ind = mssc_meta$ind,
                                diffg = diffg_10, nondiffg = nondiffg_10,
                                nde = nde_plt, nnde = nnde_plt, pnrow = pnrow,
                                logFC = 1.0)
        print(pv_10[[1]])
        print(pv_10[[2]])
        dev.off()
      } ## end of save figures in one file
      if (save_symsim_data) {
        saveRDS(object = list(simu = mysimu, dg06 = diffg_06, nondg06 = nondiffg_06,
                              dg08 = diffg_08, nondg08 = nondiffg_08,
                              dg10 = diffg_10, nondg10 = nondiffg_10),
                file = file.path(simu_data_dir, paste0(symsim_prefix, ".rds")))
      } ## end of save symsim data

      ## de analysis with mssc
      r_mssc20 <- run_mssc(
        model = mssc_20,
        symsim_umi = mysimu$symsim$umi$counts,
        mssc_meta = mssc_meta,
        save_result = save_mssc_model,
        save_path = file.path(mssc_v2_dir, str_glue("mssc2_{symsim_prefix}.rds"))
      )
      auc06_mssc20 <- mssc_20$get_auc(r_mssc20$r, c1  = diffg_06, c2 = nondiffg_06)
      auc08_mssc20 <- mssc_20$get_auc(r_mssc20$r, c1  = diffg_08, c2 = nondiffg_08)
      auc10_mssc20 <- mssc_20$get_auc(r_mssc20$r, c1  = diffg_10, c2 = nondiffg_10)
      ## de analysis with pseudobulk
      r_pseudo <- mypseudo$pseudobulk_deseq2(
        cnt_gbc = mysimu$symsim$umi$counts,
        mybatches = mssc_meta$ind,
        myconds = factor(mssc_meta$cond)
      )
      auc06_pseudo <- mypseudo$calc_auc(deseq2_res = r_pseudo,
                                        degs = diffg_06, ndegs = nondiffg_06)
      auc08_pseudo <- mypseudo$calc_auc(deseq2_res = r_pseudo,
                                        degs = diffg_08, ndegs = nondiffg_08)
      auc10_pseudo <- mypseudo$calc_auc(deseq2_res = r_pseudo,
                                        degs = diffg_10, ndegs = nondiffg_10)
      
      ## de analysis with t-test
      r_t <- apply(t.test())
      ## de analysis with wilcox

      ## save result
      auc06_i[, j] <- c(auc06_mssc20, auc06_pseudo)
      auc08_i[, j] <- c(auc08_mssc20, auc08_pseudo)
      auc10_i[, j] <- c(auc10_mssc20, auc10_pseudo)
    } ## end of ncells
    auc06[, , i] <- auc06_i
    auc08[, , i] <- auc08_i
    auc10[, , i] <- auc10_i
    message(str_glue("In {i}th turn, result is summarized below..."))
    message("when logfc is 0.6:")
    print(auc06)
    message("when logfc is 0.8:")
    print(auc08)
    message("when logfc is 1.0:")
    print(auc10)
    ## save result
    saveRDS(object = list(r06 = auc06, r08 = auc08, r10 = auc10),
            file = file.path(dea_dir, str_glue("dea_r_{symsim_prefix}.rds")))
  } ## end of rpt
}

library(optparse)
option_list <- list(
  ## make_option(c("--ncell_per_ind"), action = "store", type = "integer", default = 30),
  ## make_option(c("--ngene"), action = "store", type = "integer", default = 100),
  make_option(c("--nind_per_cond"), action = "store", type = "integer", default = 3),
  make_option(c("--brn_len"), action = "store", type = "double", default = 0.5),
  make_option(c("--bimod"), action = "store", type = "double", default = 0.1),
  make_option(c("--sigma"), action = "store", type = "double", default = 0.6),
  make_option(c("--capt_alpha"), action = "store", type = "double", default = 0.2)
)

args <- parser_args(OptionParser(option_list = option_list))
main(
  nind_per_cond = args$nind_per_cond,
  brn_len = args$brn_len,
  bimod = args$bimod,
  sigma = args$sigma,
  ncells = c(30, 50, 80, 120, 160, 240, 300),
  rpt = 5,
  save_figure = T,
  fid_width = 20,
  fig_height = 10,
  save_symsim_data = FALSE,
  save_mssc_model = FALSE
)

