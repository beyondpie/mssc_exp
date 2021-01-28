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

set_mssc_meta <- function(symsim, phyla, logfc_threshold = 0.6) {
  ## dea: DE gene analysis results
  dea <- get_symsim_de_analysis(symsim$true, popA = phyla$pop$sub1, popB = phyla$pop$sub2)

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
  pviolin_symsimdeg <- violin_plot(symsim_umi$counts, ind, diffg)
  n <- ifelse(nde > length(diffg), length(diffg), nde)
  sampled_pvd <- ggarrange(
    plotlist = pviolin_symsimdeg[sample(length(diffg), n, replace = F)],
    nrow = pnrow,
    ncol = ceiling(nde / pnrow)
  )

  ## violinplot for nondiffg
  pviolin_symsim_sndeg <- violin_plot(symsim_umi$counts, ind, nondiffg)
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
  r <- set_result_array(rpt = rpt, ncells = ncells,
                        methods = c("mssc_2-0", "pseudo", "wilcox", "t"))

  for (i in seq_len(rpt)) {
    r_i <- matrix(NA, nrow = dim(r)[1], ncol = dim(r)[2])
    rownames(r_i) <- rownames(r)
    for (j in seq_along(ncells)) {
      ncell <- ncells[j]
      ## simulate data
      symsim_data_fnm <- str_glue("{ngene}gene", "{nind_all}ind",
                                  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
                                  "{sigma}sigma", "{capt_alpha}alpha", "{i}seed",
                                  "{rpt}rpt.rds", .sep = "_")
      mysimu <- simu(ncell_per_ind = ncell,
                     nind_per_cond = nind_per_cond,
                     brn_len = brn_len,
                     bimod = bimod,
                     sigma = sigma,
                     capt_alpha = capt_alpha,
                     ngene = ngene,
                     seed = i)
    } ## end of ncells
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
