## simulate individual effect using SymSim.

## * set R environment

## How to install SymSim
## install.packages("devtools")
## library(devtools)
## devtools::install_github("JunLiuLab/SIMPLEs",
##   ref = "master",
##   build_vignettes = TRUE
## )

## Install package here to locate the root of the project
## install.packages("here")

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
  library(optparse)
}))

## warnings/errors traceback settings
options(error = traceback)
options(warn = -1)
## * meta
color3 <- c("brown", "brown2", "brown3",
  "chartreuse", "chartreuse3", "darkolivegreen1")

color5 <- c("brown", "brown2", "brown3", "brown4", "deeppink4",
  "chartreuse", "chartreuse3", "darkolivegreen1", "darkolivegreen3", "aquamarine3")

color10 <- c("chocolate", paste0("chocolate", 1:4),
  "coral", paste0("coral", 1:4),
  "deepskyblue", paste0("deepskyblue", 1:4),
  "darkslategray", paste0("darkslategray", 1:4))

## * plot functions
plot_tree_v2 <- function(phyla) {
  ## use ggtree package to plot tree in ggplot way
  p <- ggtree(phyla) + theme_tree2() +
    geom_nodepoint(color = "green", size = 5) +
    geom_rootpoint(color = "red", size = 6) +
    geom_tippoint(color = "blue", size = 4) +
    geom_tiplab(align = TRUE, linetype = 'dashed', linesize = .3, hjust =  -1) +
    geom_treescale(fontsize = 5, linesize = 2) +
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

run_symsim_simu <- function(phyla, bimod = 0, capt_alpha = 0.1,
                            gene_module_prop = 0.0,
                            ncell_per_ind = 100,
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
  ncell_total <- npop * ncell_per_ind
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

get_symsim_simu <- function(ncell_per_ind = 300, nind_per_cond = 3, brn_len = 0.5,
                            bimod = 0, sigma = 0.6, capt_alpha = 0.2, ngene = 200,
                            seed = 1, logfc_threshold = 0.8, ndg_threshold = 10, ntry = 10,
                            save_data = TRUE, save_data_path = "symsim.rds",
                            save_figure = FALSE, save_fig_path = "symsim.pdf",
                            fig_width = 20, fig_height = 10,
                            nde_plt = 20, nnde_plt = 20, pnrow = 5) {
  ## simulate individual effect using symsim
  ## get differentially expressed genes based different fold change levels
  ## 0.8 is a recommended value (0.6 or 1.0 is ok, but 1.0 seems too restricted.)
  ## the program could throw error.

  ncond <- 2
  nind <- nind_per_cond * ncond
  phyla <- set_population_struct_using_phylo(w = rep(brn_len, nind))

  symsim <- NULL
  dea <- NULL
  ndiffg <- 0
  ntry_now <- 1
  ## run symsim until it touches the needs.
  ## usually it will run only one time.
  while ((ndiffg < ndg_threshold) && (ntry_now < ntry)) {
    message(str_glue(""))
    symsim <- run_symsim_simu(phyla = phyla, bimod = bimod,
      capt_alpha = capt_alpha, gene_module_prop = 0.0,
      ncell_per_ind = ncell_per_ind,
      sigma = sigma,
      ngene = ngene,
      seed = seed,
      nevf = 10,
      n_de_evf = 7,
      vary = "s")
    dea <- get_symsim_de_analysis(true_counts_res = symsim$true,
      popA = phyla$pop$sub1,
      popB = phyla$pop$sub2)
    diffg <- which((dea$logFC_theoretical >= logfc_threshold) & (dea$nDiffEVF > 0))
    nondiffg <- setdiff(seq_along(dea$logFC_theoretical), diffg)
    ndiffg <- length(diffg)
    ntry_now  <- ntry_now + 1
  } ## end of while

  if ((ntry_now == ntry) && (ndiffg < ndg_threshold)) {
    stop(str_glue("Simulation failed: ndiffg_{ndifg} after try {ntry} times."))
  }
  r <- list(phyla = phyla, symsim = symsim, dea = dea,
    diffg = diffg, nondiffg = nondiffg,
    logfc_threshold = logfc_threshold)
  if (save_data) {
    saveRDS(object = r, file = save_data_path)
  } ## end of save data

  if (save_figure) {
    pdf(file = save_fig_path, width = fig_width, height = fig_height)
    ## phylo tree
    p_phylo <- plot_tree_v2(phyla)
    if (nind_per_cond == 3) {
      colors <- color3
    } else if (nind_per_cond == 5) {
      colors <- color5
    } else {
      colors <- color10
    }
    ## tsne of cells under population
    p_tsne <- plot_tsne(symsim_umi = symsim$umi, color_values = colors)
    grid.arrange(grobs = list(p_phylo, p_tsne), nrow  = 1, ncol = 2,
      top = "Population structure and t-SNE in SymSim")

    tryCatch(
      {
        ## violin of (non-)differentially expression genes.
        ## may have bugs when nde / nnde equals to 1
        pv_08 <- plot_de_violin(symsim_umi = symsim$umi,
          ind = symsim$umi$cell_meta$pop,
          diffg = diffg, nondiffg = nondiffg,
          nde = nde_plt, nnde = nnde_plt, pnrow = pnrow,
          logfc = logfc_threshold)
        grid.arrange(pv_08[[1]])
        grid.arrange(pv_08[[2]])
      },
      error = function(cond) {
        message(cond)
      },
      finally = {
        dev.off()
      }
    )
  } ## end of save figure.
  return(invisible(r))
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



## * main
## number of individuals in each condition
nind_per_cond <- 5
## branch length towards the root for two condiitons
## the longer the branch length, the more different the genes in each condition.
brn_len <- 0.5
bimod <- 1
## roughly the variances of the gene expressions within one individual
sigma <- 0.6
## number of cells in each individual
ncells  <-  c(100)
## capture rate of RNA molecules in scRNA-seq experiment.
capt_alpha <- 0.2
## number of genes in the simulation
ngene <- 200
## simulation repeat index
rpt <- 1
## de threshold used in symsim
logfc_threshold <- 0.8 

## where to save the figure
## By default, create a local dir
simu_data_dir <- here::here("src", "symsim",
                         paste0("simu_data_", format(Sys.time(), format = "%Y%m%d")))
if (!dir.exists(simu_data_dir)) {
  dir.create(path = simu_data_dir)
}

## simulation process
nind_all <- nind_per_cond * 2
symsim_prefix <- str_glue(
  "{ngene}gene", "{nind_all}ind", "{ncells}cell",
  "{brn_len}w", "{bimod}bimod", "{sigma}sigma", "{capt_alpha}alpha", .sep = "_")
simubulk <- get_symsim_simu(ncell_per_ind = 300, nind_per_cond = nind_per_cond,
  brn_len = brn_len, bimod = bimod, sigma = sigma,
  capt_alpha = capt_alpha, ngene = ngene,
  seed = 1, logfc_threshold = logfc_threshold, ndg_threshold = 10,
  ntry = 10, save_data = TRUE,
  save_data_path = file.path(simu_data_dir,
    paste0(symsim_prefix, ".rds")),
  save_figure = TRUE,
  save_fig_path = file.path(simu_data_dir,
    paste0(symsim_prefix, ".pdf")),
  fig_width = 20, fig_height = 10, nde_plt = 20,
  nnde_plt = 20, pnrow = 5)

## we can get the differential genes in the following way.
## diffg <- simubulk$diffg
## nondiffg <- simubulk$nondiffg


                           

