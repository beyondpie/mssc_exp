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
  library(optparse)
}))

## warnings/errors traceback settings
options(error = traceback)
options(warn = -1)
options(mc.cores = 3)

## load help functions
options("import.path" = list(here::here("rutils"),
  here::here("src", "mssc")))

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

## * pseudobulk analysis
## functions for pseudobulk analysis

get_pseudobulk <- function(cnt_gbc, mybatches) {
  ## given the cnt data, and colmeta batches,
  ## return the pseudobulk data with the colnames.

  ubatches <- sort(unique(mybatches))
  pseudobulk <- vapply(ubatches, function(i) {
    rowSums(cnt_gbc[, mybatches %in% i])
  }, FUN.VALUE = rep(0, nrow(cnt_gbc)))
  colnames(pseudobulk) <- as.character(ubatches)
  invisible(pseudobulk)
}


pseudobulk_deseq2 <- function(cnt_gbc,
                              mybatches,
                              myconds) {
  ## using deseq2 to analyze pseudobulk
  ## return the data.frame format of deseq result.

  mypseudobulk <- get_pseudobulk(cnt_gbc, mybatches)
  names(myconds) <- as.character(mybatches)
  ubatches <- colnames(mypseudobulk)
  uconds <- myconds[as.character(ubatches)]

  coldf <- data.frame(ubatches, uconds)
  exp_design <- ~ as.factor(uconds)

  dataset <- DESeq2::DESeqDataSetFromMatrix(
    countData = mypseudobulk,
    colData = coldf,
    design = exp_design
  )
  r <- data.frame(DESeq2::results(DESeq2::DESeq(dataset)))
  invisible(r)
}

calc_auc <- function(deseq2_res, degs, ndegs,
                     scorecol = "pvalue") {
  mybackend <- c(rep(TRUE, length(degs)), rep(FALSE, length(ndegs)))
  bgnms <- c(degs, ndegs)
  names(mybackend) <- as.character(bgnms)

  scores <- deseq2_res[[scorecol]]
  names(scores) <- rownames(deseq2_res)
  scores[is.na(scores)] <- 1.0
  myauc <- caTools::colAUC(scores[bgnms], mybackend)
  invisible(list(auc = myauc, sts = scores))
}


## * help functions
get_logtpm <- function(cnt, scale = 10000) {
  tpm <- scale(cnt, center = FALSE, scale = colSums(cnt)) * scale
  return(invisible(log(tpm + 1)))
}

zhu_test <- function(x, group, test = "t") {
  if (sd(x) < 1e-06) {
    1
  } else {
    if (test == "t") {
      t.test(x ~ group)$p.value
    } else if (test == "wilcox") {
      wilcox.test(x ~ group)$p.value
    }
  }
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

get_symsim_simu <- function(ncell_per_ind = 300, nind_per_cond = 3, brn_len = 0.5,
                            bimod = 0, sigma = 0.6, capt_alpha = 0.2, ngene = 200,
                            seed = 1, logfc_threshold = 0.8, ndg_threshold = 5, ntry = 5,
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
    stop(str_glue("Simulation failed:",
      "ndiffg is {ndiffg}",
      "less than threshold {ndg_threshold}",
      "after try {ntry} times."))
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

get_symsim_by_sampling <- function(simubulk,
                                   ncell_per_ind = 20,
                                   seed = 1L) {
  ## generate symsim simulation data from another symsim data
  ## by sampling the specific number of cells from each individual.

  ## Return
  ## - another symsim simulation data (list) like mysymsim
  ##   plus the sample_cell_index as one field

  ## reproduce the sampling results
  set.seed(seed = seed)
  ind <- simubulk$symsim$umi$cell_meta$pop
  sample_cell_index <- unlist(lapply(seq_len(max(ind)), function(i) {
    cell_index <- which(ind == i)
    ncell <- length(cell_index)
    if (ncell_per_ind >= ncell) {
      ## when ncell_per_ind is larger than ncell, we directly use
      ## all the cells for the individual
      return(invisible(cell_index))
    }
    r <- sample(x = cell_index, size = ncell_per_ind, replace = FALSE)
    return(invisible(r))
  }))

  r <- list(phyla = simubulk$phyla,
    symsim = list(true = NULL, umi = NULL),
    dea = simubulk$dea,
    diffg = simubulk$diffg,
    nondiffg = simubulk$nondiffg,
    logfc_threshold = simubulk$logfc_threshold,
    sample_cell_index = sample_cell_index)

  ## complete symsim true
  r$symsim$true <- list(
    counts = simubulk$symsim$true$counts[, sample_cell_index],
    ## cell_meta is a dataframe: ncell by features
    cell_meta = simubulk$symsim$true$cell_meta[sample_cell_index, ],
    gene_effects = simubulk$symsim$true$gene_effects,
    cell_meta = simubulk$symsim$true$cell_meta[sample_cell_index, ],
    ## kinetic_params: list of data.frame
    kinetic_params = lapply(simubulk$symsim$true$kinetic_params,
      function(x) {x[, sample_cell_index]}),
    in_module = simubulk$symsim$true$in_module)
  ## complete symsim umi
  r$symsim$umi <- list(
    counts = simubulk$symsim$umi$counts[, sample_cell_index],
    cell_meta = simubulk$symsim$umi$cell_meta[sample_cell_index, ],
    ## nreads_perUMI is a list (length of ncells)
    nreads_perUMI = simubulk$symsim$umi$nreads_perUMI[sample_cell_index],
    ## nUMI2seq is a vector length of ncells
    nUMI2seq = simubulk$symsim$umi$nUMI2seq[sample_cell_index])
  return(invisible(r))
}


load_mssc <- function(nind = 10, mssc_version = "mssc_2-0",
                      tol_rel_obj = 1e-06,
                      num_iter = 20000,
                      output_samples = 1000) {
  ## nind: total number individuals
  ## tol_rel_obj, num_iter, output_samples are used
  ## to tune the model
  ## return:
  ## - a mssc model
  mssc_path <- here::here("src", "mssc")
  module <- modules::import_(mssc_version)
  invisible(module$High2$new(
    stan_snb_path = file.path(mssc_path, "stan", "snb.stan"),
    stan_high2_path = file.path(mssc_path, "stan",
      paste0(mssc_version, ".stan")),
    stan_glm_path = file.path(mssc_path, "stan", "glm.stan"),
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

run_mssc <- function(model, symsim_umi, mssc_meta, diffg, nondiffg, save_result = TRUE,
                     save_path = NULL) {
  ## mssc analysis
  ## side effect
  ## - model state will be updated
  ## - save the results into rds file (save_path).
  ## return
  ## - a list of vi ranking statistics (ngene by 3 matrix)
  ##   and opt ranking statistics (ngene by 1 matrix)

  ## data from symsim
  y2c <- round(symsim_umi)
  ind <- mssc_meta$ind
  cond <- mssc_meta$cond
  sumcnt <- colSums(y2c)
  s <- sumcnt / median(sumcnt)

  ## model inference
  init_params <- model$init_params(
    cnt = y2c, s = s, cond = cond, ind = ind)
  model_data <- model$to_model_data(cnt = y2c,
                                    s = s, cond = cond, ind = ind, hp = init_params$hp)
  ## NOTE: NOT RUNNING VI
  ## model$run(data = model_data, list_wrap_ip = list(init_params$ip))
  model$run_opt(data = model_data, list_wrap_ip = list(init_params$ip))
  init_params_of_glm <- model$init_glm_params(
    cnt = y2c, s = s, cond = cond, ind = ind
  )
  model$run_glm_opt(data = model_data, list_wrap_ip = list(init_params_of_glm))

  ## get inference status
  ## mssc_vifit_state <- model$high2fit$return_codes()
  mssc_vifit_state <- 1
  mssc_optfit_state <- model$high2optfit$return_codes()
  glmfit_state <- model$glmoptfit$return_codes()

  ## get results
  ## Mannually get the samples for all the params
  ## the model, in principle, will only save the results
  ## in a temporary file since it uses cmdstan.
  if (mssc_vifit_state != 0) {
    warning(str_glue("mssc VI fit failed with codes: {mssc_vifit_state}"))
    vi_est_params <- NULL
    auc_vi <- -1
  } else {
    vi_est_params <- model$extract_draws_all(
      ngene = nrow(y2c), genenms = seq_len(nrow(y2c)), method = "vi")
    vi_mucond <- model$extract_draws(
      param = "mucond", ngene = nrow(y2c), genenms = seq_len(nrow(y2c)),
      method = "vi")
    ## three rankings in order: t, bf, m, we use the third
    vi_r <- model$get_ranking_statistics(
      mucond = vi_mucond, two_hot_vec = c(1, -1))
    auc_vi <- model$get_auc(vi_r, c1 = diffg, c2 = nondiffg)[3]
  }
  if (mssc_optfit_state != 0) {
    warning(str_glue("mssc OPT fit failed with codes: {mssc_optfit_state}"))
    opt_est_params <- NULL
    auc_opt <- -1
  } else {
    opt_est_params <- model$extract_draws_all(
      ngene = nrow(y2c), genenms = seq_len(nrow(y2c)), method = "opt"
    )
    opt_mucond <- model$extract_draws(
      param = "mucond", ngene = nrow(y2c), genenms = seq_len(nrow(y2c)),
      method = "opt")
    opt_r <- model$get_opt_ranking_statistic(
      mucond = opt_mucond, two_hot_vec = c(1, -1)
    )
    auc_opt <- model$get_auc(opt_r, c1 = diffg, c2 = nondiffg)
  }
  if (glmfit_state != 0) {
    warning(str_glue("GLM OPT fit failed with codes: {glmfit_state}"))
    auc_glm <- -1
  } else {
    glm_mucond <- model$extract_draws_from_glm(param = "mucond", ngene = nrow(y2c),
      genenms = seq_len(nrow(y2c)))
    glm_r <- model$get_opt_ranking_statistic(mucond = glm_mucond, two_hot_vec = c(1, -1))
    auc_glm <- model$get_auc(ranking_statistic = glm_r, c1 = diffg, c2 = nondiffg)
  }
  r <- list(auc_vi = auc_vi, auc_opt = auc_opt, auc_glm = auc_glm)
  if (save_result) {
    ## NOTE: check loading the results
    ## warning at opt:argparse
    saveRDS(object = list(vi_est_params = vi_est_params,
      opt_est_params = opt_est_params,
      r = r, model = model),
    file = save_path)
  }
  return(invisible(r))
}

set_result_array <- function(rpt = 5, ncells = c(20, 40, 80),
                             methods = c("mssc_2-0", "pseudo", "wilcox", "t")) {
  ## setup the result array:
  ## - num_measurement by length(ncell) by repeat_num
  invisible(array(data = NA, dim = c(length(methods), length(ncells), rpt),
    dimnames = list(methods, ncells, seq_len(rpt))))
}

## * main
main <- function(nind_per_cond,
                 brn_len,
                 bimod,
                 sigma,
                 ncells,
                 capt_alpha,
                 ngene = 200,
                 rpt = 5,
                 save_mssc_model = FALSE,
                 logfc_threshold = 0.8) {
  ## argument:
  ## - ncells: vector of ncell_per_ind
  ## set local path to save results

  ## result_dir <- here::here(
  ##   "src", "symsim",
  ##   ## use year-month-day to deside the dir
  ##   paste0("symsim_", format(Sys.time(), format = "%Y%m%d")))
  result_dir <- here::here(
    "src", "symsim", "OPT_symsim"
  )
  
  simu_data_dir <- file.path(result_dir, "simu_data")
  dea_dir <- file.path(result_dir, "dea")
  mssc_v2_dir <- file.path(result_dir, "mssc_v2")

  for (d in c(result_dir, simu_data_dir, dea_dir, mssc_v2_dir)) {
    if (!dir.exists(d)) {
      dir.create(path = d)
    }
  }

  ## - simulate datasets
  nind_all <- nind_per_cond * 2
  symsim_prefix <- str_glue(
    "{ngene}gene", "{nind_all}ind", "{brn_len}w", "{bimod}bimod",
    "{sigma}sigma", "{capt_alpha}alpha", .sep = "_")
  message(str_glue("SymSim experiment: ", "{symsim_prefix}.", .sep = "\n"))
  simubulk <- get_symsim_simu(ncell_per_ind = 300, nind_per_cond = nind_per_cond,
    brn_len = brn_len, bimod = bimod, sigma = sigma,
    capt_alpha = capt_alpha, ngene = ngene,
    seed = 1, logfc_threshold = logfc_threshold, ndg_threshold = 5,
    ntry = 10, save_data = TRUE,
    save_data_path = file.path(simu_data_dir,
      paste0(symsim_prefix, ".rds")),
    save_figure = TRUE,
    save_fig_path = file.path(simu_data_dir,
      paste0(symsim_prefix, ".pdf")),
    fig_width = 20, fig_height = 10, nde_plt = 20,
    nnde_plt = 20, pnrow = 5)
  diffg <- simubulk$diffg
  nondiffg <- simubulk$nondiffg
  ## report diffg info.
  message(str_glue("SymSim DE analysis under logFC {logfc_threshold}:",
    "{length(diffg)} diffg and {length(nondiffg)} nondiffg",
    .sep = "\n"))
  ## - load mssc model
  mssc_20 <- load_mssc(nind = nind_all, mssc_version = "mssc_2-0")

  ## - declare the result
  aucs <- set_result_array(rpt = rpt, ncells = ncells,
    methods = c("mssc_vi", "mssc_opt", "glm", "pseudo_deseq2_no_inds",
      "wilcox", "t"))

  ## - start experiment
  for (i in seq_len(rpt)) {
    auci <- matrix(NA, nrow = dim(aucs)[1], ncol = dim(aucs)[2])
    rownames(auci) <- rownames(aucs)
    for (j in seq_along(ncells)) {
      message(str_glue("when ncell per individual is {ncells[j]}:"))
      ncell <- ncells[j]
      ## simulate data
      mysimu <- get_symsim_by_sampling(simubulk, ncell_per_ind = ncell, seed = i)
      ## prepare mssc meta
      mssc_meta <- set_mssc_meta(symsim_umi = mysimu$symsim$umi, phyla = mysimu$phyla)

      tryCatch({
        ## - de analysis with mssc
        ## rarely fitting unfished, and unable to retrieve the draws.
        start_time <- Sys.time()
        r_mssc20 <- run_mssc(
          model = mssc_20,
          symsim_umi = mysimu$symsim$umi$counts,
          mssc_meta = mssc_meta,
          diffg = diffg,
          nondiffg = nondiffg,
          save_result = save_mssc_model,
          save_path = file.path(mssc_v2_dir, str_glue("mssc2_{symsim_prefix_per_rpt}.rds"))
        )
        end_time <- Sys.time()
        message(str_glue("mssc running time: ",
          "{format(end_time - start_time, nsmall = 2)}"))
        ## de analysis with pseudobulk
        ## TODO: consider the genes when p-value not NA
        ## rarely fitting might fail.
        r_pseudo_deseq2_no_inds <- pseudobulk_deseq2(
          cnt_gbc = mysimu$symsim$umi$counts,
          mybatches = mssc_meta$ind,
          myconds = factor(mssc_meta$cond)
        )

        ## de analysis with t-test
        logtpm <- get_logtpm(cnt = mysimu$symsim$umi$counts, scale = 10000)
        r_t <- apply(logtpm, 1, zhu_test, group = mssc_meta$cond, test = "t")
        r_t_adjp <- p.adjust(r_t, method = "fdr")
        ## de analysis with wilcox
        r_wilcox <- apply(logtpm, 1, zhu_test, group = mssc_meta$cond, test = "wilcox")
        r_wilcox_adjp <- p.adjust(r_wilcox, method = "fdr")

        auc_pseudo_deseq2_no_inds <- calc_auc(
          deseq2_res = r_pseudo_deseq2_no_inds, degs = diffg, ndegs = nondiffg)$auc
        auc_t <- caTools::colAUC(
          X = r_t_adjp,
          y = (seq_along(mysimu$dea$logFC_theoretical) %in% diffg))
        auc_wilcox <- caTools::colAUC(
          X = r_wilcox_adjp,
          y = (seq_along(mysimu$dea$logFC_theoretical) %in% diffg))
        ## save result
        auci[, j] <- c(r_mssc20$auc_vi, r_mssc20$auc_opt, r_mssc20$auc_glm,
          auc_pseudo_deseq2_no_inds,
          auc_t, auc_wilcox)
      },
      error = function(cond) {
        message(cond)
      },
      finally = {
        message(str_glue("In {i}th repeat under {ncell} cells per individual",
          "when logfc is {logfc_threshold}:", .sep = "\n"))
        print(auci)
      }) ## end of tryCatch for de analysis
    } ## end of ncells
    aucs[, , i] <- auci
    message(str_glue("In {i}th repeat,", "when logfc is {logfc_threshold}:", .sep = "\n"))
    print(aucs)
    ## save result even per turn
    ## in case some error happens
    symsim_prefix <- str_glue(
      "{ngene}gene", "{nind_all}ind",
      "{ncell}cell", "{brn_len}w", "{bimod}bimod",
      "{sigma}sigma", "{capt_alpha}alpha", .sep = "_")
    saveRDS(object = aucs,
      file = file.path(dea_dir, str_glue("dea_auc_{symsim_prefix}.rds")))
  } ## end of rpt
} ## end of main function

option_list <- list(
  ## make_option(c("--ncell_per_ind"), action = "store", type = "integer", default = 30),
  make_option(c("--ngene"), action = "store", type = "integer", default = 100),
  make_option(c("--nind_per_cond"), action = "store", type = "integer", default = 3),
  make_option(c("--brn_len"), action = "store", type = "double", default = 0.5),
  make_option(c("--bimod"), action = "store", type = "double", default = 1),
  make_option(c("--sigma"), action = "store", type = "double", default = 0.6),
  make_option(c("--capt_alpha"), action = "store", type = "double", default = 0.2)
)
args <- parse_args(OptionParser(option_list = option_list))

main(
  nind_per_cond = args$nind_per_cond,
  brn_len = args$brn_len,
  bimod = args$bimod,
  sigma = args$sigma,
  ncells = c(30, 50, 80, 120, 160, 240, 300),
  ngene = args$ngene,
  capt_alpha = args$capt_alpha,
  rpt = 3,
  save_mssc_model = FALSE,
  logfc_threshold = 0.8)

## * test
## test <- function() {
##   main(nind_per_cond = 5, brn_len = 0.5, bimod = 1, sigma = 0.6,
##     ncells = c(30, 50), capt_alpha = 0.2,
##     ngene = 200, rpt = 1, save_mssc_model = FALSE, logfc_threshold = 0.8)
## }

## test()
