## simulate data using symsim

## control simulation data quality
## generate a series of cells

## * set R environment
suppressWarnings(suppressMessages({
  library(SymSim)
  library(tidyverse)
  library(cmdstanr)
  library(bayesplot)
  library(posterior)
  library(grid)
  library(gtable)
  library(gridExtra)
  library(ggpubr)
}))

## warnings/errors traceback settings
options(error = traceback)
options(warn = -1)
options(mc.cores = 3)

## load help functions
options("import.path" = list(here::here("rutils"),
                             here::here("src", "mssc")))
mypseudo <- modules::import("pseudobulk")
mysymsim <- modules::import("mysymsim")

## * utils
plotviolin <- function(cnt, ind, genes,
                       myggtitle = theme(
                         plot.title = element_text(
                           size = 15, hjust = 0.5))) {
  n <- length(genes)
  gnm <- paste0("gene", genes)
  select_cnt <- cnt[genes, ]
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
  diffg <- symsim_umi$diffg
  pviolin_symsimdeg <- plotviolin(symsim_umi$obs, symsim_umi$ind, diffg)
  n <- ifelse(nde > length(diffg), length(diffg), nde)
  sampled_pvd <- ggarrange(
    plotlist = pviolin_symsimdeg[sample(length(diffg), n, replace = F)],
    nrow = pnrow,
    ncol = ceiling(nde / pnrow)
  )

  ## violinplot for nondiffg
  nondiffg <- symsim_umi$nondiffg
  pviolin_symsim_sndeg <- plotviolin(symsim_umi$obs, symsim_umi$ind, nondiffg)
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

assign_ind_for_two_cond <- function(cond, nind = 5,
                                    ncell = 100) {
  ## cond:
  ## - only two conditions are considered here
  ## - a vector of 1 and 2 with length: ncell * nind * 2
  ## nind:
  ## - is for each condition
  ## - so total number of individual is nind * 2.
  ## ncell is for each individual

  ## return the vector of ind and the condition for each ind

  num_total_ind <- nind * 2
  ind <- rep(0, ncell * num_total_ind)
  ind[cond == 1] <- rep(seq_len(nind), each = ncell)
  ind[cond == 2] <- rep((nind + 1):num_total_ind, each = ncell)
  cond_of_ind <- c(rep(1, nind), rep(2, nind))
  return(invisible(list(ind = ind, cond_of_ind = cond_of_ind)))
}

est_variation_of_ind <- function(symsim_obs, ind,
                                 dg,
                                 mssc_model,
                                 ratio_ind2cond = 0.2) {
  ## Given the simulated obs data from symsim,
  ## we estimate the individual effect size based on
  ## mssc model.

  ## Here we only use the gene-wise estimate results,
  ## and use the initial values of mucond. So all the variants of
  ## mssc model could be used here since they share the same
  ## gene-wise module.

  y2c <- symsim_obs$counts[dg, ]
  s <- colSums(y2c)
  cond <- symsim_obs$cell_meta$pop
  ## only use degenes to esimtate the variations
  init_params <- mssc_model$gwsnb$fit_mgsnb(
    cnt = y2c, s = s, cond = cond, ind = ind)

  ## init_mucond: ngene by ncond
  ## only two conditions
  init_mucond <- init_params$mgsnb[, 3:4]
  ngene <- nrow(init_mucond)
  est_mucond_sizes <- vapply(seq_len(ngene), function(i) {
    ## two conditions
    invisible(max(abs(init_mucond[i, 1] - init_mucond[i, 2])))
  }, FUN.VALUE = 0.0)
  variation_of_ind <- ratio_ind2cond * quantile(est_mucond_sizes, probs = 0.975)
  return(invisible(list(hip = init_params, voi = variation_of_ind)))
}

simu_ind_effect <- function(diffg, nondiffg,
                            low_exp_cond,
                            nondiffg_cond,
                            nind,
                            cond_of_ind,
                            variation_of_ind = 0.0001,
                            nindeff = 2,
                            scale_in_diffg = 0.5,
                            scale_in_nondiffg = 0.5) {
  ## simulate individual effect per gene per individual

  ## diffg: index of differentially expressed genes
  ## nondiffg: index of non-differentially expressed genes
  ## low_exp_cond: index of cond for diffg,
  ##   which has the lower expression level.
  ## nondiffg_cond: index of cond for nondiffg,
  ##   could be selected randomly.
  ## nind: num of ind for one condition
  ## cond_of_ind: vector, show individual condition information.
  ## variation_of_ind:
  ## - individual effect size, a scalar
  ## nindeff: num of ind, who shows individual effects
  ##   in one condition. set it as [nind/2] - 1

  get_indeff <- function(cond, scale_level) {
    ## add individual effect for the genes in cond
    ## return: matrix, length(cond) * (nind * 2)
    t_result <- vapply(cond, function(i) {
      tmp <- rep(0.0, nind * 2)
      inds <- which(cond_of_ind == i)
      ind_added_eff <- sample(inds, size = nindeff, replace = FALSE)
      tmp[ind_added_eff] <- variation_of_ind * scale_level
      return(invisible(tmp))
    }, FUN.VALUE = rep(0.0, nind * 2))
    return(invisible(t(t_result)))
  }
  ## For nondiffg
  indeff_nondiffg <- get_indeff(nondiffg_cond, scale_in_nondiffg)
  ## For diffg
  indeff_diffg <- get_indeff(low_exp_cond, scale_in_diffg)
  return(invisible(list(
    nondgeff = indeff_nondiffg,
    dgeff = indeff_diffg,
    dg = diffg,
    nondg = nondiffg
  )))
}

add_individual_effect <- function(y2c, ind,
                                  g2indeff,
                                  group_shift = TRUE) {
  ## add individual effect
  ## - group_shift is TRUE
  ##   all the counts in that individual will be added
  ##   even for zero couts
  ## - group_shift is FALSE
  ##   only non-zero counts will be added

  ## ind: vector, ncell by 1, the individual index for each cell

  ## return:
  ## - count matrix: ngene by ncell (modified by individual effects)

  r <- y2c
  ## add indeff for differential genes
  for (i in seq_along(g2indeff$dg)) {
    g <- g2indeff$dg[i]
    t <- y2c[g, ]
    tmp <- which(t == 0)
    if (group_shift & (length(tmp) > 0)) {
      t[tmp] <- 1
    }
    r[g, ] <- t * exp(g2indeff$dgeff[i, ind])
  }

  ## add indeff for non-differential genes
  for (i in seq_along(g2indeff$nondg)) {
    g <- g2indeff$nondg[i]
    t <- y2c[g, ]
    tmp <- which(t == 0)
    if (group_shift & (length(tmp) > 0)) {
      t[tmp] <- 1
    }
    r[g, ] <- t * exp(g2indeff$nondgeff[i, ind])
  }
  return(invisible(round(r)))
}

simu_symsim_with_indeffect <- function(myseed = 1,
                                       save_data_path,
                                       ncell = 100,
                                       nind = 5,
                                       ngene = 50,
                                       ratio_ind2cond = 1.0,
                                       nindeff = 2,
                                       use_group_shift = TRUE,
                                       scale_in_diffg = 1.0,
                                       scale_in_nondiffg = 1.0,
                                       mssc_model,
                                       save_data = TRUE,
                                       addgenemodule = FALSE) {
  ## ncell: num of cell per individual
  ## nind: num of ind per condition, here we only consider two conditions.
  ## simulate the true
  symsim_true <- mysymsim$sim_symsim_true(
    myseed = myseed,
    ncell = ncell * nind * 2,
    ngene = ngene,
    hasgenemodule = addgenemodule,
    npop = 2, nevf = 10, n_de_evf = 7,
    sigma = 0.2, vary = "s"
  )
  ## generate observed umi data
  symsim_umi <- mysymsim$sim_symsim_obs("UMI", symsim_true)
  symsim_umi$true_cnt <- symsim_true$counts

  ## de analysis on symsim_true data
  symsim_dea <- mysymsim$symsim_de_analysis(symsim_true,
    popA_idx = which(symsim_true$cell_meta$pop == 1),
    popB_idx = which(symsim_true$cell_meta$pop == 2)
  )

  symsim_degenes <- mysymsim$get_symsim_degenes(symsim_dea) %>%
    which(. == T)
  symsim_ndegs <- setdiff(seq_len(ngene), symsim_degenes)

  diffg_low_exp_cond <- 1 +
    (symsim_dea$logFC_theoretical[symsim_degenes] > 0)
  nondiffg_cond <- sample(c(1, 2),
    size = length(symsim_ndegs),
    replace = TRUE
  )

  cond <- symsim_umi$cell_meta$pop
  ind2cond <- assign_ind_for_two_cond(cond,
    nind = nind, ncell = ncell
  )
  ind <- ind2cond$ind
  cond_of_ind <- ind2cond$cond_of_ind

  voi <- est_variation_of_ind(symsim_umi,
    ind, symsim_degenes, mssc_model = mssc_model,
    ratio_ind2cond = ratio_ind2cond
  )

  ## set individual effect for each gene and each individual
  g2indeff <- simu_ind_effect(
    diffg = symsim_degenes,
    nondiffg = symsim_ndegs,
    low_exp_cond = diffg_low_exp_cond,
    nondiffg_cond = nondiffg_cond,
    nind = nind,
    cond_of_ind = cond_of_ind,
    variation_of_ind = voi$voi,
    nindeff = nindeff,
    scale_in_diffg = scale_in_diffg,
    scale_in_nondiffg = scale_in_nondiffg
  )

  ## add individual effect to the observed counts directly
  symsim_umi$obs <- add_individual_effect(
    y2c = symsim_umi$counts, ind = ind,
    g2indeff = g2indeff,
    group_shift = use_group_shift
  )

  symsim_umi$ind <- ind
  symsim_umi$cond <- cond
  symsim_umi$voi <- voi
  symsim_umi$diffg <- symsim_degenes
  symsim_umi$nondiffg <- symsim_ndegs
  symsim_umi$dea <- symsim_dea

  if (save_data) {
    saveRDS(object = symsim_umi, file = save_data_path)
  }
  return(invisible(symsim_umi))
}

get_symsim_by_sampling <- function(symsim_umi,
                                   cellnum = 20) {
  ## create symsim by sampling the cells
  ## from another symsim simulation.
  ## If cellnum is larger than the num of cell per individual,
  ## we use all the cells in that individual.
  ## return
  ## - a new symsim_umi for the sampled cells,
  ##   with an attribute for the sampled cell index.

  ind <- symsim_umi$ind
  sampled_cell_index <- unlist(lapply(seq_len(max(ind)), function(i) {
    cell_index <- which(ind == i)
    num <- length(cell_index)
    if (cellnum >= num) {
      invisible(cell_index)
    }
    invisible(sample(x = cell_index, size = cellnum, replace = FALSE))
  }))
  r <- symsim_umi
  r$obs <- symsim_umi$obs[, sampled_cell_index]
  r$ind <- ind[, sampled_cell_index]
  r$cond <- symsim_umi$cond[, sampled_cell_index]
  r$sampled_cell_index <- sampled_cell_index
  return(invisible(r))
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
                             methods = c("mssc_2-0", "mssc_2-1")) {
  ## setup the result array:
  ## - num_measurement by length(ncell) by repeat_num
  measures <- c("t", "bf", "m", "pbf", "pm")
  nms <- c(unlist(lapply(methods, function(m) {
    invisible(paste(m, measures, sep = "_"))
  })), "pseudo")
  invisible(array(data = NA, dim = c(length(nms), length(ncells), rpt),
    dimnames = list(nms, ncells, seq_len(rpt))))
}

run_mssc <- function(model, symsim, save_result = TRUE,
                     save_path = NULL) {
  ## mssc analysis
  ## side effect
  ## - model state will be updated
  ## - save the results into rds file (save_path).
  ## return
  ## - ranking statistics

  ## data from symsim
  y2c <- round(symsim$obs)
  ind <- symsim$ind
  cond <- symsim$cond
  sumcnt <- colSums(y2c)
  s <- sumcnt / median(sumcnt)

  ## model inference
  init_params <- model$init_params(
    cnt = y2c, s = s, cond = cond, ind = ind)
  data <- model$to_model_data(cnt = y2c,
    s = s, cond = cond, ind = ind, hp = init_params$hp)
  model$run(data = data, list_wrap_ip = list(init_params$ip))

  ## get results
  ## Mannually get the samples for all the params
  ## the model, in principle, will only save the results
  ## in a temporary file since it uses cmdstan.
  est_params <- model$extract_draws_all(
    ngene = nrow(y2c), genenms = seq_len(nrow(y2c)))

  ## mucond: nsample by ngene
  mucond <- model$extract_draws(
    param = "mucond", ngene = nrow(y2c), genenms = seq_len(nrow(y2c)))
  ## three rankings in order: t, bf, m
  ## ngene by 3
  raw_rankings <- model$get_ranking_statistics(
    mucond = mucond, two_hot_vec = c(1, -1))
  ## use PSIS (importance sampling) to further correct the bias
  capture.output(psis <- model$psis())
  ## two rankings in order: bf, m
  ## ngene by 2
  psis_rankings <- model$get_psis_ranking_statistics(
    mucond = mucond, two_hot_vec = c(1, -1), normweights = psis$normweights)

  if (save_result) {
    ## TODO: check loading the results
    ## warning at opt:argparse
    saveRDS(object = list(est_params = est_params, raw_rankings = raw_rankings,
      psis_rankings = psis_rankings, model = model),
    file = save_path)
  }

  return(invisible(
    list(raw_rankings = raw_rankings, psis_rankings = psis_rankings)))
}

main <- function(nind = 5,
                 nindeff = 2,
                 use_group_shift = FALSE,
                 ratio_ind2cond = 0.5,
                 scale_in_diffg = 1.0,
                 scale_in_nondiffg = 1.0,
                 ngene = 80,
                 rpt = 3,
                 ncells = c(20, 40, 80, 100, 120, 160),
                 save_figure = TRUE,
                 width = 20,
                 height = 10,
                 save_symsim_data = TRUE,
                 save_mssc_model = TRUE) {
  ## nind: num of individuals in one condition
  ## nindeff: num of individuals with the side effect
  ## - in symsim simulation, we choose one condition, and
  ##   add individual effects on *nindeff* individuals in that
  ##   condition.
  ## use_group_shift or not
  ## - when TRUE, zero counts will also be added the individual effects.
  ## - when FALSE, only non-zero counts will be added the individual effects.
  ##
  ## return:
  ## - the symsim simulation data
  ## - the violin plot for the genes in the simulation
  ## - the auc results for different methods
  ##   - two versions of mssc and the pseudobulk are considered.

  ## load mssc model
  ## for test firstly
  mssc_20 <- load_mssc(nind = nind * 2, mssc_version = "mssc_2-0")
  mssc_21 <- load_mssc(nind = nind * 2, mssc_version = "mssc_2-1")

  ## path to save symsim simulation data
  symsim_save_path <- here::here("src", "symsim")
  ## experiment name for the usages of saving files
  expnm <- ifelse(use_group_shift, "groupshift", "nonzeroshift")

  ## declare the result
  r <- set_result_array(rpt = rpt, ncells = ncells,
    methods = c("mssc_2-0", "mssc_2-1"))

  for (i in seq_len(rpt)) {
    ## init the record for current repeat
    r_i <- matrix(NA, nrow = dim(r)[1], ncol = dim(r)[2])
    rownames(r_i) <- rownames(r)
    for (j in seq_len(length(ncells))) {
      ncell <- ncells[j]
      ## simulate symsim dataset
      symsim_data_fnm <- stringr::str_glue(
        "{ngene}gene", "{nind*2}ind", "{nindeff}indeff", "{ncell}cell",
        expnm, "seed-{i}", "rpt-{rpt}.rds", .sep = "_")
      save_symsim_data_path <- file.path(
        symsim_save_path, "data", symsim_data_fnm)
      symsim_umi <- simu_symsim_with_indeffect(
        ## seed is set to the current repeat number
        myseed = i,
        save_data_path = save_symsim_data_path,
        nind = nind, nindeff = nindeff,
        use_group_shift = use_group_shift,
        ncell = ncell, ngene = ngene,
        ratio_ind2cond = ratio_ind2cond,
        scale_in_diffg = scale_in_diffg,
        scale_in_nondiffg = scale_in_nondiffg,
        mssc_model = mssc_21,
        save_data = save_symsim_data)

      ## get diff and non-diff genes
      diffg <- symsim_umi$diffg
      nondiffg <- symsim_umi$nondiffg
      message(stringr::str_glue("ncell: {ncell}",
        "diffg: {length(diffg)}",
        "nondiffg: {length(nondiffg)}",
        .sep = "\n"))
      ## draw violin plot
      plotfnm_prefix <- stringr::str_glue(
        "symsim", "{ngene}gene", "{nind*2}ind", "{nindeff}indeff", expnm,
        "{ncell}cell", "rpt-{i}", .sep = "_")
      vln <- plot_genes_after_batcheffect(symsim_umi,
        nde = length(diffg),
        nnde = length(nondiffg),
        save_path = file.path(symsim_save_path, "figures"),
        save_figure = save_figure,
        fnm_prefix = plotfnm_prefix,
        width = width,
        height = height
      )

      ## pseudo + deseq2 analysis
      pseudo_deseq2_res <- mypseudo$pseudobulk_deseq2(
        symsim_umi$obs, symsim_umi$ind, factor(symsim_umi$cond))
      pseudo_auc <- mypseudo$calc_auc(pseudo_deseq2_res, diffg, nondiffg)

      ## mssc model
      msscfnm_prefix <- stringr::str_glue(
        "symsim", "{ngene}gene", "{nind*2}ind", "{nindeff}indeff", expnm,
        "{ncell}cell", "rpt-{i}", .sep = "_")

      ## mssc20
      r_mssc20 <- run_mssc(
        model = mssc_20, symsim = symsim_umi, save_result = save_mssc_model,
        save_path = file.path(
          symsim_save_path, "models",
          stringr::str_glue({msscfnm_prefix}, "_mssc20.rds")))
      raw_auc_mssc20 <- mssc_20$get_auc(r_mssc20$raw_rankings,
        c1 = diffg, c2 = nondiffg)
      psis_auc_mssc20 <- mssc_20$get_auc(r_mssc20$psis_rankings,
        c1 = diffg, c2 = nondiffg)
      ## mssc21
      r_mssc21 <- run_mssc(
        model = mssc_20, symsim = symsim_umi, save_result = save_mssc_model,
        save_path = file.path(
          symsim_save_path, "models",
          stringr::str_glue({msscfnm_prefix}, "_mssc21.rds")))
      raw_auc_mssc21 <- mssc_21$get_auc(r_mssc21$raw_rankings,
        c1 = diffg, c2 = nondiffg)
      psis_auc_mssc21 <- mssc_21$get_auc(r_mssc21$psis_rankings,
        c1 = diffg, c2 = nondiffg)
      ## merge results
      r_i[, j] <- c(raw_auc_mssc20, psis_auc_mssc20,
        raw_auc_mssc21, psis_auc_mssc21,
        pseudo_auc$auc)
      print(r_i)
    } ## end of ncells
    r[, , i] <- r_i
    ## show the matrix for current result
    print(r)
    ## save result per repeat for safe.
    fnm <- stringr::str_glue(
      "{ngene}gene", "{nind*2}ind", "{nindeff}indeff", expnm,
      ratio_ind2cond, scale_in_diffg, "{scale_in_nondiffg}.rds", .sep = "_")
    saveRDS(object = r, file = file.path(symsim_save_path, "results", fnm))
  } ## end of repeat
} ## end of main

library(optparse)
option_list <- list(
  make_option(c("--ratio_ind2cond"),
    action = "store",
    type = "double",
    default = 0.3
  ),
  make_option(c("--scale_in_diffg"),
    action = "store",
    type = "double",
    default = 1.0
  ),
  make_option(c("--scale_in_nondiffg"),
    action = "store",
    type = "double",
    default = 1.0
  ),
  make_option(c("--ngene"),
    action = "store",
    type = "integer",
    default = 80
  ),
  make_option(c("--rpt"),
    action = "store",
    type = "integer",
    default = 3
  ),
  make_option(c("--groupshift"), action = "store",
    type = "integer", default = 1),
  make_option(c("--nind"), action = "store",
    type = "integer", default = 3),
  make_option(c("--nindeff"), action = "store",
    type = "integer", default = 1)
)

args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()

main(
  nind = args$nind,
  nindeff = args$nindeff,
  use_group_shift = ifelse(args$groupshift > 0, T, F),
  ratio_ind2cond = args$ratio_ind2cond,
  scale_in_diffg = args$scale_in_diffg,
  scale_in_nondiffg = args$scale_in_nondiffg,
  ngene = args$ngene,
  rpt = 3,
  ## for test
  ncells = c(20, 40, 80, 120, 160, 240, 300),
  save_figure = TRUE,
  width = 20,
  height = 10,
  save_symsim_data = TRUE,
  save_mssc_model = FALSE
)
