## simulate data using symsim

## control simulation data quality
## generate a series of cells

## * set R environment
suppressPackageStartupMessages(library(SymSim))
suppressPackageStartupMessages(library(tidyverse))
## library(MCMCpack)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(grid)
library(gtable)
library(gridExtra)
suppressPackageStartupMessages(library(ggpubr))

## warnings/errors traceback settings
options(error = traceback)
options(warn = -1)
options(mc.cores = 3)

## load model
mssc_path <- here::here("src", "mssc")
mssc_20 <- modules::import(file.path(mssc_path, "mssc_2-0", "mssc.R"))
mssc_21 <- modules::import(file.path(mssc_path, "mssc_2-1", "mssc.R"))

## load help functions
options("import.path" = here::here("rutils"))
mypseudo <- modules::import("pseudobulk")
mysymsim <- modules::import("mysymsim")

## * utils
plotviolin <- function(cnt, ind, genes) {
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
                                         pnrow = 5) {
  ## show batch effect using violin plot for different genes.
  ## nde: num of differential genes to show
  ## nnde: num of non-differential genes to show
  diffg <- symsim_umi$diffg
  pviolin_symsimdeg <- plotviolin(symsim_umi$obs, symsim_umi$ind, diffg)
  n <- ifelse(nde > length(diffg), length(diffg), nde)
  sampled_pvd <- ggarrange(
    plotlist = pviolin_symsimdeg[sample(length(diffg), n, replace = F)],
    nrow = pnrow,
    ncol = ceiling(nde / pnrow)
  )

  nondiffg <- symsim_umi$nondiffg
  pviolin_symsim_sndeg <- plotviolin(symsim_umi$obs, symsim_umi$ind, nondiffg)
  m <- ifelse(nnde > length(nondiffg), length(nondiffg), nnde)
  sampled_pvsnd <- ggarrange(
    plotlist = pviolin_symsim_sndeg[sample(length(nondiffg), m, replace = F)],
    nrow = pnrow,
    ncol = ceiling(nnde / pnrow)
  )
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
  ind[cond == 2] <- rep((nind + 1): num_total_ind, each = ncell)
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
    cnt = cnt, s = s, cond = cond, ind = ind)

  ## init_mucond: ngene by ncond
  ## only two conditions
  init_mucond <- init_params$mgsnb[, 3:4]
  ngene <- nrow(init_mucond)
  est_mucond_sizes <- vapply(seq_len(ngene), function(i) {
    ## two conditions
    invisible(max(abs(init_mucond[i, 1] - init_mucond[i, 2])))
  })
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
  invisible(lapply(seq_along(g2indeff$dg), function(i) {
    g <- g2indeff$dg[i]
    t <- y2c[g, ]
    if (group_shift) {
      t[t == 0] <- 1
    }
    r[g, ] <- t * exp(g2indeff$dgeff[i, ind])
  }))

  ## add indeff for non-differential genes
  invisible(lapply(seq_along(g2indeff$nondg), function(i) {
    g <- g2indeff$nondg[i]
    t <- y2c[g, ]
    if (group_shift) {
      t[t == 0] <- 1
    }
    r[g, ] <- t * exp(g2indeff$nondgeff[i, ind])
  }))
  return(invisible(round(r)))
}

simu_symsim_with_indeffect <- function(myseed = 1,
                                       save_data_path,
                                       vary = "s",
                                       ncell = 100,
                                       nind = 5,
                                       ngene = 50,
                                       nevf = 10,
                                       n_de_evf = 6,
                                       sigma = 0.2,
                                       ratio_ind2cond = 1.0,
                                       nindeff = 2,
                                       group_shift = TRUE,
                                       scale_in_diffg = 1.0,
                                       scale_in_nondiffg = 1.0,
                                       mssc_model,
                                       save_data = TRUE) {
  ## ncell: num of cell per individual
  ## nind: num of ind per condition, here we only consider two conditions.
  ## simulate the true
  symsim_true <- mysymsim$sim_symsim_true(
    myseed = myseed,
    ncell = ncell * nind * 2,
    ngene = ngene,
    hasgenemodule = F, minmodn = 50,
    npop = 2, nevf = nevf, n_de_evf = n_de_evf,
    sigma = sigma, vary = vary
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
    symsim_umi$counts, ind,
    g2indeff,
    add_on_diffg = group_shift
  )

  symsim_umi$ind <- ind
  symsim_umi$cond <- cond
  symsim_umi$voi <- voi
  symsim_umi$diffg <- symsim_degenes
  symsim_umi$nondiffg <- symsim_ndegs
  symsim_umi$dea <- symsim_dea

  if (save_data) {
    saveRDS(
      object = symsim_umi,
      file = file.path(
        save_data_path,
        stringr::str_glue("{ngene}gene_{nind*2}ind_{ncell}cell_{myseed}.rds")
      )
    )}
  return(invisible(symsim_umi))
}

init_params_and_data <- function(symsim_umi) {
  y2c <- round(symsim_umi$obs)
  ind <- symsim_umi$ind
  cond <- symsim_umi$cond
  sumcnt <- colSums(y2c)
  s <- sumcnt / median(sumcnt)

  hi_params <- high$set_hi_params(
    k = max(ind),
    j = max(cond),
    g = nrow(y2c),
    cnt = y2c,
    s = s,
    cond = cond, ind = ind,
    scale = 1.96^2
  )

  data <- high$to_hbnb_data(y2c, ind, cond, s, hi_params$hp)
  return(invisible(list(hip = hi_params, data = data)))
}

get_auc_hbnb <- function(vifit, data, degs, ndegs) {
  mu_cond <- high$extract_vifit(vifit, data, "mu_cond")
  est_params <- lapply(high$nm_params, function(nm) {
    high$extract_vifit(vifit, data, nm)
  })
  names(est_params) <- high$nm_params

  std_cond <- mean(sqrt(est_params$varofcond))

  rank_stats <- high$get_rank_statistics(mu_cond,
    c1 = 1, c2 = 2,
    std_cond = std_cond
  )
  auc <- high$get_auc(rank_stats, degs, ndegs)
  return(invisible(auc))
}


main <- function(ratio_ind2cond = 0.5,
                 scale_in_diffg = 1.0,
                 scale_in_nondiffg = 1.0,
                 ngene = 50,
                 rpt = 20,
                 adapt_engaged = FALSE,
                 eta = 0.5,
                 adapt_iter = 1000,
                 algorithm  =  "meanfield") {
  ## * configs
  rpt <- rpt
  ## rpt <- 2
  ngene <- ngene
  ## nind for one condition
  nind <- 5
  nindeff <- 2
  ## ncond <- 2
  ncells <- c(160, 200, 300)
  ## ncells <- c(10, 20)
  ## symsim related
  nevf <- 10
  n_de_evf <- 7
  sigma <- 0.2
  vary <- "s"
  ## evf_center <- 1 # should always fix as 1
  ## evf_type <- "discrete"
  ratio_ind2cond <- ratio_ind2cond
  add_on_diffg <- TRUE
  scale_in_diffg <- scale_in_diffg
  scale_in_nondiffg <- scale_in_nondiffg

  ## * configs
  ## myggtitle <- theme(plot.title = element_text(size = 15, hjust = 0.5))
  save_data_path <- here::here("src", "modelcheck", "high2", "symsim")

  ## * calc auc

  nms_auc <- c(
    "t", "bf", "m", "t_rsis", "bf_rsis", "m_rsis", "pseudo")

  result <- array(
    data = NA, dim = c(
      length(nms_auc),
      length(ncells),
      rpt
    ),
    dimnames = list(nms_auc, ncells, NULL)
  )
  for (i in seq_len(rpt)) {
    rpt_slice <- matrix(NA,
      nrow = length(nms_auc),
      ncol = length(ncells)
    )
    for (j in seq_len(length(ncells))) {
      ncell <- ncells[j]
      mssc <- high$High2$new(
        stan_snb_path = here::here(
          "src", "modelcheck", "high2",
          "stan", "snb.stan"
        ),
        stan_snb_cond_path = here::here(
          "src", "modelcheck", "high2",
          "stan", "snb_cond.stan"
        ),
        stan_high2_path = here::here(
          "src", "modelcheck", "high2",
          "stan", "high2.stan"
        ),
        nind = nind * 2,
        tol_rel_obj = 1e-06,
        adapt_iter = adapt_iter,
        adapt_engaged = ifelse(adapt_engaged > 0, TRUE, FALSE),
        algorithm = algorithm,
        eta = eta
      )

      symsim_umi <- simu_symsim_with_indeffect(
        myseed = i,
        save_data_path = save_data_path,
        vary = vary,
        ncell = ncell,
        nind = nind,
        ngene = ngene,
        nevf = nevf,
        n_de_evf = n_de_evf,
        sigma = sigma,
        ratio_ind2cond = ratio_ind2cond,
        nindeff = nindeff,
        add_on_diffg = add_on_diffg,
        scale_in_diffg = scale_in_diffg,
        scale_in_nondiffg = scale_in_nondiffg,
        mssc_model = mssc
      )
      ## ** get symsim violin plots
      diffg <- symsim_umi$diffg
      nondiffg <- symsim_umi$nondiffg
      message(stringr::str_glue("ncell: {ncell}"))
      message(stringr::str_glue("diffg: {length(diffg)}"))
      message(stringr::str_glue("nondiffg: {length(nondiffg)}"))
      ## vln <- plot_genes_after_batcheffect(symsim_umi,
      ##   nde = length(diffg),
      ##   nnde = length(nondiffg)
      ## )
      ## vln$spvln_degs
      ## vln$spvln_ndegs
      pseudo_deseq2_res <- mypseudo$pseudobulk_deseq2(
        symsim_umi$obs,
        symsim_umi$ind,
        factor(symsim_umi$cond)
      )
      pseudo_auc <- mypseudo$calc_auc(
        pseudo_deseq2_res, diffg,
        nondiffg
      )

      y2c <- round(symsim_umi$obs)
      ind <- symsim_umi$ind
      cond <- symsim_umi$cond
      sumcnt <- colSums(y2c)
      s <- sumcnt / median(sumcnt)

      init_all_params <- mssc$init_params(
        cnt = y2c,
        s = s,
        cond = cond,
        ind = ind)

      data <- mssc$to_model_data(cnt = y2c,
                                 s = s,
                                 cond = cond,
                                 ind = ind,
                                 hp = init_all_params$hp)
      mssc$run(data = data, list_wrap_ip = list(init_all_params$ip))

      mucond <- mssc$extract_draws(
        param = "mucond", ngene = nrow(y2c),
        genenms = NULL)
      rankings <- mssc$get_ranking_statistics(
        mucond = mucond,
        two_hot_vec = c(1, -1)
      )

      aucs <- high$get_auc(rankings, c1 = diffg, c2 = nondiffg)

      psis <- mssc$psis()
      print(psis$psis$diagnostics$pareto_k)
      rsis_rankings <- mssc$get_rsis_ranking_statistics(
        ngene = nrow(y2c),
        two_hot_vec = c(1,-1),
        genenms = NULL,
        normweights = psis$normweights
      )

      aucs_rsis <- high$get_auc(rsis_rankings, c1 = diffg, c2 = nondiffg)
      rpt_slice[, j] <- c(aucs, aucs_rsis, pseudo_auc$auc)

      ## show result
      print(rpt_slice)
      ## save symsim-related object
      ## saveRDS(
      ##   object = list(
      ##     symsim = symsim_umi,
      ##     pd = pd,
      ##     vifit = vifit_symsim,
      ##     est_params = est_params,
      ##     diffg = diffg,
      ##     nondiffg = nondiffg
      ##   ),
      ##   file = file.path(
      ##     save_data_path,
      ##     stringr::str_glue("{ngene}gene_{nind*2}ind_{ncell}cell_{ratio_ind2cond}_{scale_in_diffg}_{scale_in_nondiffg}_{i}.rds")
      ##   )
      ## )
    } ## end of ncells
    result[, , i] <- rpt_slice
    print(result)
    ## saveRDS(
    ##   object = result,
    ##   file = file.path(save_data_path, stringr::str_glue("{ngene}gene_{nind*2}ind_{ratio_ind2cond}_{scale_in_diffg}_{scale_in_nondiffg}.rds"))
    ## )
  }
}

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
    default = 0.2
  ),
  make_option(c("--scale_in_nondiffg"),
    action = "store",
    type = "double",
    default = 1.0
  ),
  make_option(c("--ngene"),
    action = "store",
    type = "integer",
    default = 50
  ),
  make_option(c("--rpt"),
    action = "store",
    type = "integer",
    default = 3
    ),
  ## for high method's parameter
  make_option(c("--algorithm"),
              action = "store",
              type = "character",
              default = "meanfield"),
  make_option(c("--eta"),
              action = "store",
              type = "double",
              default = 0.1),
  make_option(c("--adapt_engaged"),
              action = "store",
              type = "integer",
              default = 1),
  make_option(c("--adapt_iter"),
              action = "store",
              type = "integer",
              default = 500)
)

args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()

adapt_engaged <- ifelse(test = (args$adapt_engaged > 0),
                        yes = TRUE,
                        no = FALSE)
main(
  ratio_ind2cond = args$ratio_ind2cond,
  scale_in_diffg = args$scale_in_diffg,
  scale_in_nondiffg = args$scale_in_nondiffg,
  ngene = args$ngene,
  rpt = args$rpt,
  adapt_engaged = adapt_engaged,
  adapt_iter = args$adapt_iter,
  eta = args$eta,
  algorithm = args$algorithm
)
