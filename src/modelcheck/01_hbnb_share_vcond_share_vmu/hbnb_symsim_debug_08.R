## simulate data using symsim

## control simulation data quality
## generate a series of cells

## * set R environment
source("hbnb_set_r_lib_env_01.R")
suppressPackageStartupMessages(library(SymSim))
suppressPackageStartupMessages(library(ggpubr))
library(truncnorm)

hbnbm <- modules::import("hbnb_mssc_03")

options("import.path" = here::here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")
mysymsim <- modules::import("mysymsim")

options(warn = -1)
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
  diffg <- symsim_umi$diffg
  pviolin_symsimdeg <- plotviolin(symsim_umi$obs, symsim_umi$ind, diffg)
  n <- ifelse(nde > length(diffg), length(diffg), nde)
  sampled_pvd <- ggarrange(
    plotlist = pviolin_symsimdeg[sample(length(diffg), n, replace = F)],
    nrow = pnrow,
    ncol = ceiling(nde / pnrow)
  )
  nondiffg <- symsim_umi$nondiffg
  pviolin_symsim_sndeg <- plotviolin(
    symsim_umi$obs,
    symsim_umi$ind, nondiffg
  )
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
  ## nind is for each ond
  ## ncell is for each individual
  ## return the vector of ind and the condition for each ind
  ind <- rep(0, ncell * nind * 2)
  ind[cond == 1] <- rep(seq_len(nind), each = ncell)
  ind[cond == 2] <- rep((nind + 1):(nind * 2), each = ncell)
  cond_of_ind <- c(rep(1, nind), rep(2, nind))
  return(invisible(list(ind = ind, cond_of_ind = cond_of_ind)))
}

est_variation_of_ind <- function(symsim_obs, ind,
                                 symsim_degenes,
                                 ratio_ind2cond = 0.2) {
  ## given the simulated obs data from symsim,
  ## we estimate the variation of ind for hbnb model.
  y2c <- symsim_obs$counts[symsim_degenes, ]
  k <- max(ind)
  j <- 2
  g <- nrow(y2c)
  s <- colSums(y2c)
  cond <- symsim_obs$cell_meta$pop
  ## only use degenes to esimtate the variations
  hi_params <- hbnbm$set_hi_params(
    k = k, j = j, g = g, cnt = y2c,
    s = s, cond = cond, ind = ind,
    scale = 1.0, to_ind = FALSE
  )
  varofcond <- hi_params$ip$varofcond
  variation_of_ind <- ratio_ind2cond * sqrt(varofcond)
  return(invisible(list(hip = hi_params, voi = variation_of_ind)))
}

simu_ind_effect <- function(diffg, nondiffg,
                            low_exp_cond,
                            nondiffg_cond,
                            nind,
                            cond_of_ind,
                            variation_of_ind = 0.0001,
                            nindeff = 2,
                            scale_in_diffg = 0.1,
                            scale_in_nondiffg = 1.0) {
  ## simulate individual effect per gene per individual
  ## diffg: index of differentially expressed genes
  ## nondiffg: index of non-differentially expressed genes
  ## low_exp_cond: index of cond for diffg,
  ##   which has the lower expression level.
  ## nondiffg_cond: index of cond for nondiffg,
  ##   could be selected randomly.
  ## nind: num of ind for one condition
  ## nindeff: num of ind, who shows individual effects
  ##   in one condition. set it as [nind/2] - 1
  ## cond_of_ind: vector, show individual condition information.

  ## return two matrix:
  ##   diffg_indeff: diffg by ind
  ##   nondiffg_indeff: nondiffg by ind

  get_indeff <- function(cond, scale_level) {
    t_result <- vapply(cond, function(i) {
      tmp <- rep(0.0, nind * 2)
      inds <- which(cond_of_ind == i)
      ind_added_eff <- sample(inds, size = nindeff, replace = FALSE)
      ## tmp[ind_added_eff] <- truncnorm::rtruncnorm(nindeff,
      ##   a = 0.0,
      ##   b = 1.0,
      ##   mean = 0.0,
      ##   sd = variation_of_ind
      ##   )
      tmp[ind_added_eff] <- variation_of_ind * scale_level
      return(invisible(tmp))
    }, FUN.VALUE = rep(0.0, nind * 2))

    return(invisible(t(t_result)))
  }
  ## For ndiffg
  indeff_nondiffg <- get_indeff(nondiffg_cond, scale_in_nondiffg)
  ## For diffg
  indeff_diffg <- get_indeff(low_exp_cond, scale_in_diffg)
  ## g2i <- outer(on_gene, on_ind, FUN = "+")
  ## g2i <- g2i * sample(c(-1, 1), ngene * nind, replace = T)
  return(invisible(list(
    nondgeff = indeff_nondiffg,
    dgeff = indeff_diffg,
    dg = diffg,
    nondg = nondiffg
  )))
}

add_individual_effect <- function(y2c, ind,
                                  g2indeff,
                                  add_on_diffg = TRUE) {
  ## y2c: ngene by ncell
  ## ind: index of individual for each cell
  result <- y2c
  diffg <- g2indeff$dg
  nondiffg <- g2indeff$nondg
  if (add_on_diffg) {
    for (i in seq_len(length(diffg))) {
      g <- diffg[i]
      t <- y2c[g, ]
      ## t[t == 0] <- 1
      result[g, ] <- t * exp(g2indeff$dgeff[i, ind])
    }
  }

  for (i in seq_len(length(nondiffg))) {
    g <- nondiffg[i]
    t <- y2c[g, ]
    ## so that when the observed count is zero,
    ## we will increase the counts.
    ## t[t == 0] <- 1
    result[g, ] <- t * exp(g2indeff$nondgeff[i, ind])
  }

  return(invisible(round(result)))
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
                                       add_on_diffg = TRUE,
                                       scale_in_diffg = 1.0,
                                       scale_in_nondiffg = 1.0) {
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
    ind, symsim_degenes,
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
    add_on_diffg = add_on_diffg
  )
  symsim_umi$ind <- ind
  symsim_umi$cond <- cond
  symsim_umi$voi <- voi
  symsim_umi$diffg <- symsim_degenes
  symsim_umi$nondiffg <- symsim_ndegs
  symsim_umi$dea <- symsim_dea
  saveRDS(
    object = symsim_umi,
    file = file.path(
      save_data_path,
      stringr::str_glue("{ngene}gene_{nind*2}ind_{ncell}cell_{myseed}.rds")
    )
  )
  return(invisible(symsim_umi))
}

init_params_and_data <- function(symsim_umi) {
  y2c <- round(symsim_umi$obs)
  ind <- symsim_umi$ind
  cond <- symsim_umi$cond
  sumcnt <- colSums(y2c)
  s <- sumcnt / median(sumcnt)

  hi_params <- hbnbm$set_hi_params(
    k = max(ind),
    j = max(cond),
    g = nrow(y2c),
    cnt = y2c,
    s = s,
    cond = cond, ind = ind,
    scale = 1.96^2
  )

  data <- hbnbm$to_hbnb_data(y2c, ind, cond, s, hi_params$hp)
  return(invisible(list(hip = hi_params, data = data)))
}

get_auc_hbnb <- function(vifit, data, degs, ndegs) {
  mu_cond <- hbnbm$extract_vifit(vifit, data, "mu_cond")
  est_params <- lapply(hbnbm$nm_params, function(nm) {
    hbnbm$extract_vifit(vifit, data, nm)
  })
  names(est_params) <- hbnbm$nm_params

  std_cond <- mean(sqrt(est_params$varofcond))

  rank_stats <- hbnbm$get_rank_statistics(mu_cond,
    c1 = 1, c2 = 2,
    std_cond = std_cond
  )
  auc <- hbnbm$get_auc(rank_stats, degs, ndegs)
  return(invisible(auc))
}


main <- function() {
  ## * configs
  ## rpt <- 5
  rpt <- 2
  ## ngene <- 80
  ngene <- 20
  nind <- 5
  ## ncond <- 2
  ## ncells <- c(10, 20, 40, 80, 100, 200)
  ncells <- c(10, 20)
  ## symsim related
  nevf <- 10
  n_de_evf <- 8
  sigma <- 0.2
  vary <- "s"
  ## evf_center <- 1 # should always fix as 1
  ## evf_type <- "discrete"
  ratio_ind2cond <- 0.3
  add_on_diffg <- TRUE
  scale_in_diffg <- 0.01
  scale_in_nondiffg <- 1.0
  nindeff <- 2
  add_on_diffg <- TRUE
  scale_in_diffg <- 1.0
  scale_in_nondiffg <- 1.0


  ## * configs
  ## myggtitle <- theme(plot.title = element_text(size = 15, hjust = 0.5))
  save_data_path <- here::here("data", "symsim")

  ## * calc auc

  nms_auc <- c(
    "auc_z", "auc_m", "auc_p", "auc_p10",
    "auc_p05", "auc_p025", "auc_bf", "auc_pseudo"
  )


  result <- array(
    data = NA, dim = c(
      length(nms_auc),
      length(ncells),
      rpt
    ),
    dimnames = list(nms_auc, ncells, NULL)
  )
  for (i in seq_len(rpt)) {
    for (j in seq_len(length(ncells))) {
      ncell <- ncells[j]
      rpt_slice <- matrix(NA,
        nrow = length(nms_auc),
        ncol = length(ncells)
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
        scale_in_nondiffg = scale_in_nondiffg
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

      pd <- init_params_and_data(symsim_umi)
      vifit_symsim <- hbnbm$run_hbnb_vi(data = pd$data, ip = pd$hip$ip)

      ## examine the fitted parameters
      est_params <- lapply(hbnbm$nm_params, function(nm) {
        hbnbm$extract_vifit(vifit_symsim, pd$data, nm)
      })
      names(est_params) <- hbnbm$nm_params

      hbnb_auc <- get_auc_hbnb(vifit_symsim,
        data = pd$data,
        diffg, nondiffg
      )

      rpt_slice[, j] <- c(unlist(hbnb_auc), pseudo_auc$auc)
      ## show result
      print(rpt_slice)
      ## save symsim-related object
      saveRDS(
        object = list(
          symsim = symsim_umi,
          pd = pd,
          vifit = vifit_symsim,
          est_params = est_params,
          diffg = diffg,
          nondiffg = nondiffg
        ),
        file = file.path(
          save_data_path,
          stringr::str_glue("{ngene}gene_{nind*2}ind_{ncell}cell_{ratio_ind2cond}_{scale_in_diffg}_{scale_in_nondiffg}_{i}.rds")
        )
      )
    } ## end of ncells
    result[, , i] <- rpt_slice
    saveRDS(
      object = result,
      file = file.path(save_data_path, stringr::str_glue("{ngene}gene_{nind*2}ind_{ratio_ind2cond}_{scale_in_diffg}_{scale_in_nondiffg}.rds"))
    )
  }
}

main()
