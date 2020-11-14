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

## * configs
rpt <- 5
ngene <- 50
nind <- 5
ncond <- 2
ncells <- c(10, 20, 40, 80, 100, 200)
nfp <- 20
## symsim related
nevf <- 10
n_de_evf <- 7
sigma <- 0.2
vary <- "all"
evf_center <- 1 # should always fix as 1
evf_type <- "discrete"
ratio_ind2cond <- 0.3

## * configs
save_data_path <- here::here("data", "symsim")

## * utils
plot_genes_after_batcheffect <- function(symsimumibe, degs, ndegs,
                                         nde = 40, nnde = 40,
                                         pnrow = 5) {
  pviolin_symsimdeg <- mysymsim$plotviolin(symsimumibe, degs)
  n <- ifelse(nde > length(degs), length(degs), nde)
  sampled_pvd <- ggarrange(
    plotlist = pviolin_symsimdeg[sample(length(degs), n, replace = F)],
    nrow = pnrow,
    ncol = ceiling(nde / pnrow)
  )

  pviolin_symsim_sndeg <- mysymsim$plotviolin(symsimumibe, ndegs)
  m <- ifelse(nnde > length(ndegs), length(ndegs), nnde)
  sampled_pvsnd <- ggarrange(
    plotlist = pviolin_symsim_sndeg[sample(length(ndegs), m, replace = F)],
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
                            nindeff = 2) {
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

  get_indeff <- function(cond) {
    t_result <- vapply(cond, function(i) {
      tmp <- rep(0.0, nind * 2)
      inds <- which(cond_of_ind == i)
      ind_added_eff <- sample(inds, size = nindeff, replace = FALSE)
      tmp[ind_added_eff] <- truncnorm::rtruncnorm(nindeff,
        a = 0.0,
        b = 1.0,
        mean = 0.0,
        sd = variation_of_ind
      )
      return(invisible(tmp))
    }, FUN.VALUE = rep(0.0, nind * 2))

    return(invisible(t(t_result)))
  }
  ## For ndiffg
  indeff_nondiffg <- get_indeff(nondiffg_cond)
  ## For diffg
  indeff_diffg <- get_indeff(low_exp_cond)
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
                                  g2indeff) {
  ## y2c: ngene by ncell
  ## ind: index of individual for each cell
  result <- y2c
  s <- colSums(y2c)
  diffg <- g2indeff$dg
  nondiffg <- g2indeff$nondg
  for (i in seq_len(length(diffg))) {
    g <- diffg[i]
    result[g, ] <- y2c[g, ] * (s * exp(g2indeff$dgeff[i, ind]))
  }

  for (i in seq_len(length(nondiffg))) {
    g <- nondiffg[i]
    result[g, ] <- y2c[g, ] * (s * exp(g2indeff$nondgeff[i, ind]))
  }

  return(invisible(result))
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
                                       ratio_ind2cond = 0.2) {
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
    (symsim_dea$logFC_theoretical[symsim_degenes] < 0)
  nondiffg_cond <- sample(c(1, 2),
    size = length(symsim_ndegs),
    replace = TRUE
  )


  ind2cond <- assign_ind_for_two_cond(symsim_umi$cell_meta$pop,
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
    nindeff = 3
  )
  ## add individual effect to the observed counts directly
  symsim_umi$obs <- add_individual_effect(
    symsim_umi$counts, ind,
    g2indeff
  )

  saveRDS(
    object = symsim_umi,
    file = file.path(
      save_data_path,
      stringr::str_glue("{ngene}gene_{nind*2}ind_{ncell}cell_{myseed}.rds")
    )
  )

  ## ## *** quick view the batch effect on genes
  ## pvln_be <- plot_genes_after_batcheffect(
  ##   symsimumibe, symsim_degenes,
  ##   symsim_degenes
  ## )
  ## ggsave(
  ##   filename = str_glue("{plotprefix}_be_deg_{myseed}.pdf"),
  ##   plot = pvln_be$spvln_degs,
  ##   width = 25, height = 10
  ## )
  ## ggsave(
  ##   filename = str_glue("{plotprefix}_be_ndeg_{myseed}.pdf"),
  ##   plot = pvln_be$spvln_ndegs,
  ##   width = 25, height = 10
  ## )
  return(invisible(symsim_umi))
}


## * debug
myseed <- 1
save_data_path
vary <- "s"
ncell <- 100
nind <- 5
ngene <- 50
nevf <- 10
n_de_evf <- 6
sigma <- 0.2
ratio_ind2cond <- 0.2
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
  (symsim_dea$logFC_theoretical[symsim_degenes] < 0)
nondiffg_cond <- sample(c(1, 2),
  size = length(symsim_ndegs),
  replace = TRUE
)


ind2cond <- assign_ind_for_two_cond(symsim_umi$cell_meta$pop,
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
  nindeff = 2
)
## add individual effect to the observed counts directly
symsim_umi$obs <- add_individual_effect(
  symsim_umi$counts, ind,
  g2indeff
)

saveRDS(
  object = symsim_umi,
  file = file.path(
    save_data_path,
    stringr::str_glue("{ngene}gene_{nind*2}ind_{ncell}cell_{myseed}.rds")
  )
)

## ## *** quick view the batch effect on genes
## pvln_be <- plot_genes_after_batcheffect(
##   symsimumibe, symsim_degenes,
##   symsim_degenes
## )
## ggsave(
##   filename = str_glue("{plotprefix}_be_deg_{myseed}.pdf"),
##   plot = pvln_be$spvln_degs,
##   width = 25, height = 10
## )
## ggsave(
##   filename = str_glue("{plotprefix}_be_ndeg_{myseed}.pdf"),
##   plot = pvln_be$spvln_ndegs,
##   width = 25, height = 10
## )
return(invisible(symsim_umi))


## * generate data
## j = 1
## lapply(seq_len(myrep), FUN = function(i) {
##   simu_symsim_with_indeffect(myseed = i,
##                              save_data_path = save_data_path,
##                              vary = vary,
##                              ncell = ncells[j],
##                              nind = nind,
##                              ngene = ngene,
##                              nevf = nevf,
##                              n_de_evf = n_de_evf,
##                              sigma = sigma,
##                              ratio_ind2cond = ratio_ind2cond
##   )
## })
