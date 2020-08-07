suppressPackageStartupMessages(library(SymSim))
suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)
import::from(stringr, str_glue)
library(ggpubr)
options("import.path" = here("rutils"))
myt <- modules::import("transform")
mysymsim <- modules::import("mysymsim")

## * configs
nbatch <- 10

## * utils
my_sim_umi <- function(seed) {
  symsimtrue <- mysymsim$sim_symsim_true(seed)
  symsimumi <- mysymsim$sim_symsim_obs("UMI", symsimtrue)
  suffix <- str_glue("{ncell}_{ngene}_{npop}_{minmodn}_{seed}.rds")
  saveRDS(
    symsimtrue,
    here("src", "simu", "symsimdata", str_glue("true_{suffix}"))
  )
  saveRDS(
    symsimumi,
    here("src", "simu", "symsimdata", str_glue("umi_{suffix}"))
  )
}

plot_genes_after_batcheffect <- function(symsimumibe, degs, ndegs,
                                         nde = 40, nnde = 40,
                                         nrow = 5) {

    pviolin_symsimdeg <- mysymsim$plotviolin(symsimumibe,degs)
    n <- ifelse(nde > length(degs), length(degs), nde)
    sampled_pvd <- ggarrange(
        plotlist = pviolin_symsimdeg[sample(length(degs), n, replace = F)],
        nrow = nrow,
        nrow = ceiling(nde / nrow))

    pviolin_symsim_sndeg <- mysymsim$plotviolin(symsimumibe, ndegs)
    m <- ifelse(nnde > length(ndegs), length(ndegs), nnde)
    sampled_pvsnd <- ggarrange(
        plotlist = pviolin_symsim_sndeg[sample(length(ndegs), m, replace = F)],
        nrow = 5,
        ncol = ceiling(nnde / 5) )
    invisible(list(pvln_alldegs = pviolin_symsimdeg,
         pvln_allndegs = pviolin_symsim_sndeg,
         spvln_degs = sampled_pvd,
         spvln_ndegs = sampled_pvsnd))
}

## * get repeats of simulations
# lapply(seq_len(10), mysymsim)

## * eval simulation process
## ** add individual effects
symsimtrue <- mysymsim$sim_symsim_true(vary = "s")

## de analysis on symsim_true data
symsim_dea <- mysymsim$symsim_de_analysis(symsimtrue,
  popA_idx = which(symsimtrue$cell_meta$pop == 1),
  popB_idx = which(symsimtrue$cell_meta$pop == 2)
)

symsim_degenes <- mysymsim$get_symsim_degenes(symsim_dea) %>% which(. == T)
symsim_strict_ndegenes <- mysymsim$get_symsim_strict_ndegnes(
                                       symsim_dea, logFC = 0.5) %>% which(. == T)

## generate observed umi data
symsimumi <- mysymsim$sim_symsim_obs("UMI", symsimtrue)
## set batches for the two populations
batchids <- mysymsim$assign_batches_2pop(symsimumi, nbatch)
## add batch effect for all the genes and all the batches
symsimumibe <- mysymsim$add_batch_effect(
                            symsimumi, nbatch = nbatch,
                            onbatches = seq_len(1, nbatch),
                            batchids = batchids,
                            batch_effect_size = rep(1.0, nbatch),
                            batch_factor_sd = 0.01,
                            gene_mean_sd = 1,
                            isg2brandom = T
                        )
## *** quick view the batch effect on genes
pvln_be <- plot_genes_after_batcheffect(symsimumibe, symsim_degenes,
                                        symsim_strict_ndegenes)


## ** create false positive genes
## hbe short for high batch effect; lbe short for low batch effect
fphbe <- 1.0
fplbe <- -1.0
fpgenes <- sample(symsim_strict_ndegenes, 20, replace = F)
symsimumi2be <- mysymsim$add_batch_effect(
  symsimumibe,
  nbatch = nbatch,
  ongenes = fpgenes,
  onbatches = c(1, 2, 6),
  batchids = batchids,
  batch_effect_size = c(rep(fphbe, 2), rep(fplbe, 3), fphbe, rep(fplbe, 4)),
  batch_factor_sd = 0.01,
  gene_mean_sd = 0.2,
  isg2brandom = F
)

## ** eval de/nde genes after batch effect
pvln_2be <- plot_genes_after_batcheffect(symsimumi2be, symsim_degenes,
                                         symsim_strict_ndegenes)
