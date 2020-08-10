options(error = traceback)
options(warn = -1)
library(optparse)
suppressPackageStartupMessages(library(SymSim))
suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)
import::from(stringr, str_glue)
suppressPackageStartupMessages(library(ggpubr))

options("import.path" = here("rutils"))
myt <- modules::import("transform")
mysymsim <- modules::import("mysymsim")

## * options from outside
option_list <- list(
  make_option(c("--rep"),
    action = "store",
    type = "integer",
    default = 3
  ),
  make_option(c("--ncell"),
    action = "store",
    type = "integer",
    default = 2000
  ),
  make_option(c("--ngene"),
    action = "store",
    type = "integer",
    default = 300),
  make_option(c("--nbatch"),
    action = "store",
    type = "integer",
    default = 10),
  make_option(c("--nevf"),
    action = "store",
    type = "integer",
    default = 10),
  make_option(c("--n_de_evf"),
    action = "store",
    type = "integer",
    default = 6),
  make_option(c("--sigma"),
    action = "store",
    type = "double",
    default = 0.2),
  make_option(c("--vary"),
    action = "store",
    type = "character",
    default = "s"),
  make_option(c("--fpnum"),
    action = "store",
    type = "integer",
    default = 20)
)

args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()

## * configs
resultdir <- here("data", "symsim", "twostage_be_symsim")
symsim_data_dir <- paste(resultdir, "data", sep = "/")
fprefix <- "symsim_2be"

ncell <- args$ncell
ngene <- args$ngene
nbatch <- args$nbatch
nevf <- args$nevf
n_de_evf <- args$n_de_evf
sigma <- args$sigma
vary <- args$vary
myrep <- args$rep
fpnum <- args$fpnum
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

twostage_batcheffect_symsim <- function(myseed = 1,
                                        nbatch = 10L,
                                        resultdir = resultdir,
                                        vary = "s",
                                        ncell = 2000,
                                        ngene = 300,
                                        nevf = 10,
                                        n_de_evf = 6,
                                        sigma = 0.2, fpnum = 20L
) {
  dprefix <- paste(resultdir, "data", "symsim", sep = "/")
  plotprefix <- paste(resultdir, "pdf", "sample", sep = "/")
  ## simulate the true
  symsimtrue <- mysymsim$sim_symsim_true(myseed = myseed,
    ncell = ncell,
    ngene = ngene,
    hasgenemodule = F, minmodn = 50,
    npop = 2, nevf = nevf, n_de_evf = n_de_evf,
    sigma = sigma, vary = vary
  )
  ## generate observed umi data
  symsimumi <- mysymsim$sim_symsim_obs("UMI", symsimtrue)
  ## set batches for the two populations
  batchids <- mysymsim$assign_batches_2pop(symsimumi, nbatch)
  ## add batch effect for all the genes and all the batches
  symsimumibe <- mysymsim$add_batch_effect(
    symsimumi,
    nbatch = nbatch,
    onbatches = seq_len(1, nbatch),
    batchids = batchids,
    batch_effect_size = rep(1.0, nbatch),
    batch_factor_sd = 0.01,
    gene_mean_sd = 1,
    isg2brandom = T
  )
  ## de analysis on symsim_true data
  symsim_dea <- mysymsim$symsim_de_analysis(symsimtrue,
    popA_idx = which(symsimtrue$cell_meta$pop == 1),
    popB_idx = which(symsimtrue$cell_meta$pop == 2)
  )

  symsim_degenes <- mysymsim$get_symsim_degenes(symsim_dea) %>% which(. == T)
  symsim_strict_ndegenes <- mysymsim$get_symsim_strict_ndegenes(
    symsim_dea,
    logFC = 0.5
  ) %>% which(. == T)

  ## ** create false positive genes
  ## hbe short for high batch effect; lbe short for low batch effect
  fphbe <- 1.0
  fplbe <- -1.0
  nfp <- ifelse(length(symsim_strict_ndegenes) > fpnum, fpnum,
    length(symsim_strict_ndegenes)
  )
  fpgenes <- sample(symsim_strict_ndegenes, nfp, replace = F)
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

  saveRDS(object = symsimtrue, file = str_glue("{dprefix}_true_{myseed}.rds"))
  saveRDS(object = symsimumi, file = str_glue("{dprefix}_umi_{myseed}.rds"))
  saveRDS(object = symsimumibe, file = str_glue("{dprefix}_be_{myseed}.rds"))
  saveRDS(object = symsimumi2be, file = str_glue("{dprefix}_2be_{myseed}.rds"))

  ## *** quick view the batch effect on genes
  pvln_be <- plot_genes_after_batcheffect(
    symsimumibe, symsim_degenes,
    symsim_degenes
  )

  ## ** eval de/nde genes after batch effect
  pvln_2be <- plot_genes_after_batcheffect(
    symsimumi2be, symsim_degenes,
    symsim_strict_ndegenes
  )

  ggsave(
    filename = str_glue("{plotprefix}_be_deg_{myseed}.pdf"),
    plot = pvln_be$spvln_degs,
    width = 25, height = 10
  )
  ggsave(
    filename = str_glue("{plotprefix}_be_ndeg_{myseed}.pdf"),
    plot = pvln_be$spvln_ndegs,
    width = 25, height = 10
  )
  ggsave(
    filename = str_glue("{plotprefix}_2be_deg_{myseed}.pdf"),
    plot = pvln_2be$spvln_degs,
    width = 25, height = 10
  )
  ggsave(
    filename = str_glue("{plotprefix}_2be_ndeg_{myseed}.pdf"),
    plot = pvln_2be$spvln_ndegs,
    width = 25, height = 10
  )
}

## * generate data
message("save args configuration ... ")
saveRDS(file = paste(resultdir, "args.rds", sep = "/"),
  object = args)

lapply(seq_len(myrep), FUN = function(i) {
  message(str_glue("generating the {i}th symsim data ..."))
  twostage_batcheffect_symsim(i,
    nbatch = nbatch,
    resultdir = resultdir,
    vary = vary,
    ncell = ncell,
    ngene = ngene,
    nevf = nevf,
    n_de_evf = n_de_evf,
    sigma = sigma, fpnum = fpnum
  )

  message(str_glue("transforming to the {i}th stan dump data ..."))
  symsimumi2be <- readRDS(paste(symsim_data_dir,
    str_glue("{fprefix}_{i}.rds"),
    sep = "/"
  ))

  ## * transform to rstan
  cnt <- symsimumi2be$counts
  batches <- symsimumi2be$batch_meta$batch
  conds <- symsimumi2be$cell_meta$pop
  totcntpcell <- colSums(cnt)

  myt$to_bagwiff(
    cnt, batches, conds, totcntpcell,
    outf = paste(symsim_data_dir,
      str_glue("{fprefix}_{i}.rdump"),
      sep = "/"),
    rdump = T
  )
})
