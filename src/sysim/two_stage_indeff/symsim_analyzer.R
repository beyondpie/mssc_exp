options(error = traceback)
options(warn = -1)
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(bayesplot))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(import::from(here, here))
suppressPackageStartupMessages(library(DESeq2))
import::from(optparse, make_option, OptionParser, parse_args)
import::from(stringr, str_glue)

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myroc <- modules::import("roc")
mysymsim <- modules::import("mysymsim")
mypseudo <- modules::import("pseudobulk")

## * configs
myseed <- 1
version <- "v1-1"
method <- "vi"
symsim_data_dir <- here("data", "symsim", "twostage_be_symsim", "data")
symsim_exp_dir <- here("exps", "symsim", "stan", myseed)
symsim_exp_vi_dir <- paste(symsim_exp_dir, "vi")
symsim_exp_mc_dir <- paste(symsim_exp_dir, "mc")

## * functions
myviolin <- function(genes, mymsscs, mypseudos, prow = 5,
                     symsimcnt = symsim2be,
                     myversion = version,
                     mymethod = method) {
  msscs <- mymsscs[genes]
  pseudos <- mypseudos[genes]
  violin_genes <- mysymsim$plotviolin(symsimcnt, genes)
  violin_genes_with_scores <- seq_len(length(genes)) %>% map(.f = function(i) {
    t <- str_glue(
      "{myversion}({mymethod}): {myt$fmtflt(msscs[i])}; ",
      "{myt$fmtflt(pseudos[i])}"
    )
    p <- violin_genes[[i]] + labs(title = t)
    invisible(p)
  })

  invisible(ggarrange(
    plotlist = violin_genes_with_scores,
    nrow = prow,
    ncol = ceiling(length(ndegs) / prow)
  ))
}

plot_genes <- function(symsimumibe, degs, ndegs,
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


## * get truth
symsimtrue <- readRDS(
  file =
    paste(symsim_data_dir,
      str_glue("symsim_true_{myseed}.rds"),
      sep = "/"
    )
)
symsim_dea <- mysymsim$symsim_de_analysis(symsimtrue,
  popA_idx = which(symsimtrue$cell_meta$pop == 1),
  popB_idx = which(symsimtrue$cell_meta$pop == 2)
)

numofgene <- nrow(symsimtrue$counts)

symsim_degenes <- mysymsim$get_symsim_degenes(symsim_dea,
  nDiffEVF = 1,
  logFC = 0.6
) %>% which(. == T)

symsim_ndegs <- setdiff(seq_len(numofgene), symsim_degenes)
symsim_strict_ndegenes <- mysymsim$get_symsim_strict_ndegenes(symsim_dea,
  nDiffEVF = 0,
  logFC = 0.5
) %>% which(. == T)
symsim_sampled_ndegs_1 <- sample(symsim_ndegs,
  size = length(symsim_degenes),
  replace = F
)
symsim_sampled_ndegs_2 <- sample(symsim_ndegs,
  size = length(symsim_degenes),
  replace = F
)
symsim_sampled_ndegs_3 <- sample(symsim_ndegs,
  size = length(symsim_degenes),
  replace = F
)

## * load simulated counts
symsimumi <- readRDS(
  file = paste(symsim_data_dir,
    str_glue("symsim_umi_{myseed}.rds"),
    sep = "/"
  )
)
symsimbe <- readRDS(
  file = paste(symsim_data_dir,
    str_glue("symsim_be_{myseed}.rds"),
    sep = "/"
  )
)
symsim2be <- readRDS(
  file = paste(symsim_data_dir,
    str_glue("symsim_2be_{myseed}.rds"),
    sep = "/"
  )
)

## * supplement: genes violins when no batch effect
pvln_umi <- plot_genes(
  symsimumi, symsim_degenes,
  symsim_strict_ndegenes
)

## * analyzer the results (compared with pseudobulk)
## ** mssc scores
degs <- symsim_degenes
ndegs <- symsim_strict_ndegenes
## * supplement: genes violins when no batch effect
symsimumi$batch_meta <- list(batch = symsimbe$batch_meta$batch)

pvln_umi <- plot_genes_after_batcheffect(
  symsimumi, symsim_degenes,
  symsim_strict_ndegenes
)
resultdir <- here("data", "symsim", "twostage_be_symsim")
dprefix <- paste(resultdir, "data", "symsim", sep = "/")
plotprefix <- paste(resultdir, "pdf", "sample", sep = "/")

ggsave(
  filename = str_glue("{plotprefix}_umi_deg_{myseed}.pdf"),
  plot = pvln_umi$spvln_degs,
  width = 25, height = 10
)
ggsave(
  filename = str_glue("{plotprefix}_umi_ndeg_{myseed}.pdf"),
  plot = pvln_umi$spvln_ndegs,
  width = 25, height = 10
)


stanfit <- myt$load_stan(
  dirnm = symsim_exp_dir,
  modelnm = version, method = method,
  vi_dir = "vi", mc_dir = "mc"
)
dcond_mucond <- myt$get_ctrlmnscase_par(stanfit, "MuCond")
dt_mucond <- myt$calt(dcond_mucond, fn = colMeans)

## ** pseudo bulk analysis
mycnt <- symsim2be$counts
mybatches <- symsim2be$batch_meta$batch
myconds <- symsim2be$cell_meta$pop
pseudo_deseq2_res <- mypseudo$pseudobulk_deseq2(
  mycnt,
  mybatches, myconds
) %>% .[["pvalue"]]

## ** violin plot with scores
plotprefix <- here("exps", "symsim", "stan", myseed, method, "p")
## *** on negative samples
pviolin_ndegs <- myviolin(ndegs, dt_mucond, pseudo_deseq2_res,
  myversion = version, prow = 6
)
ggsave(
  filename = str_glue("{plotprefix}_{version}_ndeg.pdf"),
  plot = pviolin_ndegs,
  width = 25, height = 10
)
## *** on positive samples
pviolin_degs <- myviolin(sample(degs, 36), dt_mucond, pseudo_deseq2_res,
  myversion = version, prow = 6
)
ggsave(
  filename = str_glue("{plotprefix}_{version}_deg.pdf"),
  plot = pviolin_degs,
  width = 25, height = 10
)

## *** hist on values of mssc
