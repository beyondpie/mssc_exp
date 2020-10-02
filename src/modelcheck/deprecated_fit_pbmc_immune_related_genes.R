## Deprecated since the experiment in logic is not direct,
## and hard to understand. e.g., what happened if the low-count
## genes are themselves not immune-related.


## Aim:
## use lower counts but differentially expressed genes from PBMC dataset.
## see if mssc can rank the genes related with immunne response higher.

## configs
cmdstan_version <- "2.24.1" ## "2.23.0"
my_cell_cluster <- c(2)
pvalue_cutoff <- 0.05
## roughly A is 1.5 times higher expressed than B or
## A is 0.33 times lower expressed than B.
## log2(1.5) ~ 0.58, 1/1.5 ~ 0.67
log2foldchange_cutoff <- log2(1.5)

## load packages
library(cmdstanr)
set_cmdstan_path(path = paste(Sys.getenv("HOME"),
  "softwares",
  paste0("cmdstan-", cmdstan_version),
  sep = "/"))
library(Seurat)
suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)

## load local scripts/modules
options("import.path" = here("rutils"))
mytrform <- modules::import("transform")
mypseudo <- modules::import("pseudobulk")
myfit <- modules::import("myfitdistr")

## for debug
options(error = traceback)
options(warn = 0)

## functions
deconv_pbmc_seurat <- function(pbmc_seurat) {
  return(list(
    cnt = as.matrix(pbmc_seurat@assays$RNA@counts),
    ind = pbmc_seurat@meta.data$patient,
    cellanno = pbmc_seurat@meta.data$seurat_clusters,
    cond = pbmc_seurat@meta.data$response))
}

get_partial_pbmc_data <- function(pbmc_seurat,
                                  my_cell_cluster = c(2),
                                  rmoutliers = F) {
  mypbmc <- deconv_pbmc_seurat(pbmc_seurat)
  selected_cells <- mypbmc$cellanno %in% my_cell_cluster
  cnt <- mypbmc$cnt[, selected_cells]

  if (rmoutliers) {
    ## FIXME: this line will remove all the cells...
    ## outliers <- myfit$grpoutliers(cnt)
    outliers <- rep(F, ncol(cnt))
  } else {
    outliers <- rep(F, ncol(cnt))
  }

  return(list(
    cnt = mypbmc$cnt[, selected_cells][, !outliers],
    ind = paste0("b", mypbmc$ind[selected_cells][!outliers]),
    cond = paste0("c", mypbmc$cond[selected_cells][!outliers]),
    cellanno = selected_cells,
    outliers = outliers
  ))
}

## main
## * load PBMC data
pbmc_seurat <- here("data", "antiPDL1_PBMC_IL8",
                    "GSE145281_nicol.RDS") %>% readRDS(file = .)

## * select genes
cytot_pbmc <- get_partial_pbmc_data(pbmc_seurat, my_cell_cluster)
cytot_pseudobulk_deseq2 <- mypseudo$pseudobulk_deseq2(
  cytot_pbmc$cnt,
  cytot_pbmc$ind,
  cytot_pbmc$cond)

## ** use DESeq2 baseMean to select low-count/expressed genes.
## use p-value instead of p-adj and log2foldchange to select the genes
cytot_pseudobulk_degenes <- intersect(
  which(cytot_pseudobulk_deseq2$pvalue <= pvalue_cutoff),
  which(abs(cytot_pseudobulk_deseq2$log2FoldChange) >=
    log2foldchange_cutoff))

## further filter high-expressed genes
cytot_pseudobulk_deres <- cytot_pseudobulk_deseq2[
  cytot_pseudobulk_degenes, ] %>%
  .[.$baseMean <= quantile(.$baseMean, 0.25), ]
