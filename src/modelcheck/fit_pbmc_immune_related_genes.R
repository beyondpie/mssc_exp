## use lower counts but differentially expressed genes from PBMC dataset.
## see if mssc can rank the genes related with immunne response higher.

## configs
cmdstan_version <- "2.24.1" ## "2.23.0"
my_cell_cluster <- c(2)
pvalue_cutoff <- 0.05

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
pbmc_datadir <- here("data", "antiPDL1_PBMC_IL8", "seurat.RDS")
pbmc_seurat <- readRDS(pbmc_datadir)

## * select genes
cytot_pbmc <- get_partial_pbmc_data(pbmc_seurat, my_cell_cluster)
cytot_pseudobulk_deseq2 <- mypseudo$pseudobulk_deseq2(
                                      cytot_pbmc$cnt,
                                      cytot_pbmc$ind,
                                      cytot_pbmc$cond)

## ** use DESeq2 baseMean to select low-count genes.






















