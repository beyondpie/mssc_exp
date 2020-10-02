## check by raising #cells/patient, if the top-ranked degenes
## are bias towards the sparse genes.
## Sparse genes are defined as having much zeros in cells.
## Bias means that pseudobulk may tend to select the
## sparse genes if we add more cells.


library(Seurat)
suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)

## local modules
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")

## warnings/errors traceback settings
options(error = traceback)
options(warn = 0)

## * configs
cell_type  <- "Naive CD4+ T"

## * functions
get_celltype_specific_scdata <- function(cntmat, mresponses,
                                         mpatients,
                                         mcelltypes, whichcelltypes) {
  select_cols <- which(mcelltypes == whichcelltypes)
  invisible(list(cnt = cntmat[, select_cols],
                 resp = mresponses[select_cols],
                 inds = mpatients[select_cols],
                 ct = mcelltypes[select_cols]))
}

show_sc_stat <- function(inds, resp){
  message("The individuals stats")
  print(table(as.factor(inds)))
  message("The responses stats")
  print(table(as.factor(resp)))
}

## * load data
pbmc_seurat <- readRDS(here("data", "antiPDL1_PBMC_IL8",
                            "GSE145281_nicol.RDS"))

cntmat <- as.matrix(pbmc_seurat@assays$RNA@counts)
mresponses <- pbmc_seurat@meta.data$response
mpatients <- pbmc_seurat@meta.data$patient
mcelltypes <- pbmc_seurat@meta.data$cellType

uniqcelltypes <- levels(factor(mcelltypes))
cellstats <- data.frame(table(factor(mcelltypes)))
colnames(cellstats) <- c("celltype", "freq")

## * handling a specific cell type
myscdata <- get_celltype_specific_scdata(cntmat, mresponses,
                                         mpatients, mcelltypes,
                                         cell_type)
show_sc_stat(myscdata$inds, myscdata$resp)















