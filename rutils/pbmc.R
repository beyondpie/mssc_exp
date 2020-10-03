## functions for pbmc seurat dataset.

## > uniqcelltypes
##  [1] "B cells"                 "CD14- INFG- Monocytes"
##  [3] "CD36+ Monocytes"         "CD68- Monocytes"
##  [5] "CD8A T"                  "CXCR4- NK"
##  [7] "CXCR4+ NK"               "Effector Memory T"
##  [9] "Low-density basophils"   "Monocytes"
## [11] "Myeloid DC"              "Naive CD4+ T"
## [13] "Non classical monocytes" "RUNX3+ NK"
## [15] "Th1"

get_celltype_specific_scdata <- function(cntmat, mresponses,
                                         mpatients,
                                         mcelltypes, whichcelltypes) {
  ## select the subset data based on a specific cell types.
  ## input are the cnt, meta data, and the cell type

  select_cols <- which(mcelltypes == whichcelltypes)
  invisible(list(cnt = cntmat[, select_cols],
    resp = mresponses[select_cols],
    inds = mpatients[select_cols],
    ct = mcelltypes[select_cols]))
}

show_sc_stat <- function(inds, resp) {
  ## show the numbers of cells in the individuals or
  ## responses groups.

  message("The individuals stats")
  print(table(as.factor(inds)))
  message("The responses stats")
  print(table(as.factor(resp)))
}

sample_cells_per_ind <- function(inds, cellnum = 10) {
  ## for each individual, sample the cells given a cellnum.
  ## index's names are ordered as sort(unique(inds))
  ## return: a vector

  res <- lapply(sort(unique(inds)), FUN = function(name) {
    index <- which(inds == name)
    num <- length(index)
    if (cellnum >= num) {
      invisible(index)
    }
    invisible(sample(x = index, size = cellnum, replace = F))
  })
  invisible(unlist(res, use.names = F))
}

load_pbmc_seurat <- function() {
  invisible(readRDS(here::here("data", "antiPDL1_PBMC_IL8",
    "GSE145281_nicol.RDS")))
}

extract_from_seurat <- function(pbmc_seurat) {
  cntmat <- as.matrix(pbmc_seurat@assays$RNA@counts)
  mresponses <- pbmc_seurat@meta.data$response
  mpatients <- pbmc_seurat@meta.data$patient
  mcelltypes <- pbmc_seurat@meta.data$cellType

  uniqcelltypes <- levels(factor(mcelltypes))
  cellstats <- data.frame(table(factor(mcelltypes)))
  colnames(cellstats) <- c("celltype", "freq")

  invisible(list(cnt = cntmat,
    resp = mresponses,
    inds = mpatients,
    ctyps = mcelltypes,
    uctyps = uniqcelltypes,
    cstat = cellstats))
}
