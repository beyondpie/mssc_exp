library(Seurat)
library(magrittr)
library(ggplot2)
import::from(harmony, RunHarmony)
import::from(data.table, data.table)
import::from(mltools, one_hot)
import::from(rstan, stan_rdump)
import::from(dplyr, group_by, top_n)

## * Load scRNAseq data.
gse145281 <- readRDS("../from_avi/20200504/seurat.RDS")

scmeta <- gse145281@meta.data
sccluster <- scmeta$seurat_clusters
# cytotoxic T cell clusters
mycluster <- 2
scresp <- scmeta$response
scind <- scmeta$patient
scdata <- gse145281@assays$RNA@counts
sc_tcperi <- colSums(as.matrix(scdata))

## * get selected genes
# get deseq.dt object
load("../from_avi/20200504/deseq.dt.RData")

cells <- c("SNHG16", "OASL", "NAMPT", "NFKB1", "NA.499")
cells <- deseq.dt$gene[which(deseq.dt$padj < 1e-2)]
VlnPlot(object = gse145281, features = cells,
        group.by = "patient", idents = mycluster)

sc_gene_extract <- function(genm = "HBB", mycluster = 1) {
  mycells <- which(sccluster == mycluster)
  x_cg <- scdata[which(rownames(scdata) == genm), mycells]
  x_ <- sc_tcperi[mycells]

  sc_gd <- data.table(
    x_cg = x_cg,
    x_ = x_,
    ic = factor(scind[mycells]),
    di = factor(scresp[mycells])
  )
  data <- one_hot(sc_gd)
  N <- nrow(data)
  K <- ncol(data) - 2 - 2
  scale <- 10000
  J <- 2
  ic <- as.matrix(data[, 3:(K + 2)])
  di <- as.matrix(data[, (ncol(data) - 1):ncol(data)])
  stan_rdump(c("N", "K", "J", "scale", "di", "ic", "x_", "x_cg"),
    file = paste0("./sc", genm, ".rdump")
  )
  return(sc_gd)
}
