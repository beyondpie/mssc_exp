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
top_rank <- 100
cells <- deseq.dt$gene[top_rank]
VlnPlot(object = gse145281, features = cells[1:10],
        group.by = "patient", idents = mycluster)

# * get harmony delta mean for each individual
orig_pca <- gse145281@reductions$pca@cell.embeddings
harmony_correct_pca <- gse145281@reductions$harmony@cell.embeddings
delta_pca <- harmony_correct_pca - orig_pca

indhay <- aggregate(delta_pca, list(gse145281@meta.data$patient), mean)
colnames(indhay)[1] <- "patient"
## rownames(indhay) <- indhay$Group.1
## indhay  <- indhay[, -1]


sc_gene_extract <- function(genm = "HBB", mycluster = 1) {
  mycells <- which(sccluster == mycluster)
  x_cg <- scdata[which(rownames(scdata) == genm), mycells]
  x_ <- sc_tcperi[mycells]

  indhays <- merge(x_cg, indhay, by="patient")

  sc_gd <- data.table(
    x_cg = x_cg,
    x_ = x_,
    ic = factor(scind[mycells]),
    di = factor(scresp[mycells]),
    indhc = indhays[, 2:22]
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
