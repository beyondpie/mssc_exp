import::from(
  Seurat, DimPlot, VlnPlot, Embeddings,
  RunUMAP, FindNeighbors, FindClusters,
  FindAllMarkers, FindMarkers
)
import::from(magrittr, "%>%", "%<>%")
import::from(harmony, RunHarmony)
import::from(data.table, data.table)
import::from(mltools, one_hot)
import::from(
  ggplot2, ggplot, geom_point,
  geom_smooth, ggsave, geom_violin,
  geom_jitter, geom_boxplot, geom_dotplot, aes,
  position_jitter
)

import::from(
  rstan, stan, extract, read_rdump, read_stan_csv,
  check_hmc_diagnostics, sampling, stan_rdump
)
import::from(dplyr, group_by, top_n)

# From Avi'data, samples cannot find. So use all the gse145281 meta data first.
load("../from_avi/GSE145281.RData")
gse145281 <- gse145281 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

scmeta <- gse145281@meta.data
sccluster <- scmeta$seurat_clusters
# table(sccluster) shows cluster 0 has 6850 cells,
# but only 955 cells not zero on HBB
# cluster 1 has 5237 cells, and 1276 cells are not zero reads on HBB.
mycluster <- 1
gse145281.markers.mycluster <- FindMarkers(gse145281, only.pos = TRUE,
                                           min.pct = 0.25, logfc.threshold = 0.25, 
                                           ident.1 = mycluster)
# gse145281.markers <- FindAllMarkers(gse145281, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- gse145281.markers.mycluster %>% top_n(n = 10, wt = avg_logFC)
# LYZ, CXCL8, GPX1, S100A0, FCN1, S100A12, S100A8, TIMP1

p <- Seurat::DoHeatmap(gse145281, features = top10)

scresp <- scmeta$response
scind <- scmeta$patient
scdata <- gse145281@assays$RNA@counts
# total read counts per individual
sc_tcperi <- colSums(as.matrix(scdata))


# should no DE: HBB, HBA2, HBA1
# should DE: CCL4L1, CCL3L1, CCL3L3
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

draw_sc_gene_data <- function(dt) {
  p <- ggplot(dt, aes(x = ic, y = x_cg, color = ic)) +
    geom_violin(trim = FALSE) +
    geom_jitter(shape = 16, position = position_jitter(0.2))
  return(p)
}


scHBB <- sc_gene_extract(genm = "HBB", mycluster = mycluster)
pHBB <- draw_sc_gene_data(scHBB)
pHBB

scHBA2 <- sc_gene_extract(genm = "HBA2", mycluster = mycluster)
pHBA2 <- draw_sc_gene_data(scHBA2)
pHBA2

scHBA1 <- sc_gene_extract(genm = "HBA1", mycluster = mycluster)
pHBA1 <- draw_sc_gene_data(scHBA1)
pHBA1

scCXCL8 <- sc_gene_extract(genm = "CXCL8", mycluster = mycluster)
pCXCL8 <- draw_sc_gene_data(scCXCL8)
pCXCL8

scGPX1 <- sc_gene_extract(genm = "GPX1", mycluster = mycluster)
pGPX1 <- draw_sc_gene_data(scGPX1)
pGPX1

scS100A9 <- sc_gene_extract(genm = "S100A9", mycluster = mycluster)
pS100A9 <- draw_sc_gene_data(scS100A9)
pS100A9

scLYZ <- sc_gene_extract(genm = "LYZ", mycluster = mycluster)
pLYZ <- draw_sc_gene_data(scLYZ)
pLYZ

scCST3 <- sc_gene_extract(genm = "CST3", mycluster = mycluster)
pCST3 <- draw_sc_gene_data(scCST3)
pCST3


# HBB, LYZ, GPX1
