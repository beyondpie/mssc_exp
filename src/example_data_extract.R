import::from(
  Seurat, DimPlot, VlnPlot, Embeddings,
  RunUMAP, FindNeighbors, FindClusters,
  FindMarkers
)
import::from(magrittr, "%>%", "%<>%")
import::from(harmony, RunHarmony)
import::from(data.table, data.table)
import::from(mltools, one_hot)
import::from(
  ggplot2, ggplot, geom_point,
  geom_smooth, ggsave, geom_violin,
  geom_jitter, geom_boxplot, geom_dotplot
)
import::from(lme4, lmer)

# From Avi'data, samples cannot find. So use all the gse145281 meta data first.
load("../from_avi/GSE145281.RData")

scmeta <- gse145281@meta.data
sccluster <- scmeta$seurat_clusters
# table(sccluster) shows cluster 0 has 6850 cells,
# but only 955 cells not zero on HBB
# cluster 1 has 5237 cells, and 1276 cells are not zero reads on HBB.
mycluster <- 1
mycells <- which(sccluster == mycluster)

scresp <- scmeta$response
scind <- scmeta$patient
scdata <- gse145281@assays$RNA@counts

scHBB <- data.table(
  x = scdata[which(rownames(scdata) == "HBB"), mycells],
  i = factor(scind[mycells]),
  d = factor(scresp[mycells])
)

p <- ggplot(scHBB, aes(x = i, y = x, color = i)) +
  geom_violin(trim = FALSE) +
  geom_jitter(shape = 16, position = position_jitter(0.2))

# TODO: t.test,  DESeq and others on scHBB

# * for rstan
scHBBd  <- one_hot(scHBB)
