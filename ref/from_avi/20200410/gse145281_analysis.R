library(Matrix)
library(Seurat)
library(data.table)
library(magrittr)
library(cowplot)
library(harmony)
library(uwot)
library(parallel)
library(ggplot2)

load("gse145281.RData")
GSE145281.meta <- data.table(gse145281@meta.data)

GSE145281.meta.match <- GSE145281.meta[match(rownames(gse145281@meta.data), samples)]
gse145281@meta.data %<>% cbind(., GSE145281.meta.match[, .(patient, response)]) %>%
  as.data.frame()

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = gse145281, reduction = "pca", pt.size = .1, group.by = "patient", do.return = TRUE)
p2 <- VlnPlot(object = gse145281, features = "PC_1", group.by = "patient", do.return = TRUE, pt.size = .1)
pdf("gse145281.without.harmony.pdf")
# plot_grid(p1,p2)
print(p1)
dev.off()

gse145281 <- gse145281 %>%
  RunHarmony("patient")


harmony_embeddings <- Embeddings(gse145281, "harmony")
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = gse145281, reduction = "harmony", pt.size = .1, group.by = "patient", do.return = TRUE)
p2 <- VlnPlot(object = gse145281, features = "harmony_1", group.by = "patient", do.return = TRUE, pt.size = .1)

pdf(".figs/gse145281.harmony.pdf")
# plot_grid(p1,p2)
print(p1)
dev.off()

gse145281[["umapunorm"]] <- gse145281@reductions$pca@cell.embeddings %>%
  umap(., n_threads = 2, metric = "cosine", pca = NULL, n_neighbors = 50) %>%
  set_rownames(colnames(gse145281)) %>%
  set_colnames(paste0("umapunorm_", 1:2)) %>%
  CreateDimReducObject(embeddings = ., key = "umapunorm_", assay = DefaultAssay(gse145281))

DimPlot(gse145281, reduction = "umapunorm", group.by = "patient", pt.size = 0.5) %>%
  ggsave(filename = "./gse145281.umap.without.harmony.pdf", plot = .)

gse145281 <- gse145281 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

DimPlot(gse145281, reduction = "umap", group.by = "patient", pt.size = .1, split.by = "patient") %>%
  ggsave("./gse145281.harmony.umap.patient.pdf", plot = .)

DimPlot(gse145281, reduction = "umap", group.by = "patient", label = TRUE, pt.size = .1) %>%
  ggsave("./gse145281.harmony.umap.pdf", plot = .)

DimPlot(gse145281, reduction = "umap", label = TRUE, pt.size = .1) %>%
  ggsave("./gse145281.harmony.umap.cluster.pdf", plot = .)

##
p1 <- DimPlot(gse145281, reduction = "umap", group.by = "response", pt.size = .1, split.by = "response")
pdf(".figs/gse145281.harmony.umap.response.pdf")
plot_grid(p1)
dev.off()


## response

if (F) {
  library(dplyr)
  gse145281.markers <- FindAllMarkers(gse145281, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  gse145281.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_logFC)

  top10 <- gse145281.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_logFC)
  p <- DoHeatmap(gse145281, features = top10$gene) + NoLegend()
  ggsave(file = ".figs/gse145281.harmony.markers.pdf", p, width = 15, height = 10)

  p <- FeaturePlot(gse145281, features = c(
    "MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A", "GZMB", "PRF1", "IL7R", "CCR7", "CD4"
  ))
  ggsave(file = ".figs/gse145281.harmony.features.pdf", p, width = 15, height = 10)
}
# find differential expression between responders and non responder in a cell population.
# For a gene differntially expressed in responders show that it is specific to a patient i.e.  R1 vs. rest-of-responders
# argue that this is a confounding effect and need to be corrected in every analysis
# there is need a tool for this
## cluster 6 is B-cells

out <- gse145281@meta.data %>%
  data.table() %>%
  .[, .N, by = list(seurat_clusters, patient)] %>%
  dcast(data = ., formula = seurat_clusters ~ patient, value.var = "N", fun.aggregate = sum)


# clust 1 confounders
gse145281.clust <- gse145281[, gse145281$seurat_clusters == 1]
markers <- FindMarkers(gse145281.clust, ident.1 = "0", group.by = "response", min.pct = 0.05, logfc.threshold = 0.15)
## confounding markers
gse145281.confound <- gse145281.clust[, gse145281.clust$response == 0]
markers.confound <- FindMarkers(gse145281.confound, ident.1 = "NR1", group.by = "patient", min.pct = 0.05, logfc.threshold = 0.15)
markers.intersect <- intersect(rownames(markers[markers$avg_logFC > 0, ])[1:50], rownames(markers.confound)[1:200])
VlnPlot(gse145281.clust, markers.intersect, group.by = "patient") %>%
  ggsave(file = ".figs/gse145281.clust1.markers.pdf", ., width = 15, height = 10)
markers.common <- intersect(rownames(markers), rownames(markers.confound))
markers.merge <- markers.confound %>%
  set_colnames(paste0(colnames(markers.confound), ".nr1")) %>%
  .[markers.common, ] %>%
  cbind(., markers[markers.common, ])
markers.merge$genes <- rownames(markers.merge)
markers.merge <- data.table(markers.merge)
markers.merge[, avg_logFC_r := rank(avg_logFC)][, avg_logFC.nr1_r := rank(avg_logFC.nr1)]
library(ggpubr)
p <- ggplot(markers.merge, aes(x = avg_logFC, y = avg_logFC.nr1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman") +
  xlab("R vs. NR") +
  ylab("NR1 vs. other-NR")
ggsave(file = ".figs/gse145281.clust1.confound.pdf", p)


## implement the control for patients
## https://www.biorxiv.org/content/10.1101/2020.01.15.906248v1.full.pdf
# make sure there is issue with
one_hot.mat <- gse145281@meta.data %>%
  data.table(ID = rownames(.)) %>%
  .[, list(ID, patiient = factor(patient))] %>%
  mltools::one_hot(.)
gse145281@meta.data %<>% cbind(., one_hot.mat[, -1, with = F]) %>%
  as.data.frame()
patient.latent.vars <- colnames(one_hot.mat)[-1]

# clust 2
gse145281.clust <- gse145281[, gse145281$seurat_clusters == 2]
markers <- FindMarkers(gse145281.clust, ident.1 = "0", group.by = "response", min.pct = 0.25, logfc.threshold = 0.15)
markers.mast <- FindMarkers(gse145281.clust, ident.1 = "0", group.by = "response", min.pct = 0.25, logfc.threshold = 0.15, test.use = "MAST")
markers.mast.patient <- FindMarkers(gse145281.clust, ident.1 = "0", group.by = "response", min.pct = 0.25, logfc.threshold = 0.15, test.use = "MAST", latent.vars = patient.latent.vars)
gse145281.confound <- gse145281.clust[, gse145281.clust$response == 1]
markers.confound <- FindMarkers(gse145281.confound, ident.1 = "R3", group.by = "patient", min.pct = 0.25, logfc.threshold = 0.15)
markers.intersect <- intersect(rownames(markers[markers$avg_logFC > 0, ])[1:50], rownames(markers.confound)[1:30])
VlnPlot(gse145281.clust, markers.intersect, group.by = "patient") %>%
  ggsave(file = ".figs/gse145281.clust2.markers.pdf", ., width = 15, height = 10)

markers.common <- intersect(rownames(markers), rownames(markers.confound))
markers.merge <- markers.confound %>%
  set_colnames(paste0(colnames(markers.confound), ".r3")) %>%
  .[markers.common, ] %>%
  cbind(., markers[markers.common, ])
markers.merge$genes <- rownames(markers.merge)
markers.merge <- data.table(markers.merge)
markers.merge[, avg_logFC_r := rank(avg_logFC)][, avg_logFC.r3_r := rank(avg_logFC.r3)]
library(ggpubr)
p <- ggplot(markers.merge, aes(x = avg_logFC, y = avg_logFC.r3)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman") +
  xlab("R vs. NR") +
  ylab("R3 vs. other-NR")
ggsave(file = ".figs/gse145281.clust2.confound.pdf", p)

##
markers.common <- intersect(rownames(markers.mast.patient), rownames(markers.confound))
markers.merge <- markers.confound %>%
  set_colnames(paste0(colnames(markers.confound), ".r3")) %>%
  .[markers.common, ] %>%
  cbind(., markers.mast.patient[markers.common, ])
markers.merge$genes <- rownames(markers.merge)
markers.merge <- data.table(markers.merge)
# markers.merge[,avg_logFC_r:=rank(avg_logFC)][,avg_logFC.r3_r:=rank(avg_logFC.r3)]
markers.merge[, deg := log(p_val + 1E-30) * sign(avg_logFC)] %>%
  .[, deg.r3 := log(p_val.r3 + 1E-30) * sign(avg_logFC.r3)]

p <- ggplot(markers.merge, aes(x = deg, y = deg.r3)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman") +
  xlab("R vs. NR") +
  ylab("R3 vs. other-NR")
ggsave(file = ".figs/gse145281.clust2.pseudoreplication.pdf", p)

## Proposing the solution of patient specific confounding effects.
# 1. Leverage cluster level in DEG. The cluster DEG could be used to
# Y ~ Beta|cluster + 1|patient 1|cluster

genes.curr <- gse145281[["RNA"]]@var.features
exp.curr <- gse145281[["RNA"]][genes.curr, ]
meta.dt <- gse145281@meta.data %>%
  data.table(ID = rownames(.)) %>%
  .[, .(ID, response, patient = as.factor(patient), cluster = as.factor(seurat_clusters))]


eval.lmer.solution <- function(dt.curr, val) {
  tryCatch(
    {
      dt.curr <- dt.curr %>%
        mutate(exp = val)
      aa <- lme4::lmer(data = dt.curr, exp ~ response + (1 + response | cluster) + (1 | patient))
      bb <- coef(aa)
      c(unlist(bb$cluster), bb$patient[, 2])
    },
    error = function(e) rep(NA, 42)
  )
}

library(parallel)
out <- mclapply(seq(nrow(exp.curr)), function(tt) {
  eval.lmer.solution(meta.dt, exp.curr[tt, ])
},
mc.cores = 48
)

out.dt <- do.call(rbind, out)

cluster.markers <- (out.dt[, 1:16] - rowMeans(out.dt[, 1:16])) %>%
  set_colnames(paste0("cluster_", 0:15)) %>%
  data.table() %>%
  mutate(gene = rownames(exp.curr)) %>%
  data.table()
cluster.markers <- cluster.markers[order(cluster_0, decreasing = T)]

responder.markers <- (out.dt[, 17:32]) %>%
  set_colnames(paste0("cluster_", 0:15)) %>%
  data.table() %>%
  mutate(gene = rownames(exp.curr)) %>%
  data.table()
head(responder.markers[order(cluster_2, decreasing = T)]$gene)

gse145281.clust <- gse145281[, gse145281$seurat_clusters == 2] %>%
  .[gse145281[["RNA"]]@var.features, ]
markers <- FindMarkers(gse145281.clust, ident.1 = "0", group.by = "response", min.pct = -1, logfc.threshold = -1)
gse145281.confound <- gse145281.clust[, gse145281.clust$response == 1]
markers.confound <- FindMarkers(gse145281.confound, ident.1 = "R3", group.by = "patient", min.pct = -1, logfc.threshold = -1)
markers.common <- intersect(rownames(markers), rownames(markers.confound))
markers.merge <- markers.confound %>%
  set_colnames(paste0(colnames(markers.confound), ".r3")) %>%
  .[markers.common, ] %>%
  cbind(., markers[markers.common, ])
markers.merge$genes <- rownames(markers.merge)
markers.merge <- data.table(markers.merge)
markers.merge[, avg_logFC_r := rank(avg_logFC)][, avg_logFC.r3_r := rank(avg_logFC.r3)]
markers.merge$model1 <- responder.markers[match(markers.merge$genes, gene)]$cluster_2

library(ggpubr)
p <- ggplot(markers.merge, aes(x = -model1, y = avg_logFC.r3)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman") +
  xlab("model1 (R vs. NR)") +
  ylab("R3 vs. other-NR")
ggsave(file = ".figs/gse145281.clust2.model1.pdf", p)

p <- ggplot(markers.merge, aes(x = -model1, y = avg_logFC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman") +
  xlab("model1 (R vs. NR)") +
  ylab("Seurat (R vs. NR)")
ggsave(file = ".figs/gse145281.clust2.model1.comparison.seurat.pdf", p)

markers.deseq2 <- FindMarkers(gse145281.clust, ident.1 = "0", group.by = "response", min.pct = -1, logfc.threshold = -1, test.use = "DESeq2")


markers$avg_logFC.effect <- -markers$avg_logFC

library(EnhancedVolcano)
EnhancedVolcano(markers,
  lab = rownames(markers),
  x = "avg_logFC.effect",
  y = "p_val_adj",
  pCutoff = 1e-2,
  FCcutoff = .5,
  # ylim = c(0,5.2),
  # xlim = c(-1, 1),
  # pointSize = 4.0,
  pointSize = c(ifelse(markers$p_val_adj < 1E-10, 1, 0.2)),
  labSize = 4.0,
  legend = c(
    "NS", "Log (base 2) fold-change", "Adj.P value",
    "Adj.P value & Log (base 2) fold-change"
  ),
  legendPosition = "right",
  legendLabSize = 8,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  colAlpha = 0.8,
  colConnectors = "grey30"
) %>%
  ggsave(".figs/gse145281.clust2.volcano.wilcox.pdf", .)




## Using the psuedo-bulk as prior (model4)

## run deseq for each cluster

deseq.all.out <- lapply(sort(unique(gse145281$seurat_clusters)), function(tt) {
  tryCatch(
    {
      gse145281.clust <- gse145281[, gse145281$seurat_clusters == tt]
      patient.cell.one.hot <- gse145281.clust@meta.data %>%
        data.table(ID = rownames(.)) %>%
        .[, .(ID, patient = as.factor(patient))] %>%
        .[, .(ID, patient)] %>%
        mltools::one_hot(.) %>%
        .[, -1, with = F] %>%
        as.matrix()

      psuedo.bulk <- gse145281.clust@assays$RNA@counts %*% patient.cell.one.hot

      cells.1 <- grep("_R", colnames(psuedo.bulk), value = T)
      cells.2 <- grep("_NR", colnames(psuedo.bulk), value = T)
      deseq.out <- DESeq2DETest(data.use = psuedo.bulk, cells.1 = cells.2, cells.2 = cells.1)
      deseq.dt <- deseq.out %>%
        as.data.frame() %>%
        mutate(gene = rownames(.)) %>%
        data.table()
      deseq.dt
    },
    error = function(e) NA
  )
})

names(deseq.all.out) <- sort(unique(gse145281$seurat_clusters))
# remove NA cluster data from deseq
is.notincluded <- names(deseq.all.out)[which(is.na(deseq.all.out))]
is.included <- setdiff(names(deseq.all.out), is.notincluded)
samples.included <- which(!(gse145281$seurat_clusters %in% is.notincluded))

psuedobulk.lfc <- lapply(is.included, function(tt) deseq.all.out[[tt]]$log2FoldChange) %>%
  do.call(cbind, .)
psuedobulk.lfc[is.na(psuedobulk.lfc)] <- 0

gse145281.curr <- gse145281[, samples.included]
meta.model4.dt <- gse145281.curr@meta.data %>%
  data.table(ID = rownames(.)) %>%
  .[, .(ID, response,
    patient = as.factor(patient), cluster = as.factor(seurat_clusters),
    nCount_RNA = avinash::normalize.std(nCount_RNA),
    nFeature_RNA = avinash::normalize.std(nFeature_RNA)
  )]
meta.model4.dt.curr <- meta.model4.dt[, .(ID, cluster)] %>%
  mltools::one_hot(.) %>%
  .[, -1, with = F] %>%
  cbind(meta.model4.dt, .)

library(lme4)
eval.lmer.psuedobulk.priors <- function(dt.curr, val, fmla, psuedobulk.priors) {
  tryCatch(
    {
      val.res <- lm(val ~ psuedobulk.priors) %>%
        stats::residuals()
      dt.curr <- dt.curr %>%
        mutate(GENE = val.res)
      aa <- lmer(data = dt.curr, fmla)
      preds <- predict(aa, newdata = dt.curr)
      preds
    },
    error = function(e) NA
  )
}
latent.var.names <- c(
  grep("cluster_", colnames(meta.model4.dt.curr), value = T),
  "nCount_RNA", "nFeature_RNA"
)
fmla <- as.formula(object = paste(
  "GENE ~",
  paste(latent.var.names[-1], collapse = "+"),
  "+(1|patient)"
))

eval.lmer.psuedobulk.priors2 <- function(dt.curr, val, fmla, psuedobulk.priors) {
  tryCatch(
    {
      dt.curr <- dt.curr %>%
        mutate(GENE = val) %>%
        mutate(psuedobulk.priors = psuedobulk.priors)
      aa <- lmer(data = dt.curr, fmla)
      dt.curr$psuedobulk.priors <- 0
      preds <- predict(aa, newdata = dt.curr)
      preds
    },
    error = function(e) NA
  )
}
latent.var.names <- c(
  grep("cluster_", colnames(meta.model4.dt.curr), value = T),
  "nCount_RNA", "nFeature_RNA", "psuedobulk.priors"
)
fmla <- as.formula(object = paste(
  "GENE ~",
  paste(latent.var.names[-1], collapse = "+"),
  "+(1|patient)"
))

genes.curr <- gse145281.curr[["RNA"]]@var.features
exp.included.curr <- gse145281.curr[["RNA"]][genes.curr, ]
deseq.gene.order <- deseq.all.out[[1]]$gene
nonresponders.inx <- which(meta.model4.dt$response == 0)
out.model4.models <- mclapply(seq(nrow(exp.included.curr)), function(tt) {
  gene.curr <- rownames(exp.included.curr)[tt]
  psuedobulk.lfc.curr <- psuedobulk.lfc[deseq.gene.order == gene.curr, ]
  psuedobulk.lfc.matched <- psuedobulk.lfc.curr[match(meta.model4.dt.curr$cluster, is.included)]
  psuedobulk.lfc.matched[nonresponders.inx] <- 0

  eval.lmer.psuedobulk.priors2(meta.model4.dt.curr, exp.included.curr[tt, ], fmla, psuedobulk.priors = psuedobulk.lfc.matched)
},
mc.cores = 16
)


preds.patients.cluster <- do.call(rbind, out.model4.models)

exp.curr.clust <- exp.curr[, meta.model4.dt.curr$cluster == 2]
preds.patients.cluster.clust <- preds.patients.cluster[, meta.model4.dt.curr$cluster == 2]
response.clust <- meta.model4.dt[cluster == 2]$response
out.model4.final <- mclapply(seq(nrow(exp.curr)), function(tt) {
  tryCatch(
    {
      aa <- lm(exp.curr.clust[tt, ] ~ response.clust + preds.patients.cluster.clust[tt, ])
      summary(aa)$coefficients["response.clust", ]
    },
    error = function(e) NA
  )
}, mc.cores = 48)
out.model4.final.clust2 <- do.call(rbind, out.model4.final) %>%
  data.table() %>%
  set_colnames(c("estimate", "se", "t", "P")) %>%
  mutate(p.adj = p.adjust(P)) %>%
  mutate(gene = rownames(exp.curr)) %>%
  data.table() %>%
  .[order(P)]

markers.merge$model4 <- out.model4.final.clust2[match(markers.merge$genes, gene)]$estimate

library(ggpubr)
p <- ggplot(markers.merge, aes(x = -model4, y = avg_logFC.r3)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman") +
  xlab("model4 (R vs. NR)") +
  ylab("R3 vs. other-NR")
ggsave(file = ".figs/gse145281.clust2.model4.pdf", p)

p <- ggplot(markers.merge, aes(x = -model4, y = avg_logFC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman") +
  xlab("model4 (R vs. NR)") +
  ylab("Seurat (R vs. NR)")
ggsave(file = ".figs/gse145281.clust2.model4.comparison.seurat.pdf", p)

genes.curr <- out.model4.final.clust2$gene[1:16]
VlnPlot(gse145281.clust, genes.curr, group.by = "patient") %>%
  ggsave(file = ".figs/gse145281.clust2.model4.vlnplot.pdf", ., width = 15, height = 10)

# create psuedo.bulk expression
sco.curr <- gse145281[gse145281[["RNA"]]@var.features, gse145281$seurat_clusters == 2]

exp.curr.clust <- sco.curr[["RNA"]]@data
meta.dt.model5 <- sco.curr@meta.data
deseq.dt.model5 <- deseq.all.out[[3]] %>%
  .[match(rownames(sco.curr), gene)]
deseq.dt.model5[, factor := ifelse(abs(stat) < 4.1, 0, 1 - 3 / abs(stat))]
deseq.factor.curr <- deseq.dt.model5$factor

bulk.patient.order <- unique(meta.dt.model5$patient)
psuedobulk.exp.curr <- lapply(bulk.patient.order, function(tt) {
  rowMeans(exp.curr.clust[, meta.dt.model5$patient == tt], na.rm = T)
}) %>%
  do.call(cbind, .)
psuedobulk.response <- grepl("^R", bulk.patient.order) + 0

eval.lmer.psuedobulk.priors3 <- function(dt.curr, val, psuedobulk.val, psuedobulk.response, deseq.factor) {
  tryCatch(
    {
      dt1 <- data.table(psuedobulk.val = psuedobulk.val, psuedobulk.response = psuedobulk.response)
      m1 <- lm(psuedobulk.val ~ psuedobulk.response, data = dt1)
      predicted <- predict(m1, newdata = dt1) * deseq.factor
      predicted.matched <- predicted[dt.curr$map]
      dt.curr <- dt.curr %>%
        mutate(GENE = val) %>%
        mutate(GENE.res = val - predicted.matched) %>%
        mutate(preds.bulk = predicted.matched)

      aa <- lmer(data = dt.curr, GENE ~ 1 | patient)
      dt.curr$preds.bulk <- 0
      preds <- predict(aa, newdata = dt.curr)
      dt.curr <- dt.curr %>%
        mutate(GENE = val) %>%
        mutate(preds = preds)
      aa <- lm(GENE ~ response + preds, data = dt.curr)
      # bb = lm(GENE ~  preds, data = dt.curr)
      # anova(aa, bb)
      summary(aa)$coefficients["response", ]
    },
    error = function(e) NA
  )
}


meta.dt.model5$map <- match(meta.dt.model5$patient, bulk.patient.order)

out.model5.models <- mclapply(seq(nrow(exp.curr.clust)), function(tt) {
  psuedobulk.val <- psuedobulk.exp.curr[tt, ]
  eval.lmer.psuedobulk.priors3(meta.dt.model5, exp.curr.clust[tt, ], psuedobulk.val, psuedobulk.response = psuedobulk.response, deseq.factor.curr[tt])
},
mc.cores = 48
)

out.model5.final.clust2 <- do.call(rbind, out.model5.models) %>%
  data.table() %>%
  set_colnames(c("estimate", "se", "t", "P")) %>%
  mutate(p.adj = p.adjust(P)) %>%
  mutate(gene = rownames(exp.curr.clust)) %>%
  data.table() %>%
  .[order(P)]

EnhancedVolcano::EnhancedVolcano(deseq.dt.model5,
  lab = deseq.dt.model5$gene,
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 5E-4,
  FCcutoff = 1,
  # ylim = c(0,5.2),
  # xlim = c(-1, 1),
  # pointSize = 4.0,
  pointSize = c(ifelse(deseq.dt.model5$pvalue < 1E-2, 1, 0.2)),
  labSize = 4.0,
  legend = c(
    "NS", "Log (base 2) fold-change", "P value",
    "P value & Log (base 2) fold-change"
  ),
  legendPosition = "right",
  legendLabSize = 8,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  colAlpha = 0.8,
  colConnectors = "grey30"
) %>%
  ggsave(".figs/gse145281.clust2.volcano.deseq.pdf", .)




eval.lmer.psuedobulk.priors4 <- function(dt.curr, val, psuedobulk.val, psuedobulk.response, deseq.factor) {
  tryCatch(
    {
      dt1 <- data.table(psuedobulk.val = psuedobulk.val, psuedobulk.response = psuedobulk.response)
      m1 <- lm(psuedobulk.val ~ psuedobulk.response, data = dt1)
      predicted <- predict(m1, newdata = dt1) * deseq.factor
      predicted.matched <- predicted[dt.curr$map]
      dt.curr <- dt.curr %>%
        mutate(GENE = val) %>%
        mutate(GENE.res = val - predicted.matched) %>%
        mutate(preds.bulk = predicted.matched)

      aa <- lmer(data = dt.curr, GENE.res ~ 1 | patient)
      dt.curr$preds.bulk <- 0
      preds <- predict(aa, newdata = dt.curr)
      dt.curr <- dt.curr %>%
        mutate(GENE = val) %>%
        mutate(preds = preds)
      aa <- lm(GENE ~ response + preds, data = dt.curr)
      # bb = lm(GENE ~  preds, data = dt.curr)
      # anova(aa, bb)
      summary(aa)$coefficients["response", ]
    },
    error = function(e) NA
  )
}


out.model6.models <- mclapply(seq(nrow(exp.curr.clust)), function(tt) {
  psuedobulk.val <- psuedobulk.exp.curr[tt, ]
  eval.lmer.psuedobulk.priors4(meta.dt.model5, exp.curr.clust[tt, ], psuedobulk.val, psuedobulk.response = psuedobulk.response, deseq.factor.curr[tt])
},
mc.cores = 48
)

out.model6.final.clust2 <- do.call(rbind, out.model6.models) %>%
  data.table() %>%
  set_colnames(c("estimate", "se", "t", "P")) %>%
  mutate(p.adj = p.adjust(P)) %>%
  mutate(gene = rownames(sco.curr)) %>%
  data.table() %>%
  .[order(P)]

EnhancedVolcano::EnhancedVolcano(out.model6.final.clust2,
  lab = out.model6.final.clust2$gene,
  x = "estimate",
  y = "P",
  pCutoff = 1e-2,
  FCcutoff = .08,
  # ylim = c(0,5.2),
  # xlim = c(-1, 1),
  # pointSize = 4.0,
  pointSize = c(ifelse(out.model6.final.clust2$P < 1E-2, 2, 0.2)),
  labSize = 4.0,
  legend = c(
    "NS", "Log (base 2) fold-change", "P value",
    "P value & Log (base 2) fold-change"
  ),
  legendPosition = "right",
  legendLabSize = 8,
  legendIconSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  colAlpha = 0.8,
  colConnectors = "grey30"
) %>%
  ggsave(".figs/gse145281.clust2.volcano.model6.pdf", .)

genes.curr <- out.model6.final.clust2[P < 1E-2]$gene
VlnPlot(sco.curr, genes.curr, group.by = "patient") %>%
  ggsave(file = ".figs/gse145281.clust2.model6.vlnplot.pdf", ., width = 15, height = 10)

markers.merge$model6 <- out.model6.final.clust2[match(markers.merge$genes, gene)]$estimate

library(ggpubr)
p <- ggplot(markers.merge, aes(x = -model6, y = avg_logFC.r3)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  xlab("model6 (R vs. NR)") +
  ylab("R3 vs. other-R")
ggsave(file = ".figs/gse145281.clust2.model6.pdf", p)

p <- ggplot(markers.merge, aes(x = -model6, y = -avg_logFC)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  xlab("model6 (R vs. NR)") +
  ylab("Seurat (R vs. NR)")
ggsave(file = ".figs/gse145281.clust2.model6.comparison.seurat.pdf", p)



EnhancedVolcano::EnhancedVolcano(,
  lab = out.model6.final.clust2$gene,
  x = "estimate",
  y = "P",
  pCutoff = 1e-2,
  FCcutoff = .08,
  # ylim = c(0,5.2),
  # xlim = c(-1, 1),
  # pointSize = 4.0,
  pointSize = c(ifelse(out.model6.final.clust2$P < 1E-2, 2, 0.2)),
  labSize = 4.0,
  legend = c(
    "NS", "Log (base 2) fold-change", "P value",
    "P value & Log (base 2) fold-change"
  ),
  legendPosition = "right",
  legendLabSize = 8,
  legendIconSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  colAlpha = 0.8,
  colConnectors = "grey30"
) %>%
  ggsave(".figs/gse145281.clust2.volcano.model6.pdf", .)


sco.curr.confound <- sco.curr[, (sco.curr$seurat_clusters == 2 & sco.curr$response == 0)]

for (pat in unique(sco.curr.confound$patient)) {
  label <- sprintf("%s_confound", pat)
  aa <- FindMarkers(sco.curr.confound, ident.1 = pat, group.by = "patient", min.pct = -1, logfc.threshold = -1)
  markers.merge[[label]] <- aa[markers.merge$gene, ]$avg_logFC
}
save(file = "~/liulab_home/data/single_cell/GSE145281/GSE145281.RData", gse145281)

sco.curr.confound <- sco.curr[, (sco.curr$seurat_clusters == 2 & sco.curr$response == 1)]
for (pat in unique(sco.curr.confound$patient)) {
  label <- sprintf("%s_confound", pat)
  aa <- FindMarkers(sco.curr.confound, ident.1 = pat, group.by = "patient", min.pct = -1, logfc.threshold = -1)
  markers.merge[[label]] <- aa[markers.merge$gene, ]$avg_logFC
}

library(ggpubr)
avi.dt <- data.table(x = markers.merge$avg_logFC, y = unlist(markers.merge[, grep("confound", colnames(markers.merge)), with = F]), patient = rep(c(paste0("NR", 1:5), paste0("R", 1:5)), each = 2000))
p <- ggplot(avi.dt, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman") +
  xlab("model1 (R vs. NR)") +
  ylab("NR vs. other-NR / R vs. other-R") +
  facet_wrap(~patient)

ggsave(file = ".figs/gse145281.seurat.clust2.seurat.pdf", p)


avi.dt <- data.table(x = markers.merge$model6, y = unlist(markers.merge[, grep("confound", colnames(markers.merge)), with = F]), patient = rep(c(paste0("NR", 1:5), paste0("R", 1:5)), each = 2000))
p <- ggplot(avi.dt, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  xlab("model6 (R vs. NR)") +
  ylab("NR vs. other-NR / R vs. other-R") +
  facet_wrap(~patient)

ggsave(file = ".figs/gse145281.clust2.model6.all.pdf", p)
