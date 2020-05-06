library(Seurat)
library(magrittr)
library(ggplot2)
library(corrplot)
import::from(harmony, RunHarmony)
import::from(data.table, data.table)
import::from(mltools, one_hot)
import::from(rstan, stan_rdump)
import::from(dplyr, group_by, top_n)

## * Load scRNAseq data.
gse145281 <- readRDS("../from_avi/20200504/seurat.RDS")

scmeta <- gse145281@meta.data
sccluster <- scmeta$seurat_clusters
# response, i.e., condition
scresp <- scmeta$response
scind <- scmeta$patient
scdata <- gse145281@assays$RNA@counts
# total counts per cell
sc_tcpc <- colSums(as.matrix(scdata))

## * get selected genes
# get deseq.dt object
load("../from_avi/20200504/deseq.dt.RData")
topgnum <- 100
mygenes <- deseq.dt$gene[1:topgnum]

# cytototic T cell
mycluster <- 2
mycells <- which(sccluster == mycluster)

VlnPlot(
  object = gse145281, features = mygenes[1:10],
  group.by = "patient", idents = mycluster
)

## * get harmony delta mean for each individual
## ** extract the correction vectors.
orig_pca <- gse145281@reductions$pca@cell.embeddings[mycells,]
harmony_correct_pca <- gse145281@reductions$harmony@cell.embeddings[mycells,]
delta_pca <- harmony_correct_pca - orig_pca
indhay <- aggregate(delta_pca, list(gse145281@meta.data$patient[mycells]), mean)
colnames(indhay)[1] <- "patient"

## * get data for model iii
modelnm <- "model_v3"

## ** use the individual correction vectors for correlation estimation
nnegcor <- indhay %>%
  .[, -1] %>%
  t() %>%
  cor(., method = "pearson")
nnegcor[which(nnegcor < 0)] <- 0
corrplot(nnegcor)

## ** get counts matrix and design matrix
x_cg <- scdata[mygenes, mycells]
x_ <- sc_tcpc[mycells]
## merge data: merge(x_cg, indhay, by="patient")
ic <- as.matrix(one_hot(data.table(ic = factor(scind[mycells]))))
di <- as.matrix(one_hot(data.table(di = factor(scresp[mycells]))))

## ** set constants
N <- length(mycells)
K <- ncol(ic)
G <- topgnum
J <- 2
scale <- 10000

## ** save data for cmdstan
stan_rdump(c("N", "K", "J", "G","scale", "di", "ic", "x_", "x_cg", "nnegcor"),
  file = paste0("./", modelnm, ".rdump")
)
