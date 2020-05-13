library(Seurat)
library(magrittr)
library(ggplot2)
library(corrplot)
library(Matrix)
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
topgnum <- 1000
mygenes <- deseq.dt$gene[1:topgnum]

## * select cell cluster
# cytototic T cell
mycluster <- 2
mycells <- which(sccluster == mycluster)

## * PCA analysis for genes in the cell cluster
## normalized data
mynmldata <- gse145281@assays$RNA@data[mygenes, mycells] %>% as.data.frame()
myinds <- gse145281@meta.data$patient[mycells]

## ** center mean per individual
old_colnm <- colnames(mynmldata)
colnames(mynmldata) <- myinds

## not scale since we hope matrix can involve
## different count scale in matrix B or W.

inds <- attr(factor(myinds), "levels")
mygroups <- sapply(inds, function(x) {
  mynmldata[startsWith(names(mynmldata), x)]
}, simplify = FALSE)
## gene means across cells per individual
gmeans <- lapply(mygroups, rowMeans) %>% as.data.frame()
## center data per dividual and then merge.
mycentd <- lapply(inds, function(x) {
  sweep(mygroups[[x]], 1, gmeans[[x]])
}) %>% do.call(cbind, .)

s <- svd(mycentd)

p <- 20
myB <- mycentd %>% as.matrix %>% `%*%`(., s$v[, 1:p])

## * summarize data for stan.
modelnm <- "model_v4"
## ** get counts matrix and design matrix
x_cg <- t(as.matrix(scdata[mygenes, mycells]))
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
stan_rdump(c("N", "K", "J", "G", "scale", "di", "ic", "x_", "x_cg", "myB"),
  file = paste0("./", modelnm, ".rdump")
)





