library(Seurat)
library(magrittr)
library(ggplot2)
library(corrplot)
library(Matrix)
library(data.table)
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
topgnum <- 2000
mygenes <- deseq.dt$gene[1:topgnum]

## * select cell cluster
# cytototic T cell
mycluster <- 2
mycells <- which(sccluster == mycluster)

## * PCA analysis for genes in the cell cluster
## normalized data
## mynmldata <- gse145281@assays$RNA@data[mygenes, mycells] %>% as.data.frame()

## use all the genes for PCA
mynmldata <- gse145281@assays$RNA@data[, mycells] %>% as.data.frame()
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

## *** analyze PCA per individual
p <- 20
mycentds <- lapply(inds, function(x) {
  sweep(mygroups[[x]], 1, gmeans[[x]])
})

names(mycentds) <- inds

sVDs <- lapply(inds, function(x) {
  svd(mycentds[[x]])
})

contribs <- lapply(sVDs, function(x) {
  sum(x$d[1:p]) / sum(x$d)
})

names(contribs) <- inds

## RowVar <- function(x, ...) {
##   rowSums((x - rowMeans(x, ...))^2, ...) / (dim(x)[2] - 1)
## }

## RowVar(A)
## *** center data per dividual and then merge.
mycentd <- lapply(inds, function(x) {
  sweep(mygroups[[x]], 1, gmeans[[x]])
}) %>% do.call(cbind, .)
s <- svd(mycentd)

p <- 2
B <- mycentd %>%
  as.matrix() %>%
  `%*%`(., s$v[, 1:p])
B <- B[mygenes, ]

## * summarize data for stan.
modelnm <- "v1"
## * Redefine genes
## use top rank
utr <- 8
## use low rank
ulr <- 6
goldgenes <- data.frame(genes = c(
  "SNHG16",
  "OASL",
  "NAMPT",
  "NFKB1",
  "BCL2L11",
  "IRF8",
  "TPM4",
  "TRAF4",
  "ICAM1", "XCL2", "XCL1",
  "RPS26P11", "LOC101929876", "LOC100996747",
  "HBA1", "HBA2", "HBB", "HBD",
  "CCL3L3", "CCL3L1", "CCL3",
  "KDM6A",
  "ZNF721",
  "HDDC2",
  "YIPF5",
  "MAK16",
  "TOX"
  ), module = c(seq(1, utr), rep(utr+1, 3),
                 rep(utr+2, 3), rep(utr+3, 4),
                 rep(utr+4, 3), seq(utr+5, utr+5 + ulr-1)))
module <- factor(goldgenes$module)
GoldB <- as.matrix(one_hot(data.table(module)))
rownames(GoldB) <- goldgenes$genes

## ** get counts matrix and design matrix
mygenes <- goldgenes$genes
B <- GoldB

Xcg <- t(as.matrix(scdata[mygenes, mycells]))
## IXcg <- Xcg
Xgc <- as.matrix(scdata[mygenes, mycells])
S <- sc_tcpc[mycells]
## merge data: merge(x_cg, indhay, by="patient")
myinds <- gse145281@meta.data$patient[mycells]
inds <- attr(factor(myinds), "levels")
ic <- factor(scind[mycells], levels = inds)
XInd <- as.matrix(one_hot(data.table(ic = ic)))
## IXInd <- as.numeric(ic)

di <- factor(scresp[mycells], levels=c(0,1))
XCond <- as.matrix(one_hot(data.table(di = di)))
## IXCond <- as.numeric(di)

## ** set constants
N <- length(mycells)
K <- ncol(XInd)
G <- length(mygenes)
J <- ncol(XCond)
P <- ncol(B)

## ** save data for cmdstan
stan_rdump(c(
  "N", "K", "J", "G", "XCond", "XInd","S", "Xcg", "B", "P"),
  file = paste0("./", "my27gene18module",".rdump"))

## ** for model 1
B <- matrix(0.0, nrow=G, ncol=1)
P <- ncol(B)
stan_rdump(c(
  "N", "K", "J", "G", "XCond", "XInd","S", "Xcg", "B", "P"),
  file = paste0("./", modelnm,"_", G,".rdump"))
