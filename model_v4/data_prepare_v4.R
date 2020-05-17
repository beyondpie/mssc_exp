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
library(ggdendro)
library(reshape2)
library(grid)
library(gplots)
library(viridis)
## * Load scRNAseq data.
scale <- 10000
gse145281 <- readRDS("../from_avi/20200504/seurat.RDS")
gse145281 <- NormalizeData(object = gse145281,
                           normalization.method = "LogNormalize",
                           scale.factor=scale)

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
topgnum <- 500
mygenes <- deseq.dt$gene[1:topgnum]
## * select cell cluster
# cytototic T cell
mycluster <- 2
mycells <- which(sccluster == mycluster)
## * do clustering to find the gene modules on top ranked genes
genexp <- as.matrix(t(scale(t(gse145281@assays$RNA@data[mygenes , mycells]))))

genexphc <- hclust( as.dist(1 - cor(t(genexp), method="pearson")), method="complete")
# mycl <- cutree(genexphc, h=max(genexphc$height/1.5))
mycl <- cutree(genexphc, h=0.9)

a <- rep(0, max(mycl))
for (i in 1:max(mycl)) {
  a[i] = length(which(mycl == i))
}
clusters = which(a > 2)
# examine the cluster membership by it's order
# in the heatmap
# mycl[hr$order]
# grab a cluster
# cluster1 <- d[mycl == 1,]
# or simply add the cluster ID to your data
# foo <- cbind(d, clusterID=mycl)
# examine the data with cluster ids attached, and ordered like the heat map
# foo[hr$order,]
# too slow, and not looks well
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- rev(redgreen(75))
# heatmap.2(genexp, main="Hierarchical Cluster", Rowv=as.dendrogram(genexphc), 
#          Colv=NA, dendrogram="row", scale="row", 
#          col=myheatcol, density.info="none", trace="none", 
#          RowSideColors= myClusterSideBar)
# naive plot
genexp.dendro <- as.dendrogram(genexphc)
dendro.plot <- ggdendrogram(data = genexp.dendro, rotate = TRUE)  + 
  theme(axis.text.y = element_text(size=4))
print(dendro.plot)

gsetmp <- ScaleData(gse145281, features = mygenes)
DoHeatmap(gsetmp, features=mygenes[genexphc$order], cells = mycells,
          group.bar=T, group.by = 'patient', slot ="scale.data") + 
  theme(axis.text.y = element_text(size=2))



VlnPlot(
  object = gse145281, features = deseq.dt$gene[5005],
  group.by = "patient", idents = mycluster
)

## * PCA analysis for genes in the cell cluster
## normalized data
## mynmldata <- gse145281@assays$RNA@data[mygenes, mycells] %>% as.data.frame()

## use all the genes for PCA
## normalize the gene per cell.
mynmldata <- gse145281@assays$RNA@data[, mycells] %>% as.data.frame()
myinds <- gse145281@meta.data$patient[mycells]

## ** center mean per individual
old_colnm <- colnames(mynmldata)
colnames(mynmldata) <- myinds

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
saveRDS(s, file="SVD_Cluster2_CytotoxicTcell")
p <- 2
B <- mycentd %>%
  as.matrix() %>%
  `%*%`(., s$v[, 1:p])
B <- B[mygenes, ]

## * summarize data for stan.
modelnm <- "model_v4"
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

gsetmp <- ScaleData(gse145281, features = goldgenes$genes)
DoHeatmap(gsetmp, features=goldgenes$genes, cells = mycells,
          group.bar=T, group.by = 'patient', slot ="scale.data") + 
  theme(axis.text.y = element_text(size=9)) + scale_fill_viridis() 
## + scale_fill_gradientn(colors = c("blue", "white", "red"))

## ** get counts matrix and design matrix
Xcg <- t(as.matrix(scdata[mygenes, mycells]))
IXcg <- Xcg
Xgc <- as.matrix(scdata[mygenes, mycells])
S <- sc_tcpc[mycells]
## merge data: merge(x_cg, indhay, by="patient")
ic <- factor(scind[mycells], levels = inds)
XInd <- as.matrix(one_hot(data.table(ic = ic)))
IXInd <- as.numeric(ic)

di <- factor(scresp[mycells], levels=c(0,1))
XCond <- as.matrix(one_hot(data.table(di = di)))
IXCond <- as.numeric(di)

## ** set constants
N <- length(mycells)
K <- ncol(XInd)
G <- topgnum
J <- ncol(XCond)
P <- p

## ** save data for cmdstan
stan_rdump(c(
  "N", "K", "J", "G", "scale", "XCond", "IXCond",
  "XInd", "IXInd", "S", "Xcg", "IXcg","Xgc", "B", "P"
),
file = paste0("./", paste(c(modelnm, topgnum),collapse = "-") ,".rdump")
)
