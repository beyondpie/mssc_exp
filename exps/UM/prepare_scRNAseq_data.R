library(here)
library(tidyverse)
library(Seurat)
library(harmony)

options("import.path" = here("rutils"))
myt <- modules::import("transform")
## use modules::reload(myt) to reload the module if any updates.

## * use the UM data from Liu's lab
lab_uvm <- readRDS(here("data", "UM", "UVM_GSE139829_res.rds"))
luvm_seurat <- lab_uvm$RNA
lab_patients <- luvm_seurat@meta.data$patient

## TODO: single cell data filtering criteria.
## TODO: re cell type annotation after batch correction.

## * batch correction
luvm_seurat <- Seurat::NormalizeData(
  object = luvm_seurat, scale.factor = 1e4,
  verbose = FALSE
) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(
    pc.genes = luvm_seurat@assays$RNA@var.features,
    npcs = 20, verbose = FALSE
  )


luvm_seurat <- luvm_seurat %>% RunHarmony("patient")
luvm_seurat <- RunUMAP(object = luvm_seurat, reduction = "harmony", dims = 1:20)

luvm_seurat <- FindNeighbors(
  object = luvm_seurat,
  reduction = "harmony", dims = 1:20
)
luvm_seurat <- FindClusters(
  object = luvm_seurat,
  resolution = 0.5
) %>% identity()

saveRDS(luvm_seurat, here("data", "UM", "UVM_GSE139829_harmony.rds"))

## * reload data.
luvm_seurat <- readRDS(here("data", "UM", "UVM_GSE139829_harmony.rds"))

## * plot batch correction
umap_p <- DimPlot(
  object = luvm_seurat, reduction = "umap", group = "patient",
  label = TRUE, pt.size = 0.1
)

umap_cell <- DimPlot(
  object = luvm_seurat, reduction = "umap", group = "assign.level3_anno",
  label = TRUE, pt.size = 0.1
)

ggpubr::ggarrange(umap_p, umap_cell, nrow = 1, ncol = 2) %>%
  ggpubr::ggexport(
    filename = here(
      "exps", "UM", "result",
      "UMAP After Batch Correction.pdf"
    ), ncol = 2, nrow = 1,
    height = 10, width = 20
  )

## * re-annotate the cell clusters

## ** load patient gender information from GEO.
genders <- read.csv(here("data", "UM", "genders.csv"),
  stringsAsFactors = FALSE, header = FALSE,
  row.names = 1,
  col.names = c("pid", "gender")
)
## currently, we directly use the results
# from lab since the results seems to be resaonable.

## * transform data for gene differential expressed analysis
## ** load the genes considered in TCGA bulkRNAseq.
ensembl2symbol_bulk <- readRDS(here(
  "data", "UM",
  "tcga_bulk_ensembl2symbol.rds"
))
## ** to bagwiff model
## bagwiff: modeling batch effects on gene-wise level

the_cell <- "Malignant"

gsymbols <- rownames(luvm_seurat)
cellanno <- luvm_seurat@meta.data$assign.level3_anno
cnt <- luvm_seurat@assays$RNA@counts
## [9232,79105]
Xcg <- t(as.matrix(cnt[
  gsymbols %in% ensembl2symbol_bulk$SYMBOL,
  which(cellanno == the_cell)
]))

## TODO: if we need to further filter the genes in x_cg

## label patient ids for x_cg
patients <- gsub("_.*", "", colnames(luvm_seurat))
patient_genders <- genders[patients, 1]

XInd <- myt$to_onehot_matrix(patients)
XCond <- myt$to_onehot_matrix(patient_genders)

N <- nrow(XCond)
J <- ncol(XCond)
K <- ncol(XInd)
G <- ncol(Xcg)

S <- rowSums(Xcg)
## ** to bagmiff mdel
## bagmiff: modeling batch effects on gene-module level
## add gene module infomration.

## a trivial one
P <- 1
B <- matrix(1:G, nrow = G, ncol = P)

## * save for stan
## out of memory error for long vector

## myt$quickdump(name = here("data", "UM", "rstan", "scRNAseq_genewise.rdump"),
## myenv=environment())

## save(N, J, K, G, S, P, B, XInd, XCond, Xcg,
## file = here("data", "UM", "rstan", "sc_genewise.RData"))
data <- list(
  N = N, J = J, K = K, G = G, S = S, P = P, B = B,
  XInd = XInd, XCond = XCond, Xcg = Xcg
)
cmdstanr::write_stan_json(data, here("data", "UM", "rstan", "sc_genewise.json"))
