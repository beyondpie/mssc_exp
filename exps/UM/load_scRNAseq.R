library(here)
library(tidyverse)
library(Seurat)
library(harmony)

## * test seurat to load the scRNAseq dataset.
## deprecated below.
data_dir <- here("data", "UM", "GSE139829_RAW", "GSM4147091")
data <- Read10X(data.dir = data_dir)
seurat_object <- CreateSeuratObject(counts = data$`Gene Expression`)

## the result is the same as in scanpy.

## * use the UM data from Liu's lab
liu_uvm <- readRDS(here("data", "UM", "UVM_GSE139829_res.rds"))
luvm_seurat <- liu_uvm$RNA

liu_uvm_patients <- luvm_seurat@meta.data$patient

liu_uvm_celltypes <- luvm_seurat@meta.data$assign.curated
liu_uvm_celltypes <- luvm_seurat@meta.data$assign.level1_anno
liu_uvm_celltypes <- luvm_seurat@meta.data$assign.level3_anno

DimPlot(object = luvm_seurat, group = "patient", label = TRUE, pt.size = 0.1)

DimPlot(
  object = luvm_seurat, group = "assign.level3_anno",
  label = TRUE, pt.size = 0.1
)

## * batch correction
luvm_seurat <- Seurat::NormalizeData(object=luvm_seurat,scale.factor=1e4,
                                     verbose=FALSE) %>%
  FindVariableFeatures(selection.method="vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = luvm_seurat@assays$RNA@var.features, npcs=20, verbose=FALSE)


luvm_seurat <- luvm_seurat %>%
  RunHarmony("patient") %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

DimPlot(object = luvm_seurat, reduction = "umap", group = "patient",
        label = TRUE, pt.size = 0.1)

DimPlot(object = luvm_seurat, reduction = "umap", group = "assign.level3_anno",
        label = TRUE, pt.size = 0.1)
