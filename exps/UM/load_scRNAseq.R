library(here)
library(tidyverse)
library(Seurat)
library(harmony)

## * use the UM data from Liu's lab
liu_uvm <- readRDS(here("data", "UM", "UVM_GSE139829_res.rds"))
luvm_seurat <- liu_uvm$RNA
liu_uvm_patients <- luvm_seurat@meta.data$patient

## ** load patient gender information from GEO.
genders <- read.csv(here("data", "UM", "genders.csv"),
                    stringsAsFactors = FALSE, header = FALSE,
                    row.names = 1,
                    col.names=c("pid","gender"))

## * batch correction
luvm_seurat <- Seurat::NormalizeData(
  object = luvm_seurat, scale.factor = 1e4,
  verbose = FALSE
) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = luvm_seurat@assays$RNA@var.features, npcs = 20, verbose = FALSE)


luvm_seurat <- luvm_seurat %>% RunHarmony("patient")
luvm_seurat <- RunUMAP(object = luvm_seurat, reduction = "harmony", dims = 1:20)

luvm_seurat <- FindNeighbors(object=luvm_seurat,
                             reduction = "harmony", dims = 1:20)
luvm_seurat <- FindClusters(object = luvm_seurat,
                            resolution = 0.5) %>% identity()

saveRDS(luvm_seurat, here("data", "UM", "UVM_GSE139829_harmony.rds"))

DimPlot(
  object = luvm_seurat, reduction = "umap", group = "patient",
  label = TRUE, pt.size = 0.1
)

p <- DimPlot(
  object = luvm_seurat, reduction = "umap", group = "assign.level3_anno",
  label = TRUE, pt.size = 0.1
)

ggsave(filename = "umap_after_harmony.pdf", path = here("exps", "UM", "pdf"),
       plot=p, width = 10, height=10)

## * re-annotate the cell clusters
## currently, we directly use the results
# from liu_uvm since the results seems to be resaonable.


