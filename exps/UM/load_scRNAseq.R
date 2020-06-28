
library(Seurat)
library(here)
library(tidyverse)

# * test seurat to load the scRNAseq dataset.

data_dir <- here("data", "UM", "GSE139829_RAW", "GSM4147091")
data <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
