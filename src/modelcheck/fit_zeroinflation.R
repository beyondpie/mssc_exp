options(error = traceback)
options(warn = -1)
library(optparse)
suppressPackageStartupMessages(library(tidyverse))
library(Seurat)
import::from(here, here)
import::from(stringr, str_glue)
suppressPackageStartupMessages(library(ggpubr))
library(rtan)
options(mc.cores = parallel::detectCores())
rstan_options

## * configs
datadir <- here("data")
pbmc_IL8_dirnm <- "antiPDL1_PBMC_IL8"

## * PBMC data
## ** load seurat data
pbmcseurat <- readRDS(paste(datadir, pbmc_IL8_dirnm, "seurat.RDS", sep = "/"))
pbmccnt <- as.matrix(pbmcseurat@assays$RNA@counts)
pbmcinds <- pbmcseurat@meta.data$patient
pbmc_cellanno <- pbmcseurat@meta.data$seurat_clusters

## cytototic T cells
mycluster <- 2

ind <- "R1"

## 314 cells
oneindcells <- (pbmc_cellanno == mycluster) & (grepl(ind, colnames(pbmccnt)))
## 3885 cells
mulindcells <- pbmc_cellanno == mycluster

## ** fit log-normal poisson for each gene


## * scRNAseq benchmark
