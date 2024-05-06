library(tidyverse)
library(Seurat)
projd <- here::here()

blcvdSeu <- readRDS(
  file.path(projd, "data", "blish_covid.seu.rds"))
blcvSeu <- UpdateSeuratObject(blcvdSeu)
saveRDS(object = blcvSeu,
  file = file.path(projd, "data", "blish_covid.seu.updated.rds"))
