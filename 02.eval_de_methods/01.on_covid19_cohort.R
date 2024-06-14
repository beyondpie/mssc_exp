library(Seurat)
library(tidyverse)

# * meta
projd <- here::here()
covidir <- file.path(projd, "data", "COVID19_large_cohort")
mt <- "CD8"
seu <- readRDS(file.path(covidir, "majortype",
  str_glue("covid19.large.{mt}.seu.rds")))

# 1. set up ground truth
sampleMeta <- file.path(covidir, "GSE158055_sample_metadata.csv") |>
  data.table::fread(file = _, sep = ",", header = TRUE, data.table = FALSE)
oldNames <- colnames(sampleMeta)
simNames <- oldNames |>
  x => gsub("characteristics: +", "", x) |>
  x => gsub("\\[|\\]|\\(|\\)", "", x) |>
  x => gsub("/|-| +", "_", x) 
colnames(sampleMeta) <- simNames

# visualize data
with(sampleMeta, quantile(as.integer(Age[SARS_CoV_2 == "negative"])))
subset(sampleMeta, Age != "unknown") |>
  with(data = _,
    quantile(as.integer(Age[SARS_CoV_2 == "positive"]), na.rm = TRUE))

# select one sample per patient
dup <- duplicated(sampleMeta$Patients)
sampleMeta <- sampleMeta[!dup, ]
control_samples <- with(sampleMeta,
  Sample_name[SARS_CoV_2 == "negative"])
positive_samples <- with(sampleMeta,
  Sample_name[SARS_CoV_2 == "positive"])

# label seu
seu@meta.data$sampleID <- levels(seu@meta.data$sampleID)[seu@meta.data$sampleID]
sseu <- subset(seu, sampleID %in% sampleMeta$Sample_name)
sseu@meta.data$case <- "positive"
sseu@meta.data$case[sseu@meta.data$sampleID %in% control_samples] <-
  "negative"

# run wilcox.test on pseudo bulk levels
pseudo_sseu <- AggregateExpression(
  object = sseu, group.by = "case",
  normalization.method = "LogNormalize",
  return.seurat = TRUE
)
bulk_de <- FindMarkers(object = pseudo_sseu,
  cells.1 = "negative",
  cells.2 = "positive",
  test.use = "wilcox",
  logfc.threshold = 0.1)



# 2. perform different methods
# * EdgeR
# * DESeq2
# * Seurat
# * limma
# * Wilcon test
# * t test
# * pseudo bulk
