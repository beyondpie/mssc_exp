library(TCGAbiolinks)
library(org.Hs.eg.db)
library(data.table)
library(here)
library(tidyverse)
import::from(limma, contrasts.fit)

cancer_project <- "TCGA-UVM"

query <- TCGAbiolinks::GDCquery(
  project = cancer_project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)
sample_barcodes <- TCGAbiolinks::getResults(query, cols = c("cases"))

## same as sample_barcodes
data_tp <- TCGAbiolinks::TCGAquery_SampleTypes(
  barcode = sample_barcodes,
  typesample = "TP"
)
## below is empty
data_nt <- TCGAbiolinks::TCGAquery_SampleTypes(
  barcode = sample_barcodes,
  typesample = "NT"
)

## * download data
GDCdownload(query = query, directory = here("data", "UM", "GDCdata"))

## TODO: data is downloaded at current dir.
data_prep <- TCGAbiolinks::GDCprepare(
  query = query,
  save = T, save.filename = "TCGA_UVM_HTSeq_Countds.rds",
  directory = here("data", "UM", "GDCdata")
)

data_prep <- TCGAbiolinks::TCGAanalyze_Preprocessing(
  object = data_prep,
  cor.cut = 0.6
)
data_norm <- TCGAbiolinks::TCGAanalyze_Normalization(
  tabDF = data_prep,
  geneInfo = TCGAbiolinks::geneInfoHT,
  method = "geneLength"
)

data_fit <- TCGAbiolinks::TCGAanalyze_Filtering(
  tabDF = data_norm,
  method = "quantile",
  qnt.cut = 0.5
  )

## ** save considered genes.

ensembl2symbol_bulk <- AnnotationDbi::select(org.Hs.eg.db,
                                             keys=rownames(data_fit),
                                             column = "SYMBOL",
                                             keytype = "ENSEMBL")
saveRDS(object = ensembl2symbol_bulk, file = here("data", "UM",
                                                  "tcga_bulk_ensembl2symbol.rds"))

## * load meta data
## TODO: use library of SummarizedExperiment
## colData to get the sampel information
meta_data <- data.table::fread(here(
  "data", "UM",
  "clinical_PANCAN_patient_with_followup.tsv"
))

simple_barcodes <- substring(colnames(data_fit), 1, 12)
uvm_meta_data <- meta_data[match(simple_barcodes, bcr_patient_barcode),
                           "gender"]

uvm_degs <- TCGAbiolinks::TCGAanalyze_DEA(
  mat1 = data_fit[, which(uvm_meta_data == "FEMALE")],
  mat2 = data_fit[, which(uvm_meta_data == "MALE")],
  pipeline = "edgeR",
  ## pipeline = "limma",
  Cond1type = "FEMALE",
  Cond2type = "MALE",
  fdr.cut = 0.05,
  logFC.cut = 1,
  method = "glmLRT"
)

ensembl2symbol <- AnnotationDbi::select(org.Hs.eg.db,
  keys = rownames(uvm_degs),
  column = "SYMBOL", keytype = "ENSEMBL"
)

uvm_degs$genesymbol <- ensembl2symbol[match(rownames(uvm_degs),
                                              ensembl2symbol$ENSEMBL), ]$SYMBOL
saveRDS(object = uvm_degs, file = here("exps", "UM", "tcga_uvm_sex_deg.rds"))
