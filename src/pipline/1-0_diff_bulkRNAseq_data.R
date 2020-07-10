options(error = traceback)
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(import::from(limma, contrasts.fit))

## * load util functions.
options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * configs
cancer_project <- "TCGA-UVM"
data_dir <- "data"
subdir <- "UM"

genes_fnm <- "tcga_bulk_gsymbol.rds"
gene_filter_method <- "quantile"
gene_qnt_cut_meanreads <- 0.1

de_pipeline <- "edgeR"
de_method <- "glmLRT"
de_fdr_init_cut <- 0.7
de_fdr_cut <- 0.05
de_logfc_init_cut <- 0.1
de_logfc_cut <- 1
de_outfnm <- "tcga_diffexp_genes.rds"
fpde_outfnm <- "tcga_fp_diffexp_genes.rds"
tnde_outfnm <- "tcga_tn_diffexp_genes.rds"

message("bulkRNAseq configs: ")
args <- list(
  de_pipeline = de_pipeline, de_method = de_method,
  de_fdr_cut = de_fdr_cut, de_logfc_cut = de_logfc_cut,
  de_outfnm = de_outfnm,
  cancer = cancer_project, genes_fnm = genes_fnm,
  gene_filter_method = gene_filter_method,
  gene_qnt_cut_meanreads = gene_qnt_cut_meanreads
)
## print(args)
message(str(args))

query <- tcgabiolinks::gdcquery(
  project = cancer_project,
  data.category = "transcriptome profiling",
  data.type = "gene expression quantification",
  workflow.type = "htseq - counts"
)
sample_barcodes <- tcgabiolinks::getresults(query, cols = c("cases"))

## same as sample_barcodes
data_tp <- tcgabiolinks::tcgaquery_sampletypes(
  barcode = sample_barcodes,
  typesample = "tp"
)

## below is empty
## data_nt <- tcgabiolinks::tcgaquery_sampletypes(
## barcode = sample_barcodes,
## typesample = "nt"
## )

## * download data
gdcdownload(query = query, directory = here(data_dir, subdir, "gdcdata"))

data_prep <- tcgabiolinks::gdcprepare(
  query = query,
  save = t, save.filename = "tcga_htseq_countds.rds",
  directory = here(data_dir, subdir, "gdcdata")
)
## * save data
saverds(data_prep, here(data_dir, subdir, "tcga_gdcprepare.rds"))
## reload data: data_prep <- readrds(here(data_dir, subdir, "tcga_gdcprepare.rds"))
data_prep <- tcgabiolinks::tcgaanalyze_preprocessing(
  object = data_prep,
  cor.cut = 0.6
)
data_norm <- tcgabiolinks::tcgaanalyze_normalization(
  tabdf = data_prep,
  geneinfo = tcgabiolinks::geneinfoht,
  method = "genelength"
)

## ** mapping ensemble to symbol
message("raw bulkrnaseq: ")
myt$print_sc(nrow(data_norm), ncol(data_norm),
  row = "gene", plat = "bulkrnaseq"
)

ensembl2symbol_bulk <- annotationdbi::select(org.hs.eg.db,
  keys = rownames(data_norm),
  column = "symbol",
  keytype = "ensembl"
)
ensembl2symbol_bulk <- ensembl2symbol_bulk[
  !duplicated(ensembl2symbol_bulk$ensembl),
]
kept_rows <- which(!is.na(ensembl2symbol_bulk$symbol))
data_norm <- data_norm[kept_rows, ]
rownames(data_norm) <- ensembl2symbol_bulk[kept_rows, "symbol"]

message("after mapping to symbol: ")
myt$print_sc(nrow(data_norm), ncol(data_norm),
  row = "gene", plat = "bulkrnaseq"
)

## ** remove mt genes
data_norm <- myt$rm_mt(data_norm)

## ** remove low-reads genes
data_fit <- tcgabiolinks::tcgaanalyze_filtering(
  tabdf = data_norm,
  method = gene_filter_method,
  qnt.cut = gene_qnt_cut_meanreads
)

## ** save considered genes.
saverds(object = rownames(data_norm), file = here(
  data_dir, subdir,
  genes_fnm
))

## * load meta data
## todo: use library of summarizedexperiment
## coldata to get the sample information
meta_data <- data.table::fread(here(
  data_dir, "tcga_patient",
  "clinical_pancan_patient_with_followup.tsv"
))

simple_barcodes <- substring(colnames(data_fit), 1, 12)
the_meta_data <- meta_data[
  match(simple_barcodes, bcr_patient_barcode),
  "gender"
]
message(stringr::str_glue("{cancer_project} patient genders: "))
table(the_meta_data)

dea <- tcgabiolinks::tcgaanalyze_dea(
  mat1 = data_fit[, which(the_meta_data == "female")],
  mat2 = data_fit[, which(the_meta_data == "male")],
  pipeline = de_pipeline,
  cond1type = "female",
  cond2type = "male",
  fdr.cut = de_fdr_init_cut,
  logfc.cut = de_logfc_init_cut,
  method = de_method
)
dea <- dea[order(dea$PValue), ]
dea$genesymbol <- rownames(dea)
message(
  stringr::str_glue(
    "init fdr({de_fdr_init_cut}) ",
    " logfc({de_logfc_init_cut}): ",
    "num of genes({nrow(dea)})"
  )
)

degs <- dea[dea$FDR < de_fdr_cut, ]
degs <- degs[abs(degs$logFC) > de_fdr_cut, ]

message(
  stringr::str_glue(
    "fdr({de_fdr_cut}) and logfc({de_logfc_cut}): num of genes({nrow(degs)})"
  )
)
saveRDS(
  object = degs, file =
    here(data_dir, subdir, de_outfnm)
)

nondegs <- dea[dea$FDR > de_fdr_cut, ]
nondegs <- nondegs[order(nondegs$PValue), ]

fp_degs <- nondegs[seq_len(2 * nrow(degs)), ]

tn_from <- nrow(degs) + 100
tn_to <- tn_from + 2 * nrow(degs)
tn_degs <- nondegs[seq(tn_from, tn_to), ]

saveRDS(
  object = fp_degs, file =
    here(data_dir, subdir, fpde_outfnm)
)

saveRDS(
  object = tn_degs, file =
    here(data_dir, subdir, tnde_outfnm)
)
