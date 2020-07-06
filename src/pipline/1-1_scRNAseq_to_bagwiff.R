library(here)
library(tidyverse)
library(Seurat)
library(optparse)

option_list <- list(
  make_option(c("data_dir"),
    action = "store",
    type = "character",
    default = "data"
  ),
  make_option(
    c("tcga_tumor"),
    action = "store",
    type = "character",
    default = "UM"
  ),
  make_option(
    c("sc_file"),
    action = "store",
    type = "character",
    default = "test.rds"
  ),
  make_option(
    c("condf"),
    action = "store",
    type = "character",
    default = "gender.csv"
  ),
  make_option(
    c("celltype"),
    action = "store",
    type = "character",
    default = "malignant"
  ),
  make_option(
    c("genef"),
    action = "store",
    type = "character",
    default = "mygenes.rds"
  ),
  make_option(
    c("output"),
    action = "store",
    type = "character",
    default = "out.rds"
  )
)

args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()


## * re-annotate the cell clusters
## ** reload data.
luvm_seurat <- readRDS(here(args$data, args$tcga_tumor, args$sc_file))

## ** load patient gender information from GEO.
genders <- read.csv(here(args$data, args$tcga_tumor, args$condf),
  stringsAsFactors = FALSE, header = FALSE,
  row.names = 1,
  col.names = c("pid", "gender")
)

## currently, we directly use the results
# from lab since the results seems to be resaonable.

## * transform data for gene differential expressed analysis
## ** load the genes considered in TCGA bulkRNAseq.
ensembl2symbol_bulk <- readRDS(here(
  args.data, args.tcga_tumor,
  args.genef
))
## ** to bagwiff model
## bagwiff: modeling batch effects on gene-wise level
the_cell <- args.celltype

gsymbols <- rownames(luvm_seurat)
cellanno <- luvm_seurat@meta.data$assign.level3_anno
cnt <- luvm_seurat@assays$RNA@counts
## [9232,79105]

Xcg <- t(as.matrix(cnt[
  gsymbols %in% ensembl2symbol_bulk$SYMBOL,
  which(cellanno == the_cell)
]))

## label patient ids for x_cg
patients <- gsub("_.*", "", colnames(luvm_seurat))
patient_genders <- genders[patients, 1]

XInd <- as.data.frame(myt$to_onehot_matrix(patients))
XCond <- as.data.frame(myt$to_onehot_matrix(patient_genders))

N <- nrow(XCond)
J <- ncol(XCond)
K <- ncol(XInd)
G <- ncol(Xcg)

S <- as.data.frame(rowSums(Xcg))
## ** to bagmiff mdel
## bagmiff: modeling batch effects on gene-module level
## add gene module infomration.

## a trivial one
P <- 1
B <- as.data.frame(matrix(1:G, nrow = G, ncol = P))
## * use pystan to transform.
Xcg <- as.data.frame(Xcg)
save(N, J, K, G, S, P, B, XInd, XCond, Xcg,
  file = here(args$data_dir, args$tcga_tumor, args$output)
)
