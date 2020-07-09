## options(error = traceback)
options(warn = -1)
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
library(Seurat)
library(optparse)
suppressPackageStartupMessages(library(Matrix))

option_list <- list(
  make_option(
    c("--data_dir"),
    action = "store",
    type = "character",
    default = "data"
  ),
  make_option(
    c("--sub"),
    action = "store",
    type = "character",
    default = "UM"
  ),
  make_option(
    c("--sc_file"),
    action = "store",
    type = "character",
    default = "test.rds"
  ),
  make_option(
    c("--condf"),
    action = "store",
    type = "character",
    default = "gender.csv"
  ),
  make_option(
    c("--celltype"),
    action = "store",
    type = "character",
    default = "malignant"
  ),
  make_option(
    c("--genef"),
    action = "store",
    type = "character",
    default = "mygenes.rds"
  ),
  make_option(
    c("--output"),
    action = "store",
    type = "character",
    default = "out.rds"
  )
)

args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()

mydatadir <- args$data_dir
mysubdir <- args$sub

## * load util functions.
options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * re-annotate the cell clusters
## ** reload data.
luvm_seurat <- readRDS(here(mydatadir, mysubdir, args$sc_file))
cellanno <- luvm_seurat@meta.data$assign.level3_anno

## [9232,79105]
cnt <- floor(luvm_seurat@assays$RNA@counts)
message("Raw data:")
myt$print_sc(nrow(cnt), ncol(cnt), row = "gene")
## ** filtering scRNAseq
## *** remove ERCC related pseudo genes
erccs <- grep(pattern = "^ERCC-", x = rownames(cnt), value = FALSE)
ne <- length(erccs)
if (ne > 0) {
  message(str_glue("num of ERCC: {ne}"))
  cnt <- cnt[-ercc, ]
}

## *** remove mitochodira genes
mts <- grep(pattern = "^MT-", x = rownames(cnt), value = FALSE)
nmts <- length(mts)
if (nmts > 1) {
  message(str_glue("num of MT genes: {nmts}"))
  mtratios <- colSums(cnt[mts, ]) / colSums(cnt)
  ## num of high ratio of mitochotria reads
  nhmt <- length(which(mtratios > 0.2))
  if (nhmt > 0) {
    message(str_glue("num of cells with at least 20% mt reads: {nhmt}"))
    ## cnt <- cnt[, mtratios <= 0.2]
  }
  cnt <- cnt[-mts, ]
}

## *** remove cells with too low/high-reads and genes seldomly expressed
ncnts <- colSums(cnt)
ncells <- rowSums(cnt)

low_ncnts <- quantile(ncnts, 0.025)
nlcnt <- length(which(ncnts < low_ncnts))
message(str_glue("cells with reads below 2.5% ({low_ncnts}): {nlcnt}"))

up_ncnts <- quantile(ncnts, 0.975)
nucnt <- length(which(ncnts > up_ncnts))
message(str_glue("cells with reads above 97.5% ({up_ncnts}): {nucnt}"))

kept_cells <- (ncnts >= low_ncnts) & (ncnts <= up_ncnts)
cnt <- cnt[, kept_cells]
cellanno <- cellanno[kept_cells]

low_ncells <- 0.1 * ncol(cnt)
## number of filtered genes
nfgenes <- length(which(ncells < low_ncells))
message(str_glue("genes expressed in at most 0.1 cells: {nfgenes}"))
cnt <- cnt[ncells >= low_ncells, ]

message("After filtering: ")
myt$print_sc(nrow(cnt), ncol(cnt), row = "gene")

## ** load patient gender information from GEO.
genders <- read.csv(here(mydatadir, mysubdir, args$condf),
  stringsAsFactors = FALSE, header = FALSE,
  row.names = 1,
  col.names = c("pid", "gender")
)
message("sRNAseq individual gender summay")
print(genders)

## * transform data for gene differential expressed analysis
## ** load the genes considered in TCGA bulkRNAseq.
ensembl2symbol_bulk <- readRDS(here(
  mydatadir, mysubdir,
  args$genef
))
## ** to bagwiff model
## bagwiff: modeling batch effects on gene-wise level
the_cell <- args$celltype
gsymbols <- rownames(cnt)

Xcg <- t(as.matrix(cnt[
  gsymbols %in% ensembl2symbol_bulk$SYMBOL,
  which(cellanno == the_cell)
]))

message("After choosing cell type ", the_cell, " and genes in bulk")
myt$print_sc(nrow(Xcg), ncol(Xcg), row = "cell")

## label patient ids for x_cg
patients <- gsub("_.*", "", rownames(Xcg))
message("Cells per patients")
table(patients)

patient_genders <- genders[patients, 1]
message("cells in case and control")
table(patient_genders)

## *** Subsampling the data

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
P <- 1L
B <- matrix(1:G, nrow = G, ncol = P)
## * use pystan to transform.
Xcg <- as.data.frame(Xcg)
XInd <- as.data.frame(XInd)
XCond <- as.data.frame(XCond)
save(N, J, K, G, S, P, B, XInd, XCond, Xcg,
  file = here(mydatadir, mysubdir, args$output)
)
