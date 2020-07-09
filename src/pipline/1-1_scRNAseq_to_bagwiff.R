options(warn = -1)
library(here)
suppressPackageStartupMessages(library(tidyverse))
library(Seurat)
library(optparse)

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
## [9232,79105]
cnt <- floor(luvm_seurat@assays$RNA@counts)
message("Raw data:")
myt$print_sc(nrow(cnt), ncol(cnt), row="gene")
## ** filtering scRNAseq
## *** remove ERCC related pseudo genes
erccs <- grep(pattern = "^ERCC-", x = rownames(cnt), value = TRUE)
ne <- length(erccs)
if (ne > 0) {
  message(str_glue("num of ERCC: {ne}"))
  cnt <- cnt[-ercc, ]
}

## *** remove cells with higher mitochodia-related genes
mts <- grep(pattern = "^MT-", x = rownames(cnt), value = TRUE)
nmts <- length(mts)
if (nmts > 0) {
  message(str_glue("num of MT genes: {nmts}"))
  mtratios <- colSums(cnt[mts, ]) / colSums(cnt)
  ## num of high ratio of mitochotria reads
  nhmt <- length(which(ntratios > 0.2))
  if (num_high_mtratios > 0) {
    message(str_glue("num of cells with at least 20% mt reads: {nhmt}"))
  }
  cnt <- cnt[, ntratios <= 0.2]
}

## *** remove cells with too low/high-reads and genes seldomly expressed
ncnts <- colSums(cnt)
ncells <- rowSums(cnt)

low_ncnts <- quantile(ncnts, 0.05)
nlcnt <- length(which(ncnts < low_ncnts))
message(str_glue("cells with reads below {low_ncnts}: {nlcnt}"))

up_ncnts <- quantile(ncnts, 0.95)
nucnt <- length(which(ncnts > up_ncnts))
message(str_glue("cells with reads above {up_ncnts}: {nucnt}"))

low_ncells <- 500
## number of filtered genes
nfgenes <- length(which(ncells < low_ncells))
message(str_glue("genes expressed less than {low_ncells}:{nfgenes}"))

cnt <- cnt[ncells > 500, (ncounts)]
message("After filtering: ")
myt$print_sc(nrow(cnt), ncol(cnt), row = "gene")

## ** load patient gender information from GEO.
genders <- read.csv(here(mydatadir, mysubdir, args$condf),
  stringsAsFactors = FALSE, header = FALSE,
  row.names = 1,
  col.names = c("pid", "gender")
)

## * transform data for gene differential expressed analysis
## ** load the genes considered in TCGA bulkRNAseq.
ensembl2symbol_bulk <- readRDS(here(
  mydatadir, mysubdir,
  args$genef
))
## ** to bagwiff model
## bagwiff: modeling batch effects on gene-wise level
the_cell <- args$celltype

gsymbols <- rownames(luvm_seurat)
cellanno <- luvm_seurat@meta.data$assign.level3_anno

Xcg <- t(as.matrix(cnt[
  gsymbols %in% ensembl2symbol_bulk$SYMBOL,
  which(cellanno == the_cell)
]))

## label patient ids for x_cg
patients <- gsub("_.*", "", rownames(Xcg))
patient_genders <- genders[patients, 1]

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
