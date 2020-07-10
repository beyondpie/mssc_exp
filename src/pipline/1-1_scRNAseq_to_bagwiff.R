options(error = traceback)
options(warn = -1)
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
library(Seurat)
library(optparse)
suppressPackageStartupMessages(library(Matrix))
import::from(stringr, str_glue)

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
    default = "UVM_GSE139829_harmony.rds"
  ),
  make_option(
    c("--condf"),
    action = "store",
    type = "character",
    default = "genders.csv"
  ),
  make_option(
    c("--celltype"),
    action = "store",
    type = "character",
    default = "Malignant"
  ),
  make_option(
    c("--genef"),
    action = "store",
    type = "character",
    default = "tcga_bulk_gsymbol.rds"
  ),
  make_option(
    c("--output"),
    action = "store",
    type = "character",
    default = "UVM_scRNAseq.RData"
  ),
  make_option(
    c("--deg"),
    action = "store",
    type = "character",
    default = "tcga_diffep_genes.rds"
  ),
  make_option(
    c("--gfilteratio"),
    action = "store",
    type = "double",
    default = 0.1
  )
)

args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()

message("load arguments: ")
print(args)
message(str(args))

mydatadir <- args$data_dir
mysubdir <- args$sub
gfilteratio <- args$gfilteratio
## * load util functions.
options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * load data
luvm_seurat <- readRDS(here(mydatadir, mysubdir, args$sc_file))

## [9232,79105]
cnt <- floor(luvm_seurat@assays$RNA@counts)
message("Raw data:")
myt$print_sc(nrow(cnt), ncol(cnt), row = "gene")

## * filtering scRNAseq
## ** remove ERCC related pseudo genes
erccs <- grep(pattern = "^ERCC-", x = rownames(cnt), value = FALSE)
ne <- length(erccs)
if (ne > 0) {
  message(str_glue("num of ERCC: {ne}"))
  cnt <- cnt[-ercc, ]
}

## ** remove mitochodira genes
cnt <- myt$rm_mt(cnt)

## ** filtering cells

## *** remove cells with too low/high-reads
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
message("After filtering low-quality cells: ")
myt$print_sc(nrow(cnt), ncol(cnt), row = "gene")

## *** select the cell type
cellanno <- luvm_seurat@meta.data$assign.level3_anno
cellanno <- cellanno[kept_cells]
the_cell <- args$celltype
cnt <- cnt[, which(cellanno == the_cell)]
message(str_glue("After choosing cell type {the_cell}"))
myt$print_sc(nrow(cnt), ncol(cnt), row = "gene")

## ** filtering genes

## *** low-quality genes
low_ncells <- gfilteratio * ncol(cnt)
## number of filtered genes
nfgenes <- length(which(ncells < low_ncells))
message(str_glue("genes expressed in at most {gfilteratio} cells: {nfgenes}"))

## *** check differential expressed genes overlap
known_degs <- readRDS(here(mydatadir, mysubdir, args$deg))
ovlp_degs <- intersect(rownames(cnt), known_degs$genesymbol)
message(str_glue("num of known DE genes: {nrow(known_degs)}"))
message(str_glue("num of DE genes left in sc: {length(ovlp_degs)}"))

## *** load the genes considered in TCGA bulkRNAseq.
bulk_gsymbol <- readRDS(here(
  mydatadir, mysubdir,
  args$genef
))
message(str_glue("number of bulk genes: {length(bulk_gsymbol)}"))
ovlp_scbulk <- intersect(rownames(cnt), bulk_gsymbol)
message(
  str_glue("number of overlap between sc and bulk:{length(ovlp_scbulk)}")
)

## *** finally filter the genes
high_quality_genes <- rownames(cnt)[which(ncells >= low_ncells)]
cnt <- cnt[union(union(high_quality_genes, ovlp_degs), ovlp_scbulk), ]
message("After filtering genes")
myt$print_sc(nrow(cnt), ncol(cnt), row = "gene")

## ** load patient gender information from GEO.
genders <- read.csv(here(mydatadir, mysubdir, args$condf),
  stringsAsFactors = FALSE, header = FALSE,
  row.names = 1,
  col.names = c("pid", "gender")
)
message("sRNAseq individual gender summay")
print(genders)

patients <- gsub("_.*", "", colnames(cnt))
message("Cells per patients")
table(patients)

patient_genders <- genders[patients, 1]
message("cells in case and control")
table(patient_genders)


## * to bagwiff model
## bagwiff: modeling batch effects on gene-wise level
Xcg <- t(as.matrix(cnt))
XInd <- myt$to_onehot_matrix(patients)
XCond <- myt$to_onehot_matrix(patient_genders)

N <- nrow(XCond)
J <- ncol(XCond)
K <- ncol(XInd)
G <- ncol(Xcg)
S <- rowSums(Xcg)

## bagmiff mdel
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
