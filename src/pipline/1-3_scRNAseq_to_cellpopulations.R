options(error = traceback)
options(warn = -1)
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
library(Seurat)
library(optparse)
suppressPackageStartupMessages(library(Matrix))
import::from(stringr, str_glue)

options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * configs
## ** data related configs
datadir <- "data"
subdir <- "UM"

de_outfnm <- "tcga_diffexp_genes.rds"
fpde_outfnm <- "tcga_fp_diffexp_genes.rds"
tnde_outfnm <- "tcga_tn_diffexp_genes.rds"
tcga_genefnm <- "tcga_bulk_gsymbol.rds"
sc_file <- "UVM_GSE139829_harmony.rds"
sc_condfnm <- "genders.csv"

sc_out_no_malignant <- "scRNAseq_no_malignant.rds"
sc_out_all <- "scRNAseq_allcell.rds"

gfilteratio <- 0.1

## * load genes information
deg <- readRDS(here(datadir, subdir, de_outfnm))
fpdeg <- readRDS(here(datadir, subdir, fpde_outfnm))
tndeg <- readRDS(here(datadir, subdir, tnde_outfnm))

## * load data
luvm_seurat <- readRDS(here(datadir, subdir, sc_file))
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

## ** function for select different cell populations.
filtergenes_savedata <- function(mycnt, outfnm) {
    myt$print_sc(nrow(mycnt), ncol(mycnt), row = "gene")
    ## ** filtering genes
    ## *** low-quality genes
    low_ncells <- gfilteratio * ncol(mycnt)
    ## number of filtered genes
    nfgenes <- length(which(ncells < low_ncells))
    message(str_glue("genes expressed in at most {gfilteratio} cells: {nfgenes}"))

    ## *** check differential expressed genes overlap
    message("positve deg:")
    ovlp_deg <- myt$stat_geneset(rownames(mycnt), deg$genesymbol)

    message("false posive deg:")
    ovlp_fpdeg <- myt$stat_geneset(rownames(mycnt), fpdeg$genesymbol)

    message("true negative deg:")
    ovlp_tndeg <- myt$stat_geneset(rownames(mycnt), tndeg$genesymbol)

    nondeg <- union(ovlp_fpdeg, ovlp_tndeg)
    dea <- union(ovlp_deg, nondeg)

    ## *** load the genes considered in TCGA bulkRNAseq.
    bulk_gsymbol <- readRDS(here(
        datadir, subdir,
        tcga_genefnm
    ))
    message(str_glue("number of bulk genes: {length(bulk_gsymbol)}"))
    ovlp_scbulk <- intersect(rownames(mycnt), bulk_gsymbol)
    message(
        str_glue("number of overlap between sc and bulk:{length(ovlp_scbulk)}")
    )

    ## *** finally filter the genes
    high_quality_genes <- rownames(mycnt)[which(ncells >= low_ncells)]
    mycnt <- mycnt[union(union(high_quality_genes, dea), ovlp_scbulk), ]
    message("After filtering genes")
    myt$print_sc(nrow(mycnt), ncol(mycnt), row = "gene")

    ## ** load patient gender information from GEO.
    conds <- read.csv(here(datadir, subdir, sc_condfnm),
        stringsAsFactors = FALSE, header = FALSE,
        row.names = 1,
        col.names = c("pid", "gender")
    )

    message("sRNAseq conds")
    print(conds)

    batches <- gsub("_.*", "", colnames(mycnt))
    table(batches)

    conds <- conds[batches, 1]
    table(conds)

    totmycntpcell <- colSums(mycnt)

    ## * save data
    saveRDS(
        list(
            cnt = mycnt, batches = batches,
            conds = conds, totmycntpcell = totmycntpcell
        ),
        here(datadir, subdir, outfnm)
    )
}

## * select non-malignant cells
cellanno <- luvm_seurat@meta.data$assign.level3_anno
cellanno <- cellanno[kept_cells]
malignant <- "Malignant"
cnt_no_malignant <- cnt[, which(cellanno != malignant)]

filtergenes_savedata(cnt_no_malignant, sc_out_no_malignant)

## * select all the cells
filtergenes_savedata(cnt, sc_out_all)
