options(error = traceback)
options(warn = -1)

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(import::from(here, here))
suppressPackageStartupMessages(library(DESeq2))
import::from(stringr, str_glue)

options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * configs
## ** data related configs
data_dir <- "data"
subdir <- "UM"

de_outfnm <- "tcga_diffexp_genes.rds"
fpde_outfnm <- "tcga_fp_diffexp_genes.rds"
tnde_outfnm <- "tcga_tn_diffexp_genes.rds"
scdata_sumfnm <- "sampled_scRNAseq_summary.rds"

## ** DESeq2

## * load genes information
deg <- readRDS(here(data_dir, subdir, de_outfnm))
fpdeg <- readRDS(here(data_dir, subdir, fpde_outfnm))
tndeg <- readRDS(here(data_dir, subdir, tnde_outfnm))
sc_data_list <- readRDS(here(data_dir, subdir, scdata_sumfnm))
sc_genes <- rownames(sc_data_list$cnt)

## * retrieve de/nonde-related genes
degnms <- myt$stat_geneset(sc_genes, deg$genesymbol)
fpdegnms <- myt$stat_geneset(sc_genes, fpdeg$genesymbol)
tndegnms <- myt$stat_geneset(sc_genes, tndeg$genesymbol)

## * psuedo-bulk method
mybatches <- sc_data_list$batches
myconds <- sc_data_list$conds
names(myconds) <- mybatches

ubatches <- unique(mybatches)
uconds <- ubatches %>% myconds[.]

mypsedobulk <- ubatches %>%
    map(
        .f = function(batch) {
            rowSums(sc_data_list$cnt[, mybatches %in% batch])
        }
    ) %>%
    do.call(what = cbind, args = .)

colnames(mypsedobulk) <- ubatches

deseqds <- DESeq2::DESeqDataSetFromMatrix(
    countData = mypsedobulk,
    colData = data.frame(uconds),
    design = ~uconds
)

deseqres <- DESeq2::DESeq(deseqds) %>%
    DESeq2::results() %>%
    data.frame()

## * eval result
getndegnms <- function(myndegnm = "extreme") {
    if (myndegnm == "extreme") {
        ndegnms <- tndegnms
    }
    if (myndegnm == "nearpositive") {
        ndegnms <- fpdegnms
    }
    if (myndegnm == "all") {
        ndegnms <- c(tndegnms, fpdegnms)
    }
    message(str_glue("using {myndegnm} negatives"))
    return(ndegnms)
}

evalDESeq2 <- function(scorecol = "pvalue",
                       myndegnm = "extreme", mydegnms = degnms) {
    ndegnms <- getndegnms(myndegnm)
    mybackend <- c(rep(TRUE, length(mydegnms)), rep(FALSE, length(ndegnms)))
    bgnms <- c(mydegnms, ndegnms)
    names(mybackend) <- bgnms

    scores <- deseqres[[scorecol]]
    names(scores) <- rownames(deseqres)
    scores[is.na(scores)] <- 1.0
    myauc <- myt$fmtflt(caTools::colAUC(scores[bgnms], mybackend))
    message(str_glue("use DESeq2 results col {scorecol}"))
    message(str_glue("AUC: {myauc}"))
    return(list(auc = myauc, sts = scores))
}

for (myndeg in c("extreme", "nearpositive", "all")) {
    for (mycol in c("padj", "pvalue")) {
        tmp <- evalDESeq2(mycol, myndeg)
    }
}
