## Following Jun's suggestions
## Eval the sampler firstly by generating from the prior.

library(cmdstanr)
import::from(here, here)
suppressPackageStartupMessages(library(tidyverse))

## local modules
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## warnings/errors traceback settings
options(error = traceback)
options(warn = 0)
options(mc.cores = 3)

## set cmdstan path
set_cmdstan_path(path = paste(Sys.getenv("HOME"),
  "softwares",
  "cmdstan-2.23.0", sep = "/"))

## * load pbmc for parameter estimate and setting.
## classical genes as DE
## SNHG16, OASL, NAMPT, NFKB1, BCL2L11, TRAF4, ICAM1,
## XCL2, XCL1, CCL3L3, CCL3L1

## strong individual effect genes
## HBA1, HBA2, HBD

## possible non-DE
## TOX, YIPF5, CCL3, KDM6A, HDDC2

cell_type <- "Naive CD4+ T"
sample_cellnum <- 280
# a DE gene: DE one gene
d1g <- "NFKB1"

## the whole dataset
pbmc_seurat <- mypbmc$load_pbmc_seurat() %>%
  mypbmc$extract_from_seurat(pbmc_seurat = .)
## limit to the cell type
subscdata <- mypbmc$get_celltype_specific_scdata(pbmc_seurat$cnt,
  pbmc_seurat$resp,
  pbmc_seurat$inds,
  pbmc_seurat$ct,
  cell_type)

## ** sample cell numbers
sample_cells <- mypbmc$sample_cells_per_ind(subscdata$inds,
                                            sample_cellnum)
## note individual order is changed according to sample_cells
cnt <- subscdata$cnt[, sample_cells]
inds <- subscdata$inds[sample_cells]
resp <- subscdata$resp[sample_cells]
colsumcnt <- colSums(cnt)

## ** limit to one DE gene
## should we remove outliers firstly?
d1g_cnt <- cnt[d1g]
par(mfrow=c(1,2))
hist(d1g_cnt[resp==1])
hist(d1g_cnt[resp==0])

## fit NB dist


