options(error = traceback)
options(warn = -1)
library(optparse)
suppressPackageStartupMessages(library(tidyverse))
library(Seurat)
import::from(here, here)
import::from(stringr, str_glue)
suppressPackageStartupMessages(library(ggpubr))


## * load bulk data
## * detect DEGs in bulk.

## * load scRNAseq data
## * analyze the DEGs in scRNAseq
