library(GEOquery)
library(here)
library(tidyverse)

## * load GEO data.
gsenm <- "GSE139829"
gse <- getGEO(gsenm)

show(gse)

pdata <- pData(phenoData(gse[[1]]))
pdata[["Sex:ch1"]]
show(pdata)
