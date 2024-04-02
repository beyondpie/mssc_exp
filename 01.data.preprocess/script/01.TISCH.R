library(tidyverse)

projroot <- here::here()
TISCHdir <- file.path(projroot, "data", "TISCH")
files <- list.files(path = TISCHdir, pattern = "\\w+.rds")

a <- readRDS(file.path(TISCHdir, files[1]))
