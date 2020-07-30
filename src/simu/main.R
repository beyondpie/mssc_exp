options(error = traceback)
options(warn = -1)
library(tidyverse)
library(splatter)
library(scater)
import::from(here, here)
options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * configs
myseed <- 0L
ncelltype <- 1
nbatch <- 10
ncellbatches <- rep(100, nbatch)
group_prob <- c(0.5, 0.5)
ngene <- 300
ndiffgene <- 50
nmodule <- 5

## * test
set.seed(1)
sce <- scater::mockSCE()
params <- splatter::splatEstimate(sce)
sim <- splatSimulate(params)


## * splat
params <- splatter::newSplatParams()
slotNames(params)
params <- splatter::setParams(params,
                             update = list(nGenes = ngene,
                                           batchCells = ncellbatches,
                                           group.prob = group_prob,
                                           seed = myseed
                                           ))
sim <- splatSimulate(params)
