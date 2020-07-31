options(error = traceback)
options(warn = -1)
library(tidyverse)
library(splatter)
library(scater)
library(ggplot2)
library(wesanderson)
import::from(here, here)
options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * configs
myseed <- 0L
ncelltype <- 3
nbatch <- 10
ncellbatches <- rep(200, nbatch)
group_prob <- c(0.3, 0.3, 0.4)
## group_prob <- 1.0
ngene <- 300
nmodule <- 5

## * test
set.seed(1)
sce <- scater::mockSCE()
params <- splatter::splatEstimate(sce)
sim <- splatSimulate(params)

## * splat
params <- splatter::newSplatParams()
slotNames(params)
## or use newSplatParams(nGenes=200, nCells=100) to set params
params <- splatter::setParams(params,
                             update = list(nGenes = ngene,
                                           batchCells = ncellbatches,
                                           batch.facLoc = rep(0.01, nbatch),
                                           batch.facScale = 0.01,
                                           group.prob = group_prob,
                                           mean.shape = 0.5,
                                           de.prob = c(0.2,0.2,0.3),
                                           de.downProb = rep(0.5, 3),
                                           de.facLoc = c(0.3, 0.2, 0.01),
                                           de.facScale = c(0.2, 0.5, 0.4),
                                           seed = myseed
                                           ))
sim <- splatter::splatSimulate(params, method="groups")

## * scater check
sim <- scater::logNormCounts(sim)
sim <- scater::runTSNE(sim)
sim <- scater::runPCA(sim)
p <- scater::multiplot(plotlist = list(
    scater::plotTSNE(sim, colour_by = "Batch") + scale_fill_brewer(palette = "Paired"),
    scater::plotTSNE(sim, colour_by = "Group") + scale_fill_brewer(palette="Set1")
), cols = 2)

# plotPCA(sim, shape_by = "Batch", colour_by = "Group")
