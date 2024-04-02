options(error = traceback)
options(warn = 0)
library(optparse)
suppressPackageStartupMessages(library(tidyverse))
library(Seurat)
import::from(here, here)
import::from(stringr, str_glue)
suppressPackageStartupMessages(library(ggpubr))

## fitdist in MASS
## It supports Poisson, NB, Gamma and so on.
## NB and Gamma use optim to estimate
## we can use fitdistrplus instead of MASS
## the latter is an extension towards MASS::fitdist
library(MASS)
## library(fitdistrplus)

## Support numerical way to solve MLE of Poisson lognormal
## distribution:
## - Y ~ Poisson( S * lambda)
## - log(lambda) ~ Normal(mu, sigma)
## The MLE process needs all the counts share the same
## S and lambda
## library(poilog)
##
## sads does MLE of Poisson lognormal by poilog
## sads mainly uses mle2 from bbmle package, which
## then based on mle in stats4 for MLE.
## qunminorm-paper also uses sads to estimate poisson lognormal
library(bbmle)
library(sads)

## MLE of zero-inflated and hurdle models for count data.
## library(pscl)


## This lib is to discover the zero-inflated genes
## by MLE of poisson, negive binomial.
## library(HIPPO)

## This lib is to UMI-based scRNA-Seq DEE analysis
## with negative binomial with independent dispersions
## library(NBID)

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
# modules::reload(myfit)

## * configs
datadir <- here("data")
figdir <- here("src", "modelcheck", "figures")
pbmc_IL8_dirnm <- "antiPDL1_PBMC_IL8"

## * functions
is_outlier <- function(x, up_prob = 0.995) {
  return(x > 5 * quantile(x, up_prob))
}

## * PBMC data
## ** fit log-normal poisson for each gene

## TODO: double check those high zero inflation (near or equal to one)
##       if there are some bugs in the codes.
## TODO: fit species diversity for different cell populations.

pbmcseurat <- readRDS(paste(datadir, pbmc_IL8_dirnm, "seurat.RDS", sep = "/"))
pbmccnt <- as.matrix(pbmcseurat@assays$RNA@counts)
pbmctpm <- as.matrix(pbmcseurat@assays$RNA@data)
pbmcinds <- pbmcseurat@meta.data$patient
pbmc_cellanno <- pbmcseurat@meta.data$seurat_clusters

## cytototic T cells
mycluster <- 2
ind <- "R1"

## 314 cells
oneindcells <- (pbmc_cellanno == mycluster) & (grepl(ind, colnames(pbmccnt)))
## 3885 cells
mulindcells <- pbmc_cellanno == mycluster

DEGs <- c("CCL4L1", "CCL4L2", "CCL3L1", "CCL3L3")
heavyzeroGs <- c("MIR155HG", "TNFRSF4", "ICAM1", "NA.499", "HIST2H2AA4")
heavyindeffectGs <- c("HBB", "HBA2", "HBA1")
genes <- c(DEGs, heavyzeroGs, heavyindeffectGs)

zrs_c2_R1 <- myfit$estimate_zeroratios(pbmccnt, pbmcinds, pbmc_cellanno,
  genes, whichind = "R1",
  whichcluster = 2)
rownames(zrs_c2_R1) <- genes
saveRDS(object = zrs_c2_R1,
  file = here("src", "modelcheck", "zeroratio_R1_cluster2.RDS"))

## should remove genes when all the counts are zeros
## other wise nb fitting, poislog fitting might be errors.
zrs_c1_R1 <- myfit$estimate_zeroratios(pbmccnt, pbmcinds, pbmc_cellanno,
  genes, whichind = "R1",
  whichcluster = 1)
rownames(zrs_c1_R1) <- genes
saveRDS(object = zrs_c1_R1,
  file = here("src", "modelcheck", "zeroratio_R1_cluster1.RDS"))

zrs_c2_Rall <- myfit$estimate_zeroratios(pbmccnt, pbmcinds, pbmc_cellanno,
  genes,
  whichcluster = 2)
rownames(zrs_c2_Rall) <- genes
saveRDS(object = zrs_c2_Rall,
  file = here("src", "modelcheck", "zeroratio_Rall_cluster2.RDS"))

zrs_c1_c2_R1 <- myfit$estimate_zeroratios(pbmccnt, pbmcinds, pbmc_cellanno,
  genes,
  whichcluster = c(1,2),
  whichind = "R1")
rownames(zrs_c1_c2_R1) <- genes
saveRDS(object = zrs_c1_c2_R1,
  file = here("src", "modelcheck", "zeroratio_R1_cluster1_cluster2.RDS"))


zrs_c1_c2_Rall <- myfit$estimate_zeroratios(pbmccnt, pbmcinds, pbmc_cellanno,
  genes,
  whichcluster = c(1, 2))
rownames(zrs_c1_c2_Rall) <- genes
saveRDS(object = zrs_c1_c2_Rall,
  file = here("src", "modelcheck", "zeroratio_Rall_cluster1_cluster2.RDS"))

## put them aside.
## zero-inflated and hurdle model

## ignore cell sequence depth scaling factor design for
## each individual or each cell types.

## * check why our model cannot perform as good as pseudobulk
## SymSim
