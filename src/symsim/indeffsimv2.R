## simulate individual effect (bf) v2 using SymSim

## As individual effect simulator, this file is self-sufficient,
## which means it only needs several published R packages, but no
## extra need for other script in this project (different with v1).

## first version of simulating batch effect named symsim.R
## under the same dir.

## * set R environment
suppressWarnings(suppressMessages({
  library(SymSim)
  library(tidyverse)
  library(cmdstanr)
  library(bayesplot)
  library(posterior)
  library(grid)
  library(gtable)
  library(gridExtra)
  library(ggpubr)
  library(R6)
  library(ape)
}))

## warnings/errors traceback settings
options(error = traceback)
options(warn = -1)
options(mc.cores = 3)

## load help functions
options("import.path" = list(here::here("rutils"),
                             here::here("src", "mssc")))
mypseudo <- modules::import("pseudobulk")
mysymsim <- modules::import("mysymsim")



IndividualEffectSimulatorV2 <- R6::R6Class(
  classname = "IndividualEffectSimulatorV",
  public = list(
    seed = NULL,
    alpha = NULL,
    qtl = NULL,
    ## num of cells in one individual
    ncell = NULL,
    ngene = NULL,
    ## num of individuals in one condition
    nind = NULL,
    betas = NULL,
    symsim_true = NULL,
    symsim_umi = NULL,
    dea = NULL,
    ## symsim parameters
    nevf = 10,
    n_de_evf = 7,
    Sigma = 0.2,
    vary = "s",
    initialize = function(alpha = 0.3,
                          qtl = 0.975,
                          ncell = 100,
                          ngene = 100,
                          betas = c(-1, 0, 1),
                          seed = 1L) {
      self$alpha <- alpha
      self$qtl <- qtl
      self$ncell <- ncell
      self$ngene <- ngene
      self$betas <- betas
      self$nind <- length(self$betas)
      self$seed <- seed
    }
  ) ## end of public methods
) ## end of class IndividualEffectSimulatorV2

IndividualEffectSimulatorV2$set("private", "simu_true", function() {
  ## Simulate mRNA count in cells
  phylotree <- ape::rtree(2, tip.label = c(1, 2))
  phyla2 <- ape::compute.brlen(phylotree, 1)
  totalcells <- self$ncell * self$nind * 2
  min_popsize <- floor(totalcells / 2)
  
  self$symsim_true <- SymSim::SimulateTrueCounts(
    randseed = self$seed,
    ncells_total = totalcells,
    ngene = self$ngene,
    Sigma = self$Sigma,
    vary = self$vary,
    nevf =  self$nevf,
    n_de_evf = self$n_de_evf,
    phyla = phyla2,
    min_popsize = min_popsize,
    ## fixed parameters
    prop_hge = 0.0,
    gene_module_prop = 0.0,
    i_minpop = 1,
    evf_type = "discrete"
  )
})

IndividualEffectSimulatorV2$set("private", "simu_umi", function() {
  ## Simualte UMI-based scRNA-seq data based on the true count.
  data(gene_len_pool, package = "SymSim")
  gene_len <- sample(gene_len_pool, self$ngene, replace = FALSE)
  ## UMI settings
  depth_mean_umi <- 45000
  depth_sd_umi <- 4500
  alpha_mean_umi <- 0.1
  self$symsim_umi <- SymSim::True2ObservedCounts(
    true_counts = self$symsim_true$counts,
    meta_cell = self$symsim_true$cell_meta,
    protocol = "UMI",
    alpha_mean = alpha_mean_umi,
    alpha_sd = 0.02,
    gene_len = gene_len,
    depth_mean = depth_mean_umi,
    depth_sd = depth_sd_umi,
    nPCR1 = 14
  )
  ## make symsim_umi record everything we need for saving
  self$symsim_umi$true_cnt <- self$symsim_true$counts
})

IndividualEffectSimulatorV2$set("private", "de", function() {
  ## Differential expression analysis in SymSim

  self$dea <- SymSim::getDEgenes(true_counts_res = self$symsim_true,
                                   popA = 1, popB = 2)
  ## select differentially expressed genes
  deg_index <- (self$dea$nDiffEVF >= 1) & (abs(self$dea$logFC_theoretical) >= 0.6)
  self$dea$diffg <- which(deg_index == T)
  ## select non-differentially expressed genes
  self$dea$nondiffg <- setdiff(seq_len(self$ngene), self$dea$diffg)
  ## record condition for each cell
  self$dea$cond <- symsim_umi$cell_meta$pop
  ## assign individuals for two conditions
  ## length of total number of cells
  self$dea$ind <- rep(0, self$ncell * self$nind * 2)
  self$dea$ind[self$dea$cond == 1] <- rep(seq_len(self$nind), each = self$ncell)
  self$dea$ind[self$dea$cond == 2] <- rep(seq_len(self$nind), each = self$ncell) + self$nind
  ## condition of each individual, length of self$nind * 2
  self$dea$cond_of_ind <- c(rep(1, self$nind), rep(2, self$nind))
})

IndividualEffectSimulatorV2$set("public", "simu", function() {
  self$simu_true()
  self$simu_umi()
  self$de()
})

