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

## * modify symsim functions
## load from a local file
library(phytools)
library(MASS)

mysymsim <- modules::import("simulation_functions")
## modules::reload(mysymsim)

## * help functions
phyla <- function(w = rep(0.1, 6), plot_tree = F) {
  ## two populations (two conditiosn), while each has some subpopulations
  ## used when individual effects as subpopulations in phylogenetic tree
  ## w: vector of branch length in phylo tree, such as (0.1,0.2,0.3, 0.2,0.1,0.2)
  ##   - one minus the branch length describes the correlation
  ##     between the parent and the tip, so should be less than 1
  ##   - 1 to length(w)/2 describes the subpopulation structure in one case,
  ##     rest the other

  n <- ceiling(length(w) / 2)
  sub1 <- 1:n
  sub2 <- (n+1):length(w)
  ## generate newick tree format
  newick_tree <- paste("((", paste(sub1, w[sub1], sep =":", collapse = ","), "):1, (",
                       paste(sub2, w[sub2], sep =":", collapse = ","), "):1);")
  phyla <- ape::read.tree(text = newick_tree)
  if (plot_tree) {
    plot(phyla, show.tip.label = F, lwd = 2)
    ape::nodelabels(cex = 1)
    ape::tiplabels(cex = 2)
    ape::edgelabels(text = sprintf("%0.2f", phyla$edge.length),
                    col = "black", bg = "lightgreen", font = 1, adj = c(0.5, 1.5))
  }
  phyla$pop <- list(sub1 = sub1, sub2 = sub2)
  return(invisible(phyla))
}


## * class for simulator
IndividualEffectSimulatorV2 <- R6::R6Class(
  classname = "IndividualEffectSimulatorV",
  public = list(
    ## parameters f
    scale_alpha = NULL,
    qtl = NULL,
    betas = NULL,
    ### num of cells in one individual
    ncell = NULL,
    ngene = NULL,
    ### num of individuals in one condition
    nind = NULL,
    symsim_true = NULL,
    symsim_umi = NULL,
    dea = NULL,
    ## symsim parameters
    seed = NULL,
    bimod = NULL,
    capt_alpha = NULL,
    nevf = NULL,
    n_de_evf = NULL,
    sigma = NULL,
    vary = NULL,
    initialize = function(bimod = 0,
                          capt_alpha = 0.2,
                          ncell = 100,
                          ngene = 100,
                          nind = 3,
                          betas = c(-1, 0, 1),
                          ## parameters should take as default
                          scale_alpha = 0.3,
                          qtl = 0.975,
                          nevf = 10,
                          n_de_evf = 7,
                          vary = "s",
                          sigma = 0.2,
                          seed = 1L) {
      ## init simulator
      ## - bimod: bimodility in expressions
      ##   - In SymSim, be default, half of the genes will use this feature
      ##   - Though we can use continuous value between [0, 1], in SymSim examples,
      ##     they only use 0 or 1.
      ## - capt_alpha: capture efficiency mean
      ##   - used to generate the observed counts, default is 0.1
      ##   - we could vary this from (0.0, 0.2] by 0.05 step size.
      ## - sigma: population variability
      ##   - It determines gene expression variations whithn one population.
      ##   - default as 0.2 or 0.6
      self$scale_alpha <- scale_alpha
      self$qtl <- qtl
      self$ncell <- ncell
      self$ngene <- ngene
      self$betas <- betas
      self$nind <- nind

      ## symsim parameters
      self$bimod <- bimod
      self$capt_alpha <- capt_alpha
      self$nevf <- nevf
      self$n_de_evf <- 7
      self$vary <- "s"
      self$sigma <- sigma
      self$seed <- seed
    }
  ) ## end of public methods
) ## end of class IndividualEffectSimulatorV2

IndividualEffectSimulatorV2$set("public", "phyla2", function() {
  ## creating two cell populations corresponding to different conditions
  ## this is used when individual effects are added as batch effects in symsim
  t <- ape::rtree(2, tip.label = c(1,2))
  invisible(ape::compute.brlen(t, 1))
})

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
    Sigma = self$sigma,
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
  ## condition: either 1 or 2
  self$dea$cond <- symsim_umi$cell_meta$pop
  ## assign individuals for two conditions
  ## length of total number of cells
  self$dea$ind <- rep(0, self$ncell * self$nind * 2)
  self$dea$ind[self$dea$cond == 1] <- rep(seq_len(self$nind), each = self$ncell)
  self$dea$ind[self$dea$cond == 2] <- rep(seq_len(self$nind), each = self$ncell) + self$nind
  ## condition of each individual, length of self$nind * 2
  self$dea$cond_of_ind <- c(rep(1, self$nind), rep(2, self$nind))

  ## save dea field to symsim_umi
  self$symsim_umi$dea <- self$dea
})

IndividualEffectSimulatorV2$set("private", "add_indeff", function() {
  ## add individual-effect to different genes and cells.
  log2_fc <- self$dea$logFC_theoretical
})

IndividualEffectSimulatorV2$set("public", "simu", function() {
  self$simu_true()
  self$simu_umi()
  self$de()
  self$add_indeff()
})


## * main
simulator <- IndividualEffectSimulatorV2(
  
)
