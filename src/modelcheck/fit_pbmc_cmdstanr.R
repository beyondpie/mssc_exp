library(cmdstanr)
set_cmdstan_path(path = paste(Sys.getenv("HOME"),
  "softwares",
  "cmdstan-2.23.0", sep = "/"))
library(bayesplot)
library(posterior)
library(Seurat)
suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)

options("import.path" = here("rutils"))
myt <- modules::import("transform")

myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")

## * warnings/errors traceback settings
options(error = traceback)
options(warn = 0)

## * ggplot settings.
myaxis_text <- ggplot2::element_text(size = 13, face = "bold",
  color = "black")
myaxis_text_x <- myaxis_text
myaxis_title_x <- ggplot2::element_blank()
myaxis_text_y <- myaxis_text
myaxis_title_y <- ggplot2::element_text(size = 15, face = "bold",
  color = "black")
mytitle <- ggplot2::element_text(size = 20, hjust = 0.5,
  color = "black", face = "bold",
  family = "Helvetica")


## * load data
datadir <- here("data")
figdir <- here("src", "modelcheck", "figures")
pbmc_IL8_dirnm <- "antiPDL1_PBMC_IL8"

pbmcseurat <- readRDS(paste(datadir, pbmc_IL8_dirnm, "seurat.RDS", sep = "/"))
pbmccnt <- as.matrix(pbmcseurat@assays$RNA@counts)
pbmcinds <- pbmcseurat@meta.data$patient
pbmc_cellanno <- pbmcseurat@meta.data$seurat_clusters
pbmc_cond <- pbmcseurat@meta.data$response
## ** select genes
DEGs <- c("CCL4L1", "CCL4L2", "CCL3L1", "CCL3L3")
heavyzeroGs <- c("MIR155HG", "TNFRSF4", "ICAM1", "NA.499", "HIST2H2AA4")
heavyindeffectGs <- c("HBB", "HBA2", "HBA1")

## * utilities
get_stan_bagwiff_data <- function(whichgenes,
                                  whichcluster = mycluster,
                                  cnts,
                                  cellmeta_inds,
                                  cellmeta_conds,
                                  cellmeta_clusters,
                                  rmoutliers = T) {
  whichcells <- cellmeta_clusters %in% whichcluster
  tmp <- cnts[whichgenes, whichcells]
  if (rmoutliers) {
    outliers <- myfit$grpoutliers(tmp)
  } else {
    outliers <- rep(F, ncol(tmp))
  }
  mycnts <- tmp[, !outliers]
  mytotcnts <- colSums(cnts[, whichcells])[!outliers]
  myinds <- cellmeta_inds[whichcells][!outliers]
  myconds <- cellmeta_conds[whichcells][!outliers]
  invisible(list(
    standata = myt$to_bagwiff_r(
      mycnts, myinds, myconds, mytotcnts),
    cnts = mycnts,
    totcnts = mytotcnts,
    inds <- myinds,
    conds <- myconds))
}


get_init_nb_phi <- function(cnts) {
  get_single_init <- function(x) {
    tryCatch({
      nbfit <- MASS::fitdistr(x, densfun = "negative binomial",
        lower = c(1e-5, 1e-5))
      invisible(nbfit$estimate["size"])
    }, error = function(cond) {
      message(cond)
      invisible(1.0)
    }
    )
  }
  vapply(seq_len(nrow(cnts)), FUN = function(i) {
    get_single_init(cnts[i, ])
  }, FUN.VALUE = 1.0)
}

get_mystanmodel <- function(fnm) {
  cmdstanr::cmdstan_model(
    stan_file = here("src", "stan", fnm),
    compile = T,
    stanc_options = list(include_paths = here("src", "stan")),
    dir = here("src", "modelcheck", "stan_binary")
  )
}

sampling_mystanmodel <- function(mystanmodel,
                                 mystandata,
                                 mydelta = 0.82,
                                 mytreedepth = 11,
                                 num_chains = 4,
                                 myinit = NULL) {
  mystanmodel$sample(data = mystandata,
    seed = 355113L,
    refresh = 200,
    iter_warmup = 500,
    iter_sampling = 500,
    adapt_delta = mydelta,
    max_treedepth = mytreedepth,
    validate_csv = T,
    show_messages = T,
    chains = num_chains,
    parallel_chains = getOption("mc.cores", num_chains),
    init = myinit,
    output_dir = here("src", "modelcheck",
      "stan_sampling")
  )
}

## ** model declare
stanmodel_poiglm_gw_mc <- get_mystanmodel("v1-1.stan")
stanmodel_nbglm_gw_mc <- get_mystanmodel("v3-1.stan")

## stanmodel_nblognm_gw_mc <- get_mystanmodel("v4-1.stan")

## * sampling
## ** set data
mygenes <- c("CCL3L3", "ICAM1", "HBB")
## cytototic T cells
mycluster <- c(2)
standata_cytoxictcell <- get_stan_bagwiff_data(mygenes, mycluster,
  pbmccnt,
  pbmcinds,
  pbmc_cond,
  pbmc_cellanno, T)

## ** sampling
stanfit_poiglm_gw_mc <- sampling_mystanmodel(
  stanmodel_poiglm_gw_mc, standata_cytoxictcell$standata,
  0.83, 20, 4)
stanfit_poiglm_gw_mc$save_object(file = here("src", "modelcheck",
  "stan_save", "v1-1.rds"))

## need to assign all the variables.
## init_nb_phis <- get_init_nb_phi(standata_cytoxictcell$cnts)

stanfit_nbglm_gw_mc <- sampling_mystanmodel(
  stanmodel_nbglm_gw_mc, standata_cytoxictcell$standata,
  0.85, 20, 4)
stanfit_nbglm_gw_mc$save_object(file = here("src", "modelcheck",
  "stan_save", "v3-1.rds"))

## stanfit_nblognm_gw_mc <- sampling_mystanmodel(
##   stanmodel_nblognm_gw_mc, standata_cytoxictcell,
##   0.85, 20, 4
## )

## ** view results
stanfit_poiglm_gw_mc <- readRDS(here("src", "modelcheck",
  "stan_save", "v1-1.rds"))
stanfit_nbglm_gw_mc <- readRDS(here("src", "modelcheck",
  "stan_save", "v3-1.rds"))

## get dataframe num_sample x num_variable
mcdf_poiglm_gw <- posterior::as_draws_df(stanfit_poiglm_gw_mc$draws())
mcdf_nbglm_gw <- posterior::as_draws_df(stanfit_nbglm_gw_mc$draws())

## *** view individual effect estimate
get_indeffs <- function(stan_mcdf, mygenes, numinds = 10) {
  get_indnms <- function(geneid, numinds) {
    ind_colnms <- vapply(seq_len(numinds),
      FUN = function(j) {
        paste0("MuInd[", geneid, ",", j, "]")
      },
      FUN.VALUE = "MuInd")
  }
  geneindex <- seq_len(length(mygenes))
  res <- lapply(geneindex,
    FUN = function(i) {
      indcols <- get_indnms(i, numinds)
      tmp <- stan_mcdf[, indcols]
      colnames(tmp) <- paste0(mygenes[i], "-", "Ind", seq_len(numinds))
      invisible(tmp)
    })
  names(res) <- mygenes
  invisible(res)
}


get_indeff_mcmcarea <- function(mygenes, indeff_df, modelnm = "PoissonGLM") {
  lapply(mygenes,
    FUN = function(g) {
      bayesplot::mcmc_areas(indeff_df[[g]],
        prob = 0.95,
        point_est = "mean",
      ) +
        theme(axis.text.x = myaxis_text,
          axis.text.y = myaxis_text,
          axis.title.y = myaxis_title_y,
          plot.title = mytitle) +
        ggtitle(str_glue("{g}_{modelnm}"))
    })
}

indeff_poiglm_gw <- get_indeffs(mcdf_poiglm_gw, mygenes)
indeff_nbglm_gw <- get_indeffs(mcdf_nbglm_gw, mygenes)

mcmcarea_indeff_poiglm_gw <- get_indeff_mcmcarea(mygenes,
  indeff_poiglm_gw,
  modelnm = "PoissonGLM")
mcmcarea_indeff_nbglm_gw <- get_indeff_mcmcarea(mygenes,
  indeff_nbglm_gw,
  modelnm = "NegaBinomialGLM")
bayesplot::bayesplot_grid(
  plots = c(mcmcarea_indeff_poiglm_gw,
    mcmcarea_indeff_nbglm_gw),
  grid_args = list(nrow = 2)
)

## *** view conditional effect estimations
get_condeffs <- function(stan_mcdf, mygenes, numconds = 2) {
  get_indnms <- function(geneid, numinds) {
    ind_colnms <- vapply(seq_len(numinds),
      FUN = function(j) {
        paste0("MuCond[", geneid, ",", j, "]")
      },
      FUN.VALUE = "MuCond")
  }
  geneindex <- seq_len(length(mygenes))
  res <- lapply(geneindex,
    FUN = function(i) {
      indcols <- get_indnms(i, numconds)
      tmp <- stan_mcdf[, indcols]
      colnames(tmp) <- paste0(mygenes[i], "-", "Cond", seq_len(numconds))
      tmp$delta <- tmp[, 1] - tmp[, 2]
      invisible(tmp)
    })
  names(res) <- mygenes
  invisible(res)
}

condeff_poiglm_gw <- get_condeffs(mcdf_poiglm_gw, mygenes)
condeff_nbglm_gw <- get_condeffs(mcdf_nbglm_gw, mygenes)

mcmcarea_condeff_poiglm_gw <- get_indeff_mcmcarea(mygenes,
  condeff_poiglm_gw,
  modelnm = "PoissonGLM")
mcmcarea_condeff_nbglm_gw <- get_indeff_mcmcarea(mygenes,
  condeff_nbglm_gw,
  modelnm = "NegaBinomialGLM")
bayesplot::bayesplot_grid(
  plots = c(mcmcarea_condeff_poiglm_gw,
    mcmcarea_condeff_nbglm_gw),
  grid_args = list(nrow = 2)
)

## *** view intercetps, i.e., the overall mean
get_intercepts <- function(stan_mcdf, mygenes) {
  geneindex <- seq_len(length(mygenes))
  res <- lapply(geneindex,
    FUN = function(i) {
      tmp <- stan_mcdf[, paste0("MuG[", i, "]"), drop = F]
      colnames(tmp) <- paste0(mygenes[i], "-", "Intercept")
      invisible(tmp)
    })
  names(res) <- mygenes
  invisible(res)
}

intercept_poiglm_gw <- get_intercepts(mcdf_poiglm_gw, mygenes)
mcmcarea_intercept_poiglm_gw <- get_indeff_mcmcarea(mygenes,
  intercept_poiglm_gw,
  modelnm = "PoissonGLM")
intercept_nbglm_gw <- get_intercepts(mcdf_nbglm_gw, mygenes)
mcmcarea_intercept_nbglm_gw <- get_indeff_mcmcarea(mygenes,
  intercept_nbglm_gw,
  modelnm = "NegaBinomialGLM")

bayesplot::bayesplot_grid(
  plots = c(mcmcarea_intercept_poiglm_gw,
    mcmcarea_intercept_nbglm_gw),
  grid_args = list(nrow = 2)
)

## *** check pseudobulk with DESeq2 analysis, and t statistics from posterior
## pseudobulk
tempcluster <- c(2)
## mygenes <- c(DEGs, heavyzeroGs, heavyindeffectGs)
mygenes <- c("CCL3L3", "ICAM1", "HBB")

tempcells <- pbmc_cellanno %in% tempcluster
tempcnts <- pbmccnt[mygenes, tempcells]
tempbatches <- paste0("b", pbmcinds[tempcells])
tempconds <- paste0("c", pbmc_cond[tempcells])

temp_psedobulk_DESeq2 <- mypseudo$pseudobulk_deseq2(tempcnts,
  tempbatches,
  tempconds)

##         baseMean log2FoldChange     lfcSE        stat    pvalue      padj
## CCL3L3  47.17908   -0.008957221 0.5693854 -0.01573139 0.9874487 0.9874487
## ICAM1   14.26437    0.592224791 0.4927516  1.20187279 0.2294128 0.3441192
## HBB    316.58898   -4.196680297 2.8890210 -1.45263057 0.1463264 0.3441192


tempcells <- pbmc_cellanno %in% tempcluster
tempcnts <- pbmccnt[, tempcells]
tempbatches <- paste0("b", pbmcinds[tempcells])
tempconds <- paste0("c", pbmc_cond[tempcells])

temp_psedobulk_DESeq2 <- mypseudo$pseudobulk_deseq2(tempcnts,
  tempbatches,
  tempconds)
temp_psedobulk_DESeq2[mygenes, ]

##         baseMean log2FoldChange     lfcSE      stat       pvalue       padj
## CCL3L3  66.30936      -2.113009 0.6147236 -3.437332 5.874752e-04 0.09164613
## ICAM1   14.97278      -1.636157 0.3949447 -4.142750 3.431655e-05 0.02462556
## HBB    507.53113      -6.708014 2.8910296 -2.320285 2.032545e-02 0.43828326

## t statistics
delta_poiglm_gw <- do.call(what = cbind,
  args = list(condeff_poiglm_gw$CCL3L3$delta,
    condeff_poiglm_gw$ICAM1$delta,
    condeff_poiglm_gw$HBB$delta))
colnames(delta_poiglm_gw) <- mygenes
t_poiglm_gw <- myt$calt(delta_poiglm_gw)
##    CCL3L3     ICAM1       HBB
## 2.2532804 2.3067713 0.4124514

delta_nbglm_gw <- do.call(what = cbind,
  args = list(condeff_nbglm_gw$CCL3L3$delta,
    condeff_nbglm_gw$ICAM1$delta,
    condeff_nbglm_gw$HBB$delta))
colnames(delta_nbglm_gw) <- mygenes
t_nbglm_gw <- myt$calt(delta_nbglm_gw)

##    CCL3L3     ICAM1       HBB
## 2.2682012 2.2492147 0.3861089
