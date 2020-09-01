library(cmdstanr)
library(bayesplot)
library(posterior)
library(Seurat)
suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")

options(error = traceback)
options(warn = 0)

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


stanfit_poiglm_gw_mc$draws()

## csv_contents <- read_cmdstan_csv(fit$output_files())
## str(csv_contents)
## fit <- mod$sample(data = data_list, save_latent_dynamics = TRUE)
## fit$latent_dynamics_files()
## x <- utils::read.csv(fit$latent_dynamics_files()[1], comment.char = "#")
## head(x)
## temp_rds_file <- tempfile(fileext = ".RDS") # temporary file just for demonstration
## fit$save_object(file = temp_rds_file)
## fit <- readRDS(temp_rds_file)
## fit$summary()
