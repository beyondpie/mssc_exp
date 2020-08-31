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

mygenes <- c("HBB", "CCL3L3", "ICAM1")
## cytototic T cells
mycluster <- 2
mulindcells <- pbmc_cellanno == mycluster

mytotcnts <- colSums(pbmccnt[, mulindcells])
mycnts <- pbmccnt[mygenes, mulindcells]
myinds <- pbmcinds[mulindcells]
myconds <- pbmc_cond[mulindcells]

## TODO: reove cells if some genes have outliers
standata_bagwiff <- myt$to_bagwiff_r(
  mycnts, myinds,
  myconds,
  mytotcnts
)

## * setup model
get_mystanmodel <- function(fnm) {
  cmdstanr::cmdstan_model(
    stan_file = here("src", "stan", fnm),
    compile = T,
    stanc_options = list(include_paths = here("src", "stan")),
    dir = here("src", "modelcheck", "stan_binary")
  )
}

sampling_mystanmodel <- function(mystanmodel, num_chains = 4) {
  stanmodel_poiglm_gw_mc$sample(data = standata_bagwiff,
    seed = 355113L,
    refresh = 200,
    iter_warmup = 500,
    iter_sampling = 500,
    adapt_delta = 0.85,
    max_treedepth = 12,
    validate_csv = T,
    show_messages = T,
    chains = num_chains,
    parallel_chains = getOption("mc.cores", num_chains),
    output_dir = here("src", "modelcheck",
      "stan_sampling")
  )
}

stanmodel_poiglm_gw_mc <- get_mystanmodel("v1-1.stan")

stanmodel_nbglm_gw_mc <- get_mystanmodel("v3-1.stan")

stanmodel_nblognm_gw_mc <- get_mystanmodel("v4-1.stan")

## * sampling
num_chains <- 4
stanfit_poiglm_gw_mc <- sampling_mystanmodel(stanmodel_poiglm_gw_mc, num_chains)
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
