## run mssc on pbmc dataset

library(Seurat)
suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)

## local modules
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## warnings/errors traceback settings
options(error = traceback)
options(warn = 0)

## * config
cell_type <- "Naive CD4+ T"
num_top_gene <- 1000L
cmdstan_version <- "2.23.0"
mssc_stan_fnm <- "v3-1.stan"

## * load cmdstan
library(cmdstanr)
set_cmdstan_path(path = paste(Sys.getenv("HOME"), "softwares",
  paste0("cmdstan-", cmdstan_version), sep = "/"))

## * load data
pbmc_seurat <- mypbmc$load_pbmc_seurat() %>%
  mypbmc$extract_from_seurat(pbmc_seurat = .)

## limit to a given cell type
## TODO: check outliers
subscdata <- mypbmc$get_celltype_specific_scdata(pbmc_seurat$cnt,
  pbmc_seurat$resp,
  pbmc_seurat$inds,
  pbmc_seurat$ct,
  "CXCR4+ NK")

## * data preprocessing
## sample cell numbers to keep each ind
## has the same num of cells
sample_cells <- mypbmc$sample_cells_per_ind(subscdata$inds,
  mypbmc$get_minum_of_inds(subscdata$inds))

## note individual order is changed according to sample_cells
cnt <- subscdata$cnt[, sample_cells]
inds <- subscdata$inds[sample_cells]
resp <- subscdata$resp[sample_cells]
colsumcnt <- colSums(cnt)

## use pseudobulk + wilcox-test
pseudobulk <- mypseudo$get_pseudobulk(cnt, inds)
names(resp) <- inds
pseudoconds <- resp[colnames(pseudobulk)]

diff_wilcoxp <- apply(pseudobulk, 1,
  FUN = function(x, group) {
    p <- wilcox.test(x ~ group)$p.value
    if (is.nan(p)) {
      1
    }
    else {
      p
    }
  }, group = pseudoconds)

## select top ranked genes
top_ranked_index <- order(diff_wilcoxp, decreasing = FALSE)[1:num_top_gene]

## final data
mssc_cnt <- cnt[top_ranked_index, ]
mssc_inds <- inds
mssc_resp <- resp
mssc_totcnt <- colsumcnt

## * mssc
## TODO: setting the initial values for the models.

## pack up the data for stan usage
mssc_stan_data <- myt$to_bagwiff_r(mssc_cnt, mssc_inds,
  mssc_resp, mssc_totcnt)
## load mssc stan model
mssc_stan_model <- cmdstanr::cmdstan_model(
  stan_file = here("src", "stan", mssc_stan_fnm),
  compile = T,
  stanc_options = list(include_paths = here("src", "stan")),
  dir = here("exps", "stanbin")
  )
## using stan vi method
mssc_stan_model$variational(
                  data = mssc_stan_data,
                  seed = 355133L,
                  init = NULL,
                  refresh = 200,
                  algorithm = "meanfield",
                  output_samples = 1000,
                  save_latent_dynamics = F,
                  output_dir = here("exps", "pbmc", "vi")
                  )
## result analysis

