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
output_stan_fnm <- "pbmc_naivecd4t_top1000_stan_v3-1_vi.rds"

## * load cmdstan
library(cmdstanr)
set_cmdstan_path(path = paste(Sys.getenv("HOME"), "softwares",
  paste0("cmdstan-", cmdstan_version), sep = "/"))

library(bayesplot)
library(posterior)

## * functions
get_stan_variable_nms <- function(geneid, numinds, varnm = "MuCond") {
  ## get variable names in stan models
  ## given the gene id and the individual number.
  ## Return: vector of strings.

  invisible(vapply(seq_len(numinds),
    FUN = function(j) {
      paste0(varnm, "[", geneid, ",", j, "]")
    },
    FUN.VALUE = varnm))
}

get_posterior_condiff <- function(stan_res_df, genenms, geneindex,
                                  numconds = 2, varnm = "MuCond") {
  ## get the posterior diff (first - second) between different conditions
  ## Assume 2 in this case
  ## Return list of two elements:
  ## conds: sample_num x 3 colm(Cond1, Cond2, delta)
  ## gds: gene deltas: dataframe sample_num x genenms
  ## posterior samples.

  res <- lapply(geneindex,
    FUN = function(i) {
      indcols <- get_stan_variable_nms(i, numconds, varnm)
      tmp <- stan_res_df[, indcols]
      colnames(tmp) <- paste0(varnm, seq_len(numconds))
      tmp$delta <- (tmp[, 1] - tmp[, 2])[, 1]
      invisible(tmp)
    })
  names(res) <- genenms
  deltas <- vapply(res, FUN = function(x) {
    invisible(x$delta)
  }, FUN.VALUE = res[[1]]$delta)
  ## gd <- do.call(cbind, deltas)
  invisible(list(conds = res, gds = deltas))
}

get_posterior_indeff <- function(stan_res_df, genenms, geneindex,
                                 numinds = 10, varnm = "MuInd") {
  ## get the posterior individual effects.
  ## Return
  ## indeff: List of gene num, each is a df/tibble: sample_num x numinds.

  res <- lapply(geneindex,
    FUN = function(i) {
      indcols <- get_stan_variable_nms(i, numinds, varnm)
      tmp <- stan_res_df[, indcols]
      colnames(tmp) <- paste0(varnm, seq_len(numinds))
      invisible(tmp)
    })
  names(res) <- genenms
  invisible(res)
}


get_posterior_intercept <- function(stan_res_df, genenms,
                                    geneindex, varnm = "MuG") {
  ## get the posterior intercept, i.e. the average
  ## expression level removing
  ## individual effect and conditional effects.
  ## Return
  ## a df: sample_num x gene_num
  ## this function can be used for variance summary
  ## i.e., Kappa2G, Tau2G, Phi2G

  res <- vapply(geneindex,
    FUN = function(i) {
      stan_res_df[, paste0(varnm, "[", i, "]"), drop = T]
    }, FUN.VALUE = stan_res_df[, paste0(varnm,"[1]"), drop = T])
  colnames(res) <- genenms
  invisible(res)
}

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
mssc_stanvi_obj <- mssc_stan_model$variational(
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
mssc_res_df <- posterior::as_draws_df(mssc_stanvi_obj$draws())
mucond_sum <- get_posterior_condiff(stan_res_df = mssc_res_df,
  genenms = rownames(mssc_cnt),
  geneindex = seq_len(nrow(mssc_cnt)),
  numconds = 2,
  varnm = "MuCond")
mucond_t <- myt$calt(mucond_sum$gd)

muind_sum <- get_posterior_indeff(stan_res_df = mssc_res_df,
  genenms = rownames(mssc_cnt),
  geneindex = seq_len(nrow(mssc_cnt)),
  numinds = 10,
  varnm = "MuInd")

mug_sum <- get_posterior_intercept(stan_res_df = mssc_res_df,
  genenms = rownames(mssc_cnt),
  geneindex = seq_len(nrow(mssc_cnt)),
  varnm = "MuG")

kappa2g_sum <- get_posterior_intercept(stan_res_df = mssc_res_df,
  genenms = rownames(mssc_cnt),
  geneindex = seq_len(nrow(mssc_cnt)),
  varnm = "Kappa2G")

tau2g_sum <- get_posterior_intercept(stan_res_df = mssc_res_df,
  genenms = rownames(mssc_cnt),
  geneindex = seq_len(nrow(mssc_cnt)),
  varnm = "Tau2G")

phi2g_sum <- get_posterior_intercept(stan_res_df = mssc_res_df,
  genenms = rownames(mssc_cnt),
  geneindex = seq_len(nrow(mssc_cnt)),
  varnm = "Phi2G")

saveRDS(object = list(
  mucond = mucond_sum,
  mucondt = mucond_t,
  muind = muind_sum,
  mug = mug_sum,
  kappa2g = kappa2g_sum,
  tau2g = tau2g_sum,
  phi2g = phi2g_sum),
file = here("exps", "pbmc", "vi", output_stan_fnm))
