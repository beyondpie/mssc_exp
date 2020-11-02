## Check MSSC Hierarhical Bayesian Model

## Instead of using SymSim to simulate the dataset
## Here we firstly use the data simulated from our model.
## Hyper parameters are learned from the real dataset: PBMC.

## * set R environment
library(optparse)
import::from(here, here)
suppressPackageStartupMessages(library(tidyverse))
library(MCMCpack)

## develop version of cmdstanr
## devtools::install_github("stan-dev/cmdstanr")
library(cmdstanr)
## use the code below to reload the library
## detach("package:cmdstanr", unload = TRUE)
## library(cmdstanr)

library(grid)
library(gtable)
library(gridExtra)
library(bayesplot)
## for bayesplot plotting
## color_scheme_set("brewer-Spectral")

library(posterior)
library(bbmle)
library(sads)
## install_github("olafmersmann/truncnorm")
library(truncnorm)

## local modules
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
## use the code below to reload the modules
## modules::reload(myfit)
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## warnings/errors traceback settings

## options(error = traceback)
options(warn = 0)
options(mc.cores = 3)

## * configs
## for simulation
opts <- list(
  make_option(c("--ncell", action = "store", type = "integer", default = 200)),
  make_option(c("--nind", action = "store", type = "integer", default = 10)),
  make_option(c("--celltype",
    action = "store", type = "character",
    default = "Naive CD4+ T"
  ))
)

args <- parse_args(OptionParser(option_list = opts))
num_of_cell_per_ind <- args$ncell
num_of_ind <- args$nind
num_of_cond <- 2
num_of_ind_per_cond <- floor(num_of_ind / num_of_cond)
celltype <- args$celltype
num_top_gene <- 200

## for debug
num_of_cell_per_ind <- 200
num_of_ind <- 10
num_of_cond <- 2
num_of_ind_per_cond <- floor(num_of_ind / num_of_cond)
celltype <- "Naive CD4+ T"
num_top_gene <- 200


## * load stan models
## if compiling not work, try rebuild_cmdstan()

scale_nb_model <- cmdstan_model(
  here::here("src", "stan", "scale_nb.stan")
)
mssc_hbnb_rndeff_model <- cmdstan_model(
  here::here("src", "dirty_stan", "hbnb_rndeff.stan"),
  compile = T, quiet = FALSE
)

## * load pbmc dataset
pbmc_seurat <- mypbmc$load_pbmc_seurat() %>%
  mypbmc$extract_from_seurat(pbmc_seurat = .)
## limit to a cell type
subscdata <- mypbmc$get_celltype_specific_scdata(
  pbmc_seurat$cnt,
  pbmc_seurat$resp,
  pbmc_seurat$inds,
  pbmc_seurat$ct,
  ## limit to the cell type
  celltype
)

## * estimate hyper parameters
## select top 1000 genes based on pseudobulk analysis

## use pseudobulk + wilcox-test
pseudobulk <- mypseudo$get_pseudobulk(subscdata$cnt, subscdata$inds)
names(subscdata$resp) <- subscdata$inds
pseudoconds <- subscdata$resp[colnames(pseudobulk)]

pseudo_analysis <- mypseudo$pseudobulk_deseq2(
  subscdata$cnt,
  subscdata$inds, subscdata$resp
)
na_index_from_pseudo <- which(is.na(pseudo_analysis$pvalue) == TRUE)
pseudo_analysis$pvalue[na_index_from_pseudo] <- 1.0

top_ranked_index <- order(pseudo_analysis$pvalue,
                          decreasing = FALSE)[1:num_top_gene]
pvalue_pseudo_deseq2 <- pseudo_analysis$pvalue[top_ranked_index]

## finally used pbmc data
## use all the counts to get sumcnt
sumcnt <- colSums(subscdata$cnt)
inds <- subscdata$inds
resp <- subscdata$resp
cnt <- subscdata$cnt[top_ranked_index, ]

## * functions
fit_scalenb_stan <- function(s, y, scale_nb_model,
                             seed = 355113, numiter = 5000,
                             refresh = 500, r_default = 10) {
  ## mu in scaled log level, and minus log(s)
  ## use stan to fit

  result <- list(mu = NaN, r = NaN)
  # fit scaled negative binomial using stan
  n <- length(s)
  ## ref: MASS::fitdistr for nb
  m <- mean(y)
  v <- var(y)
  r <- if (v > m) {
    m^2 / (v - m)
  } else {
    r_default
  }
  opt <- scale_nb_model$optimize(
    data = list(n = n, s = s, y = y),
    seed = seed,
    refresh = refresh,
    iter = numiter,
    init = list(list(
      mu = log(m) - log(median(s)),
      r = r
    )),
    algorithm = "lbfgs"
    )
  if (!opt$runset$procs$is_finished()) {
    message(str_glue("Gene [gn]: optimization is not finished."))
  } else if (opt$runset$procs$is_failed()) {
    message(str_glue("Gene [gn]: optimizaiton failed."))
  } else {
    t <- opt$mle()
    result$mu <- t["mu"]
    result$r <- t["r"]
  }

  invisible(result)
}

fit_singlegene_nb <- function(gn, cnt, resp, sumcnt,
                              scale_nb_model,
                              seed = 355113,
                              id_control = 0) {
  ## fit nb distribute to get scaled mu0 and nb dispersion called r.
  ## may fail to fit due to data specificity.
  ## gn can gene name or gene id.

  result <- list(mu0 = NaN, r0 = NaN,
    mu_cond = c(NaN, NaN),
    r_cond = c(NaN, NaN))


  y <- cnt[gn, ]
  outliers <- myfit$is_outlier(y)
  y <- y[!outliers]
  conds <- resp[!outliers]
  s <- sumcnt[!outliers]

  ## fit nb distr for mu0
  mur0 <- fit_scalenb_stan(s, y, scale_nb_model, seed)
  result$mu0  <- mur0$mu
  result$r0 <- mur0$r

  ## fit nb distr for mu_cond
  y_control <- y[conds == id_control]
  s_control <- s[conds == id_control]

  y_case <- y[conds != id_control]
  s_case <- s[conds != id_control]

  mur_control <- fit_scalenb_stan(s_control, y_control, scale_nb_model, seed)
  mur_case <- fit_scalenb_stan(s_case, y_case, scale_nb_model, seed)

  result$mu_cond <- c(mur_control$mu, mur_case$mu)
  result$r_cond <- c(mur_control$r, mur_case$r)
  invisible(result)
}

fit_multgenes_nb <- function(cnt, resp, sumcnt,
                             scale_nb_model,
                             seed = 355113,
                             id_control = 0) {
  result <- lapply(seq_len(nrow(cnt)),
                   FUN = function(i) {
                     unlist(fit_singlegene_nb(i, cnt, resp, sumcnt,
                                       scale_nb_model, seed, id_control))
                   })
  matres <- do.call(rbind, result)
  rownames(matres) <- rownames(cnt)
  invisible(matres)
}

