## fit scaled negative binomial for pbcm dataset.
## used for model parameters settings.

## * load R env
source("set_r_lib_env.R")

## * load param fitting functions
source("mssc_param_fitting.R")

## local modules
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## * configs
celltype <- "Naive CD4+ T"
num_top_gene <- 500


## * load stan models
snbm_for_mur <- cmdstanr::cmdstan_model(
  here::here("src", "stan", "scale_nb.stan"),
  compile = T, quiet = F
)

snbm_for_mucond <- cmdstanr::cmdstan_model(
  here::here("src", "stan", "scale_nb_fixed_r.stan"),
  compile = T, quiet = F)

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


## * select top genes based on pseudobulk analysis
## - [NO] use wilcox-test leads the top genes have very few counts.
## - [YES] use deseq2 instead

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

## * estimate mssc parameters
mfit <- fit_multgenes_nb(cnt = cnt, s = sumcnt,
                 cond = resp + 1, ind = inds,
                 snbm = snbm_for_mur,
                 snbm_for_mucond = snbm_for_mucond)
