## fit scaled negative binomial for pbcm dataset.
## used for model parameters settings.

## * load R env
source("hbnb_set_r_lib_env_01.R")

## * load param fitting functions
source("hbnb_param_fitting_02.R")

## local modules
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## * configs
celltype <- "Naive CD4+ T"
num_top_gene <- 10


## * load stan models
snbm_for_mur <- cmdstanr::cmdstan_model(
  here::here("stanutils", "scale_nb.stan"),
  compile = T, quiet = F
)

snbm_for_mucond <- cmdstanr::cmdstan_model(
  here::here("stanutils", "scale_nb_fixed_r.stan"),
  compile = T, quiet = F
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

main <- function() {
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
    decreasing = FALSE
  )[1:num_top_gene]
  pvalue_pseudo_deseq2 <- pseudo_analysis$pvalue[top_ranked_index]

  ## finally used pbmc data
  ## use all the counts to get sumcnt
  sumcnt <- colSums(subscdata$cnt)
  inds <- subscdata$inds
  resp <- subscdata$resp
  cnt <- subscdata$cnt[top_ranked_index, ]

  ## consistent with stan
  resp <- resp + 1
  ind <- vapply(inds, function(nm) {
    for (i in 1:5) {
      if (nm == str_glue("NR{i}")) {
        return(invisible(as.integer(i)))
      }
      if (nm == str_glue("R{i}")) {
        return(invisible(as.integer(i) + 5L))
      }
    }
  }, 0L)

  ## scale the total number of count.
  s <- sumcnt / median(sumcnt)

  ## * estimate mssc parameters
  mg_snb_mat <- fit_mg_snb(
    cnt = cnt, s = s,
    cond = resp, ind = ind,
    snbm = snbm_for_mur,
    snbm_for_mucond = snbm_for_mucond
  )

  ## * save the data as a pool
  r <- list(
    sumcnt = sumcnt,
    s = s,
    ind = ind,
    cond = resp,
    y2c = cnt,
    est_snb_mat = mg_snb_mat
  )
  saveRDS(object = r, file = here::here(
    "src",
    "modelcheck",
    "snb_pool_ref_pbmc.rds"
  ))
}

## run main
main()
