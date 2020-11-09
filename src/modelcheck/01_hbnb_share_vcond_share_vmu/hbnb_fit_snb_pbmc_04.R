## fit scaled negative binomial for pbcm dataset.
## used for model parameters settings.

## * load R env
source("hbnb_set_r_lib_env_01.R")
source("hbnb_param_fitting_02.R")

## local modules
hbnb_mssc <- modules::import("hbnb_mssc_03")
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

select_data_rankby_pseudobulk <- function(num_top_gene = 10,
                                          celltype = "Naive CD4+ T") {
  ## load pbmc dataset
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
  ## select top genes based on pseudobulk analysis
  ## - [NO] use wilcox-test leads the top genes have very few counts.
  ## - [YES] use deseq2 instead

  pseudo_analysis <- mypseudo$pseudobulk_deseq2(
    subscdata$cnt,
    subscdata$inds, subscdata$resp
  )
  na_index_from_pseudo <- which(is.na(pseudo_analysis$pvalue) == TRUE)
  pseudo_analysis$pvalue[na_index_from_pseudo] <- 1.0

  top_ranked_index <- order(pseudo_analysis$pvalue,
    decreasing = FALSE
  )[1:num_top_gene]
  ## pvalue_pseudo_deseq2 <- pseudo_analysis$pvalue[top_ranked_index]

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

  ## estimate mssc parameters
  mg_snb_mat <- fit_mg_snb(
    cnt = cnt, s = s,
    cond = resp, ind = ind,
    snbm = hbnb_mssc$snbm_for_mur,
    snbm_for_mucond = hbnb_mssc$snbm_for_mucond
  )

  ## save the data as a pool
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
  return(invisible(r))
}

get_vifit <- function(pbmc, j = 2) {
  k <- max(pbmc$ind)
  g <- nrow(pbmc$y2c)
  hi_params <- hbnb_mssc$set_hi_params(
    k = k, j = j, g = g,
    cnt = pbmc$y2c, s = pbmc$s,
    cond = pbmc$cond,
    ind = pbmc$ind,
    scale = 1.96^2
  )
  data <- hbnb_mssc$to_hbnb_data(pbmc$y2c,
    ind = pbmc$ind,
    cond = pbmc$cond, s = pbmc$s,
    hp = hi_params$hp
  )
  ip <- hbnb_mssc$calibrate_init_params(hi_params$ip, data = data)
  vifit <- hbnb_mssc$run_hbnb_vi(data = data, ip = ip)
  return(invisible(list(
    vifit = vifit, data = data,
    ip = ip, hip = hi_params
  )))
}

est_hbnb_params <- function(vifit, data) {
  est_params <- lapply(hbnb_mssc$nm_params, function(nm) {
    hbnb_mssc$extract_vifit(vifit, data, nm)
  })
  names(est_params) <- hbnb_mssc$nm_params
  saveRDS(object = est_params, file = here::here(
    "src", "modelcheck", "est_hbnb_params_ref_pbmc.rds"
  ))
  return(invisible(est_params))
}

## run main
pbmc <- select_data_rankby_pseudobulk(num_top_gene = 20)
vifit <- get_vifit(pbmc, 2)
est_params <- est_hbnb_params(vifit$vifit, vifit$data)
