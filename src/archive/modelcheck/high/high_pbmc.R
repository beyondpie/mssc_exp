## fit high model for pbmc dataset

## * load R env
suppressPackageStartupMessages(library(tidyverse))
library(MCMCpack)
library(cmdstanr)
library(grid)
library(gtable)
library(gridExtra)
library(bayesplot)
library(posterior)
library(bbmle)
library(sads)
library(truncnorm)

## warnings/errors traceback settings
options(error = traceback)
options(warn = 1)
options(mc.cores = 3)

high <- modules::import("mssc_hbnb_indeff_across_gene")
options("import.path" = here::here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## * configs
cell_type <- "Naive CD4+ T"

myggtitle <- theme(
  plot.title = element_text(
    size = 15,
    hjust = 0.5, family = "Helvetica", color = "black",
    face = "bold"
  ),
  axis.title = element_text(size = 13, family = "Helvetica", color = "black"),
  axis.text = element_text(size = 13, family = "Helvetica", color = "black"),
  legend.text = element_text(size = 13, family = "Helvetica", color = "black"),
  legend.title = element_text(size = 13, family = "Helvetica", color = "black")
)

select_data_rankby_pseudobulk <- function(num_top_gene = 100,
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
    subscdata$inds, factor(subscdata$resp)
  )
  na_index_from_pseudo <- which(is.na(pseudo_analysis$pvalue) == TRUE)
  pseudo_analysis$pvalue[na_index_from_pseudo] <- 1.0

  ## ranking by pvalue
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
    "high",
    "high_snb_pool_ref_pbmc.rds"
  ))
  return(invisible(r))
}

fit_high <- function(pbmc, j = 2) {
  k <- max(pbmc$ind)
  g <- nrow(pbmc$y2c)
  hi_params <- high$set_hi_params(
                      k = k,
                      j = j,
                      g = g,
                      cnt = pbmc$y2c,
                      s = pbmc$s,
                      cond = pbmc$cond,
                      ind = pbmc$ind,
                      scale = 1.96^2
                    )
  data <- high$to_hbnb_data(pbmc$y2c, ind = pbmc$ind,
                            cond = pbmc$cond, s = pbmc$s,
                            hp = hi_params$hp)
  ip <- hi_params$ip
  vifit <- high$run_hbnb_vi(data = data, ip = ip)
  return(invisible(list(vifit = vifit, data = data,
                        ip = ip, hip = hi_params)))
}

get_pbmc_demo <- function() {
  pbmc <- select_data_rankby_pseudobulk(num_top_gene = 20)
  highfit <- fit_high(pbmc)
  est_params <- highfit$extract_all_params_from_fit(highfit, data)
  saveRDS(object = est_params, file = here::here(
                                              "src", "modelcheck",
                                              "high",
                                              "high_est_params_ref_pbmc.rds"
                                            ))
}

get_high_ranked_scores <- function(rank_stats, from = 200) {
  ranked_index <- order(rank_stats, descreasing = TRUE)
}
