## fit scaled negative binomial for pbcm dataset.
## used for model parameters settings.

## * load R env
source("hbnb_set_r_lib_env_01.R")
source("hbnb_param_fitting_02.R")
## * other packages
library(GGally)
library(hrbrthemes)
library(viridis)

## local modules
hbnb_mssc <- modules::import("hbnb_mssc_03")
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")
## * config
cell_type <- "Naive CD4+ T"
num_top_gene <- 200
output_stan_fnm <- "pbmc_naivecd4t_top1000_hbnb_vi.rds"
hbnb_p_epsilon <- 0.1

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
    subscdata$inds, factor(subscdata$resp)
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

get_vifit <- function(pbmc, j = 2, calibrate_with_opt = FALSE) {
  k <- max(pbmc$ind)
  g <- nrow(pbmc$y2c)

  hi_params <- hbnb_mssc$set_hi_params(
    k = k, j = j, g = g,
    cnt = pbmc$y2c, s = pbmc$s,
    cond = pbmc$cond,
    ind = pbmc$ind,
    scale = 1.96^2,
    mg_snb_mat = pbmc$est_snb_mat
    )

  data <- hbnb_mssc$to_hbnb_data(pbmc$y2c,
    ind = pbmc$ind,
    cond = pbmc$cond, s = pbmc$s,
    hp = hi_params$hp
  )
  if (calibrate_with_opt) {
    ip <- hbnb_mssc$calibrate_init_params(hi_params$ip, data = data)
  } else {
    ip <- hi_params$ip
  }
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

## * run main: use for 05 script
main <- function() {
  pbmc <- select_data_rankby_pseudobulk(num_top_gene = 20)
  vifit <- get_vifit(pbmc, 2)
  est_params <- est_hbnb_params(vifit$vifit, vifit$data)
}


## * hbnb pbmc rank
get_hbnb_ranked_scores <- function(rank_stat, from = num_top_gene) {
  ranked_index <- order(rank_stat, decreasing = TRUE)
  ranked_scores <- seq(from = from, to = 1, by = -1)
  names(ranked_scores) <- names(rank_stat[ranked_index])
  return(invisible(ranked_scores))
}

get_comp_figure <- function(pseudo_scores,
                            hbnbz_scores,
                            hbnbp_scores,
                            ingroupnum = 50, genenum = 150,
                            printnm = T, hjust = 1, vjust = 1,
                            text_angle = 45, text_size = 4,
                            alphaLines = 0.3) {
  pseu_genes <- names(pseudo_scores)[seq_len(genenum)]
  msscp_genes <- names(hbnbp_scores)[seq_len(genenum)]

  rank_comp <- data.frame(
    pseudo = pseudo_scores[pseu_genes],
    mssc_z = hbnbz_scores[pseu_genes],
    mssc_p = hbnbp_scores[pseu_genes],
    group = c(rep("Top_pseudo", ingroupnum),
      rep("Middle_pseudo", genenum - 2 * ingroupnum),
      rep("Last_pseudo", ingroupnum)),
    pseu_genenm = pseu_genes,
    msscp_genenm = msscp_genes
  )
  ## rownames(rank_comp) <- genes

  p <- ggparcoord(rank_comp,
    columns = 1:3,
    showPoints = TRUE,
    groupColumn = 4,
    alphaLines = alphaLines,
    title = "Naive CD4+ T cell on PBMC",
    )

  if (printnm) {
    p <- p +
      geom_text(data = rank_comp %>%
        dplyr::select(pseu_genenm, pseudo) %>%
        dplyr::mutate(x = 0.9, y=pseudo),
      aes(x = x, y = scale(pseudo), label = pseu_genenm),
      hjust = hjust,
      vjust = vjust,
      inherit.aes = FALSE,
      angle = text_angle,
      size = text_size) + 
    geom_text(data = rank_comp %>%
                dplyr::select(msscp_genenm, pseudo) %>%
                dplyr::mutate(x = 3.2, y = pseudo),
              aes(x = x, y = scale(pseudo), label = msscp_genenm),
              hjust = hjust,
              vjust = vjust,
              inherit.aes = FALSE,
              angle = text_angle,
              size = text_size)

  }
  p <- p +
    scale_color_manual(values = c("Top_pseudo" = "red",
      "Last_pseudo" = "blue",
      "Middle_pseudo" = "black")) +
    xlab("Method") +
    ylab("Inverse Rank Index") +
    myggtitle
  invisible(p)
}

## get pbmc data and pseudo bulk analysis
pbmc <- select_data_rankby_pseudobulk(num_top_gene = num_top_gene)
pseudo_ranked_scores <- seq(from = num_top_gene, to = 1, by = -1)
names(pseudo_ranked_scores) <- rownames(pbmc$y2c)

## get hbnb vi fit
hbnb_vifit <- get_vifit(pbmc, 2, calibrate_with_opt = FALSE)
est_params <- lapply(hbnb_mssc$nm_params, function(nm) {
  hbnb_mssc$extract_vifit(hbnb_vifit$vifit, hbnb_vifit$data, nm)
})
names(est_params) <- hbnb_mssc$nm_params

rank_stats <- hbnb_mssc$get_rank_statistics(est_params$mu_cond, epsilon = 0.1)

## rank the genes by hbnb statistics
hbnbz_ranked_scores <- get_hbnb_ranked_scores(rank_stats$z, from = num_top_gene)

hbnbp_ranked_scores <- get_hbnb_ranked_scores(rank_stats$p)

## get the plots
p <- get_comp_figure(pseudo_scores = pseudo_ranked_scores,
                         hbnbz_scores = hbnbz_ranked_scores,
                         hbnbp_scores = hbnbp_ranked_scores,
                         genenum = 200, vjust = 1, text_angle = 0,
                         text_size= 1.6, hjust = 1)
pdf(str_glue("pbmc_naivecd4t_top_{num_top_gene}.pdf"),
    width = 10,
    height = 15)
p
dev.off()

## ggsave(path = here("exps", "pbmc", "vi"),
##   filename = str_glue("pbmc_naivecd4t_top_{num_top_gene}.pdf"),
##   device = "pdf",
##   plot = p,
##   width = 7,
##   height = 30)
