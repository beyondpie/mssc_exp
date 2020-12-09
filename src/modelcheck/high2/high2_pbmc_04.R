## fit scaled negative binomial for pbmc dataset.
## used for model parameters settings.

## * load R env
library(optparse)
import::from(here, here)
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
options(warn = 0)
options(mc.cores = 3)

## * other packages
library(GGally)
library(hrbrthemes)
library(viridis)

## local modules
high <- modules::import("high2")
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")
## * config
## cell_type <- "Naive CD4+ T"
cell_type <- "CD8A T"
num_top_gene <- 100

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

  ## save the data as a pool
  r <- list(
    sumcnt = sumcnt,
    s = s,
    ind = ind,
    cond = resp,
    y2c = cnt
  )
  return(invisible(r))
}

## * run main: use for 05 script


## * hbnb pbmc rank
get_hbnb_ranked_scores <- function(rank_stat, from = num_top_gene) {
  ranked_index <- order(rank_stat, decreasing = TRUE)
  ranked_scores <- seq(from = from, to = 1, by = -1)
  names(ranked_scores) <- names(rank_stat[ranked_index])
  return(invisible(ranked_scores))
}

get_comp_figure <- function(pseudo_scores,
                            mssc,
                            mssc_rsis,
                            ingroupnum = 50, genenum = 150,
                            printnm = T, hjust = 1, vjust = 1,
                            text_angle = 45, text_size = 4,
                            alphaLines = 0.3) {
  pseu_genes <- names(pseudo_scores)[seq_len(genenum)]
  msscp_genes <- names(mssc)[seq_len(genenum)]

  rank_comp <- data.frame(
    pseudo = pseudo_scores[pseu_genes],
    mssc_orig = mssc[pseu_genes],
    mssc_with_rsis = mssc_rsis[pseu_genes],
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

## high2
## nind <- 10
nind <- max(pbmc$ind)
ncond <- 2
ngene <- nrow(pbmc$y2c)
cnt <- pbmc$y2c
s <- pbmc$s
cond <- pbmc$cond
ind <- pbmc$ind

mssc <- high$High2$new(
  stan_snb_path = here::here(
    "src", "modelcheck", "high2",
    "stan", "snb.stan"
  ),
  stan_snb_cond_path = here::here(
    "src", "modelcheck", "high2",
    "stan", "snb_cond.stan"
  ),
  stan_high2_path = here::here(
    "src", "modelcheck", "high2",
    "stan", "high2.stan"
  ),
  nind = nind,
  tol_rel_obj = 1e-06,
  adapt_iter = 400,
  adapt_engaged = 0,
  algorithm = "meanfield",
  eta = 0.1
)
init_all_params <- mssc$init_params(
  cnt = cnt,
  s = s,
  cond = cond,
  ind = ind)
data <- mssc$to_model_data(cnt = cnt,
                           s = s,
                           cond = cond,
                           ind = ind,
                           hp = init_all_params$hp)
mssc$run(data = data, list_wrap_ip = list(init_all_params$ip))

mucond <- mssc$extract_draws(
  param = "mucond", ngene = nrow(cnt),
  genenms = rownames(cnt))

rankings <- mssc$get_ranking_statistics(
  mucond = mucond,
  two_hot_vec = c(1, -1)
)
## rank the genes by hbnb statistics
### ranked by abs(mean(diff))
hbnbz_ranked_scores <- get_hbnb_ranked_scores(rankings[, 3], from = num_top_gene)

### ranked by rsis
psis <- mssc$psis()
print(psis$psis$diagnostics$pareto_k)
rsis_rankings <- mssc$get_rsis_ranking_statistics(
  ngene = nrow(cnt),
  two_hot_vec = c(1,-1),
  genenms = rownames(cnt),
  normweights = psis$normweights
)
hbnbp_ranked_scores <- get_hbnb_ranked_scores(rsis_rankings[,3], from = num_top_gene)

## get the plots
p <- get_comp_figure(pseudo_scores = pseudo_ranked_scores,
                         mssc = hbnbz_ranked_scores,
                         mssc_rsis = hbnbp_ranked_scores,
                         genenum = num_top_gene, vjust = 1, text_angle = 0,
                         text_size= 1.6, hjust = 1)
pdf(str_glue("pbmc_naivecd4t_top_{num_top_gene}.pdf"),
    width = 10,
    height = 15)
p
dev.off()
