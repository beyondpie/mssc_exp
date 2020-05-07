library(rstan)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(Seurat)
library(bayesplot)

## * meta settings
bayesplot::color_scheme_set("red")

## * get gene names
load("../from_avi/20200504/deseq.dt.RData")
topgnum <- 100
mygenes <- deseq.dt$gene[1:topgnum]

## * view gene read counts in cells
gse145281 <- readRDS("../from_avi/20200504/seurat.RDS")
mycluster <- 2
VlnPlot(
  object = gse145281, features = "HBA1",
  group.by = "patient", idents = mycluster
)

## * utilies
deg <- function(logc = c(0.0, 1), thres = 0) {
  if (thres < logc[1]) {
    return(-1)
  } else if (thres > logc[2]) {
    return(1)
  } else {
    return(0)
  }
}

stan_analyzer <- function(model, prefix = "./05062300/") {
  model_fnm <- paste0(prefix, "result/", model, ".csv")
  fit <- read_stan_csv(model_fnm)
  ## plot(fit, pars = "mu_g_di")
  mu_ic_plot <- mcmc_areas(fit,
    par = paste0(
      "mu_g_ic",
      paste0(paste0("[", seq(1, 10)), "]")
    ),
    prob = 0.8
  ) +
    theme(axis.text = element_text(size = 12, face = "bold", color = "black")) +
    ggtitle("Individual effect posterior", "with medians and 80% intervals.")

  fitdf <- as.data.frame(fit)

  a <- rstan::extract(fit, pars = "mu_g_di")$mu_g_di
  log_mu_d_change <- as.data.frame(a[, , 1] - a[, , 2])

  colnames(log_mu_d_change) <- mygenes
  quants <- c(0.025, 0.975)
  q <- as.data.frame(apply(log_mu_d_change, 2, quantile, probs = quants))

  degs <- apply(q, 2, deg)
  degnms <- mygenes[degs != 0]

  dhist <- ggplot(gather(log_mu_d_change[degnms]), aes(value)) +
    geom_histogram(bins = 10, color="red", fill="red") +
    facet_wrap(~key, scales = "free_x") +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_text(size = 14, face="bold")
    ) +
    ggtitle("Genes with zeros out of the 95% quantitle intervals.",
            "Hist of the differences of log fold changes.")


  p <- ggarrange(mu_ic_plot, dhist, nrow=1, ncol=2, widths=c(1,2))
  p <- annotate_figure(p,
    top = text_grob(model, color = "blue", face = "bold", size = 16)
    )
  ggsave(paste0(model, ".pdf"), p)
  return(list(
    p = p, pic = mu_ic_plot,
    dhist = dhist, qt = q,
    degnms = degnms
  ))
}

## * analyze the results
r3_0 <- stan_analyzer("v3_0")
r3_1 <- stan_analyzer("v3_1")
r3_2 <- stan_analyzer("v3_2")
r3_3 <- stan_analyzer("v3_3")
r3_4 <- stan_analyzer("v3_4")
r3_5 <- stan_analyzer("v3_5")
