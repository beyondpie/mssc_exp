library(rstan)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(Seurat)
library(bayesplot)
library(data.table)
import::from(mltools, one_hot)
## * meta settings
bayesplot::color_scheme_set("blue")

## * get gene names
load("../from_avi/20200504/deseq.dt.RData")
## topgnum <- 100
## mygenes <- deseq.dt$gene[1:topgnum]

## * get experimental genes
utr <- 8
## use low rank
ulr <- 6
goldgenes <- data.frame(genes = c(
  "SNHG16",
  "OASL",
  "NAMPT",
  "NFKB1",
  "BCL2L11",
  "IRF8",
  "TPM4",
  "TRAF4",
  "ICAM1", "XCL2", "XCL1",
  "RPS26P11", "LOC101929876", "LOC100996747",
  "HBA1", "HBA2", "HBB", "HBD",
  "CCL3L3", "CCL3L1", "CCL3",
  "KDM6A",
  "ZNF721",
  "HDDC2",
  "YIPF5",
  "MAK16",
  "TOX"
), module = c(
  seq(1, utr), rep(utr + 1, 3),
  rep(utr + 2, 3), rep(utr + 3, 4),
  rep(utr + 4, 3), seq(utr + 5, utr + 5 + ulr - 1)
))
module <- factor(goldgenes$module)
GoldB <- as.matrix(one_hot(data.table(module)))
rownames(GoldB) <- goldgenes$genes

## ** get counts matrix and design matrix
mygenes <- goldgenes$genes
B <- GoldB


## * view gene read counts in cells
gse145281 <- readRDS("../from_avi/20200504/seurat.RDS")
mycluster <- 2
VlnPlot(
  object = gse145281, features = "HBA1",
  group.by = "patient", idents = mycluster
)

## * utilies
draw_invgamma <- function(alpha, rate = 1 / alpha, from = 0.0001, to = 20) {
  library(invgamma)
  library(ggplot2)
  theme_set(theme_bw())
  x <- seq(from, to, 0.01)
  qplot(x, dinvgamma(x, alpha, rate), geom = "line")
}

deg <- function(logc = c(0.0, 1), thres = 0) {
  if (thres < logc[1]) {
    return(-1)
  } else if (thres > logc[2]) {
    return(1)
  } else {
    return(0)
  }
}

stan_analyzer <- function(model, prefix = "./result/", 
                          desc="", color="red") {
  model_fnm <- paste0(prefix, model, ".csv")
  csvfiles <- dir(path="./result/", pattern=paste0(model,'[0-9].csv'),full.names=TRUE)
  fit <- read_stan_csv(csvfiles)
  ## plot(fit, pars = "mu_g_di")
  ## mu_ic_plot <- mcmc_areas(fit,
  ##   par = c("MuCond[15,1]", "MuCond[15,2]"),
  ##   prob = 0.8
  ## ) +
  ##   theme(axis.text = element_text(size = 12, face = "bold", color = "black")) +
  ##   ggtitle("Individual effect posterior", "with medians and 80% intervals.")

  fitdf <- as.data.frame(fit)

  a <- rstan::extract(fit, pars = "MuCond")$MuCond
  log_mu_d_change <- as.data.frame(a[, , 1] - a[, , 2])

  colnames(log_mu_d_change) <- mygenes
  quants <- c(0.025, 0.975)
  q <- as.data.frame(apply(log_mu_d_change, 2, quantile, probs = quants))

  degs <- apply(q, 2, deg)
  degnms <- mygenes[degs != 0]

  dhist <- ggplot(gather(log_mu_d_change[degnms], factor_key = TRUE), aes(value)) +
    geom_histogram(bins = 20, color=color, fill=color) +
    facet_wrap(~key, scales = "free_x", nrow =2) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10),
      strip.text.x = element_text(size = 12, face="bold")
    ) +
    ggtitle( paste0("Differential genes detected by ", desc),
            "Histogram of the differences of log fold changes between control and case.") +
    theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size=12))


  ## p <- ggarrange(mu_ic_plot, dhist, nrow=1, ncol=2, widths=c(1,2))
  ## p <- annotate_figure(p,
  ##   top = text_grob(model, color = "blue", face = "bold", size = 16)
  ##   )
  ## ggsave(paste0(model, ".pdf"), p)
  return(list(
    ## p = p, pic = mu_ic_plot,
    fit = fit,
    dhist = dhist, qt = q,
    degnms = degnms
  ))
}

## * analyze the results
v11_stan <- stan_analyzer("v1-1",
                          desc="modeling gene-wise invididual effect indepdently",
                          color="red")
v12_stan <- stan_analyzer("v1-2",
                          desc="modeling gene-wise invididual effect, and sharing hyper prior for variances.",
                          color="blue")

p1 <- ggarrange(v11_stan$dhist, v12_stan$dhist, nrow = 2, ncol=1)

v21_stan <- stan_analyzer("v2-1",
                          desc="modeling gene-module invididual effect indepdently",
                          color="red")
v22_stan <- stan_analyzer("v2-2",
                          desc="modeling gene-module invididual effect, and sharing hyper prior for variances.",
                          color="blue")

p2 <- ggarrange(v21_stan$dhist, v22_stan$dhist, nrow = 2, ncol=1)

v22_stan_rep <- stan_analyzer("v2-2",
                              prefix = "./rep_result/", 
                          desc="modeling gene-module invididual effect, and sharing hyper prior for variances.",
                          color="blue")
## * vln plot
pos_sgenes <- c("SNHG16",
  "OASL",
  "NAMPT",
  "NFKB1",
  "BCL2L11",
  "IRF8",
  "TPM4",
  "TRAF4")

IXXmodule <- c("ICAM1","XCL2","XCL1")
RLLmodule <- c("RPS26P11","LOC101929876","LOC100996747")

HBmodule <- c("HBA1", "HBA2", "HBB","HBD")

CCCmodule <- c("CCL3L3","CCL3L1","CCL3")

neg_sgenes <- c( "KDM6A",
                           "ZNF721",
                           "HDDC2",
                          "YIPF5",
                           "MAK16",
                           "TOX")

gse145281 <- readRDS("../from_avi/20200504/seurat.RDS")
scmeta <- gse145281@meta.data
sccluster <- scmeta$seurat_clusters
mycluster <- 2
mycells <- which(sccluster == mycluster)

pos_svln <-VlnPlot(object = gse145281, features = pos_sgenes,
                         group.by = "patient", idents = mycluster)
neg_svln <- VlnPlot(object = gse145281, features = neg_sgenes,
                          group.by = "patient", idents = mycluster)

ixx_vln <-VlnPlot(object = gse145281, features = IXXmodule,
                         group.by = "patient", idents = mycluster)
rll_vln <-VlnPlot(object = gse145281, features = RLLmodule,
                         group.by = "patient", idents = mycluster)
hb_vln <-VlnPlot(object = gse145281, features = HBmodule,
                  group.by = "patient", idents = mycluster)
ccc_vln <- VlnPlot(object = gse145281, features = CCCmodule,
                   group.by = "patient", idents = mycluster)

## * test vi

v22vi <- stan_analyzer("v2-2",
                              prefix = "./", 
                              desc="modeling gene-module invididual effect, and sharing hyper prior for variances.",
                              color="blue")

v21vi <- stan_analyzer("v2-1",
                              prefix = "./result/vi/", 
                              desc="modeling gene-module invididual effect, and sharing hyper prior for variances.",
                              color="blue")
