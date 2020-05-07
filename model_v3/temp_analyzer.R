library(rstan)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(Seurat)

prefix <- "./05062300/"
resfnm <- "v3_1.csv"

params <- c("mu_g_ic", "mu", "mu_g_di")

model_fnm <- paste0(prefix, "result/", resfnm)
fit <- read_stan_csv(model_fnm)
plot(fit, pars="mu_g_di")
fitdf <- as.data.frame(fit)

#allparams <- fit@par_dims
#mu_ic <- allparams$

a <- extract(fit, pars="mu_g_di")$mu_g_di
log_mu_d_change <- as.data.frame(a[, , 1] - a[,,2])

load("../from_avi/20200504/deseq.dt.RData")
topgnum <- 100
mygenes <- deseq.dt$gene[1:topgnum]

colnames(log_mu_d_change) <- mygenes
# log_mu_d_change %>% gather() %>% dim()
ggplot(gather(log_mu_d_change), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

quants <- c(0.025, 0.975)
q <- as.data.frame(apply(log_mu_d_change, 2, quantile, probs=quants))

deg <- function(logc=c(0.0, 1), thres= 0) {
  if (thres < logc[1]) {
    return(-1)
  } else if (thres > logc[2]) {
    return(1)
  } else {
    return(0)
  }
}

degs <- apply(q, 2, deg)

degnms <- colnames(q)[degs != 0]

ggplot(gather(log_mu_d_change[degnms]), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

gse145281 <- readRDS("../from_avi/20200504/seurat.RDS")
mycluster <- 2
VlnPlot(
  object = gse145281, features = "HBA1",
  group.by = "patient", idents = mycluster
)
