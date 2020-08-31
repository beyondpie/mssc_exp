## * check if the model reflects the assumptions about the data.
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
# modules::reload(myfit)

## ** PBMC data
## PBMC data are the same as above.
pbmcseurat <- readRDS(paste(datadir, pbmc_IL8_dirnm, "seurat.RDS", sep = "/"))
pbmccnt <- as.matrix(pbmcseurat@assays$RNA@counts)
pbmcinds <- pbmcseurat@meta.data$patient
pbmc_cellanno <- pbmcseurat@meta.data$seurat_clusters
pbmc_cond <- pbmcseurat@meta.data$response

mygenes <- c("HBB", "CCL3L3", "ICAM1")
## cytototic T cells
mycluster <- 2
mulindcells <- pbmc_cellanno == mycluster

mytotcnts <- colSums(pbmccnt[, mulindcells])
mycnts <- pbmccnt[mygenes, mulindcells]
myinds <- pbmcinds[mulindcells]
myconds <- pbmc_cond[mulindcells]
## TODO: remove cells if some genes having too much outliers

d_bagwiff <- myt$to_bagwiff_r(mycnts, myinds, myconds,
  mytotcnts)

myfit_poiglm_mc <- stan(file = "v1-1.stan",
  data = d_bagwiff,
  iter = 1000,
  warmup = 500,
  thin = 1,
  chains = 2,
  refresh = 50,
  seed = 1,
  control = list(adapt_delta = 0.9,
    max_treedepth = 20))

myfit_nbglm_mc <- stan(file = "v3-1.stan",
  data = d_bagwiff,
  iter = 1000,
  warmup = 500,
  thin = 1,
  chains = 2,
  refresh = 50,
  seed = 1,
  control = list(adapt_delta = 0.9,
    max_treedepth = 20)
)

myfit_nblognm_mc <- stan(file = "v4-1.stan",
  data = d_bagwiff,
  iter = 1000,
  warmup = 500,
  thin = 1,
  chains = 2,
  refresh = 50,
  seed = 1,
  control = list(adapt_delta = 0.9,
    max_treedepth = 20))
