## Eval hbnb in SymSim simulation.

## * set R environment
source("hbnb_set_r_lib_env_01.R")

hbnbm <- modules::import("hbnb_mssc_03")

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")
mysymsim <- modules::import("mysymsim")


## * symsim data
symsim_data_path <- here::here("data", "symsim", "twostage_be_symsim", "data")
symsim_exp_pah <- here::here("exps", "symsim", "stan")

## * configs
tag <- 3

## * load and calculate ground truth
symsimtrue <- readRDS(file = file.path(
  symsim_data_path,
  str_glue("symsim_true_{tag}.rds")
))

symsim_dea <- mysymsim$symsim_de_analysis(symsimtrue,
  popA_idx = which(symsimtrue$cell_meta$pop == 1),
  popB_idx = which(symsimtrue$cell_meta$pop == 2)
)

ngene <- nrow(symsimtrue$counts)

symsim_degenes <- mysymsim$get_symsim_degenes(symsim_dea,
  nDiffEVF = 1,
  logFC = 0.6
) %>%
  which(. == T)

symsim_ndegs <- setdiff(seq_len(ngene), symsim_degenes)

symsim_strict_ndegenes <- mysymsim$get_symsim_strict_ndegenes(symsim_dea,
  nDiffEVF = 0,
  logFC = 0.5
) %>%
  which(. == T)

symsim_zerodiffevf_genes <- mysymsim$get_symsim_ndiffevf_genes(symsim_dea) %>%
  which(. == T)

## * load simulated read counts
symsimumi <- readRDS(file = file.path(
  symsim_data_path,
  str_glue("symsim_umi_{tag}.rds")
))
symsimbe <- readRDS(file = file.path(
  symsim_data_path,
  str_glue("symsim_be_{tag}.rds")
))
symsim2be <- readRDS(file = file.path(
  symsim_data_path,
  str_glue("symsim_2be_{tag}.rds")
))

## * hbnb analysis
init_params_and_data <- function(symsim2be) {
  y2c <- symsim2be$counts
  ind <- symsim2be$batch_meta$batch
  cond <- symsim2be$cell_meta$pop
  sumcnt <- colSums(y2c)
  s <- sumcnt / median(sumcnt)

  hi_params <- hbnbm$set_hi_params(
    k = max(ind),
    j = max(cond),
    g = nrow(y2c),
    cnt = y2c,
    s = s,
    cond = cond, ind = ind,
    scale = 1.96^2
  )

  data <- hbnbm$to_hbnb_data(y2c, ind, cond, s, hi_params$hp)
  return(invisible(list(hip = hi_params, data = data)))
}

get_auc_hbnb <- function(vifit, data, degs, ndegs, epsilon = 0.1) {
  mu_cond <- hbnbm$extract_vifit(vifit, data, "mu_cond")
  rank_stats <- hbnbm$get_rank_statistics(mu_cond,
    c1 = 1, c2 = 2,
    epsilon = epsilon
  )
  auc <- hbnbm$get_auc(rank_stats, degs, ndegs)
  return(invisible(auc))
}

pd <- init_params_and_data(symsim2be)
saveRDS(object = pd, file = str_glue("symsim2be_{tag}.rds"))
pd <- readRDS(str_glue("symsim2be_{tag}.rds"))

symsim2be_vifit <- hbnbm$run_hbnb_vi(data = pd$data, ip = pd$hip$ip)
hbnb_auc <- get_auc_hbnb(symsim2be_vifit, data = pd$data,
                         symsim_degenes, symsim_ndegs)

message(hbnb_auc$auc_z)
message(hbnb_auc$auc_p)

hbnb_auc_strict <- get_auc_hbnb(symsim2be_vifit, data = pd$data,
                                symsim_degenes, symsim_strict_ndegenes)
message(hbnb_auc_strict$auc_z)
message(hbnb_auc_strict$auc_p)

mycnt <- symsim2be$counts
mybatches <- symsim2be$batch_meta$batch
myconds <- symsim2be$cell_meta$pop
pseudo_deseq2_res <- mypseudo$pseudobulk_deseq2(mycnt, mybatches, myconds)
tmp <- mypseudo$calc_auc(
  pseudo_deseq2_res, symsim_degenes,
  symsim_ndegs
)
message(tmp$auc)
tmp_strict <- mypseudo$calc_auc(
  pseudo_deseq2_res, symsim_degenes,
  symsim_strict_ndegenes
)
message(tmp_strict$auc)
