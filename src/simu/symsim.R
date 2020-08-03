library(SymSim)
library(tidyverse)
import::from(here, here)
options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * configs
myseed <- 0L
nbatch <- 10
add_batch_effect <- T
batch_effect_size <- 1
ncell <- 5000
ngene <- 270

## ** cell evf settings
npopulation <- 5
myphyla <- Phyla5()
min_popsize <- 500
i_minpop <- 1
evf_type <- "discrete"
nevf <- 30
n_de_evf <- 18
## ** gene modules
minmodulesize <- 50
genemoduleprop <- minmodulesize * npopulation / ngene

## ** get gene length
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngene, replace = FALSE)

## ** UMI settings
depth_mean_umi <- 45000
depth_sd_umi <- 4500
alpha_mean_umi <- 0.1

## ** fullength settings
depth_mean_fullength <- 1e+05
depth_sd_fullength <- 10000
alpha_mean_fullength <- 0.4

## ** utils
my_sim_true <- function(myseed = 0) {
    SimulateTrueCounts(
        ncells_total = ncell, min_popsize = min_popsize,
        i_minpop = i_minpop, ngenes = ngene,
        evf_type = evf_type, nevf = nevf, n_de_evf = n_de_evf,
        phyla = myphyla,
        randseed = myseed,
        gene_module_prop = genemoduleprop, min_module_size = minmodulesize,
        Sigma = 0.2
    )
}

my_sim_obs <- function(protocol, true_data, add_batch_effect = T,
                       nbatch = 1, batch_effect_size = 1) {
    depth_mean <- ifelse(protocol == "UMI", depth_mean_umi, depth_mean_fullength)
    depth_sd <- ifelse(protocol == "UMI", depth_sd_umi, depth_sd_fullength)
    alpha_mean <- ifelse(protocol == "UMI", alpha_mean_umi, alpha_mean_fullength)
    tmp <- True2ObservedCounts(
        true_counts = true_data[[1]],
        meta_cell = true_data[3], protocol = protocol, alpha_mean = alpha_mean,
        alpha_sd = 0.02, gene_len = gene_len, depth_mean = depth_mean,
        depth_sd = depth_sd, nPCR1 = 14
    )
    if (add_batch_effect) {
        tmp <- DivideBatches(
            observed_counts_res = tmp, nbatch = nbatch,
            batch_effect_size = batch_effect_size
        )
    }
    intc <- tmp
    intc[[1]] <- apply(tmp[[1]], c(1, 2), function(x) {
        ifelse(x > 0, as.integer(x + 1), 0L)
    })
    return(intc)
}


## * main
symsim_true <- my_sim_true(myseed)
symsim_umi <- my_sim_obs("UMI", symsim_true,
    nbatch = nbatch,
    batch_effect_size = batch_effect_size
    )
## nonUMI, i.e., fullength, here nonUMI is needed by SymSim
symsim_fullen <- my_sim_obs("nonUMI", symsim_true,
    nbatch = nbatch,
    batch_effect_size = batch_effect_size
)

## ** save data for plotting
saveRDS(symsim_true, "symsim_true.rds")
saveRDS(symsim_umi, "symsim_umi.rds")
saveRDS(symsim_fullen, "symsim_fullen.rds")

symsim_ptsne <- function(protocol, obs_data, label = "cell_meta.pop") {
    PlotTsne(
        meta = obs_data[[2]], data = log2(obs_data[[1]] +
            1), evf_type = "discrete", n_pc = 20, label = label,
        saving = F, plotname = protocol
    )
}

mysymsim <- function(seed) {
    ## * scRNAseq data simulation ** get true expressions
    true_counts <- my_sim_true(seed)

    ## ** get observed expressions
    obs_umi <- my_sim_obs("UMI", true_counts,
        add_batch_effect = add_batch_effect,
        nbatch = nbatch, batch_effect_size = batch_effect_size
    )
    obs_fullength <- my_sim_obs("fullength", true_counts,
        add_batch_effect = add_batch_effect,
        nbatch = nbatch, batch_effect_size = batch_effect_size
    )
    save(true_counts,
        obs_umi, obs_fullength,
        file = stringr::str_c("symsim_data/sim",
            ncell, ngene, seed, ".RData",
            sep = "_"
        )
    )
}

## main <- function(myrep = 10) {
##     for (i in 1:myrep) {
##         mysymsim(i)
##     }
## }
