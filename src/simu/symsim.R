library(SymSim)
library(tidyverse)
import::from(here, here)
import::from(stringr, str_glue)
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
nevf <- 10
n_de_evf <- 5
sigma <- 0.15

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
        Sigma = sigma
    )
}

my_sim_obs <- function(protocol, true_data, add_batch_effect = T,
                       nbatch = 1, batch_effect_size = 1) {
    depth_mean <- ifelse(protocol == "UMI",
        depth_mean_umi, depth_mean_fullength
    )
    depth_sd <- ifelse(protocol == "UMI", depth_sd_umi, depth_sd_fullength)
    alpha_mean <- ifelse(protocol == "UMI",
        alpha_mean_umi, alpha_mean_fullength
    )
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


mysymsim <- function(seed) {
    true_counts <- my_sim_true(seed)
    obs_umi <- my_sim_obs("UMI", true_counts,
        add_batch_effect = add_batch_effect,
        nbatch = nbatch, batch_effect_size = batch_effect_size
    )
    suffix <- str_glue("{ncell}_{ngene}_{npopulation}_{minmodulesize}_{seed}.rds")
    saveRDS(
        true_counts,
        here("src", "simu", "symsimdata", str_glue("true_{suffix}"))
    )
    saveRDS(
        obs_umi,
        here("src", "simu", "symsimdata", str_glue("umi_{suffix}"))
    )
}

## * main
lapply(seq_len(10), mysymsim)
