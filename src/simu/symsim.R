suppressPackageStartupMessages(library(SymSim))
suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)
import::from(stringr, str_glue)
options("import.path" = here("rutils"))
myt <- modules::import("transform")
mysymsim <- modules::import("mysymsim")

## * configs
myseed <- 0L
nbatch <- 10
add_batch_effect <- T
batch_effect_size <- 1
ncell <- 2000
ngene <- 270
hasgenemoudle <- F

## ** cell evf settings
npop <- 2
myphyla <- mysymsim$phyla2()

min_popsize <- 1000
i_minpop <- 1
evf_type <- "discrete"
nevf <- 10
n_de_evf <- 9
sigma <- 0.2

## ** gene modules
minmodn <- 50
gmodprop <- ifelse(hasgenemoudle,
                         minmodn * npop / ngene, 0.0)

## ** utils
my_sim_true <- function(myseed = 0) {
    SimulateTrueCounts(
        ncells_total = ncell, min_popsize = min_popsize,
        i_minpop = i_minpop, ngenes = ngene,
        evf_type = evf_type, nevf = nevf, n_de_evf = n_de_evf,
        phyla = myphyla,
        randseed = myseed,
        gene_module_prop = gmodprop, min_module_size = minmodn,
        Sigma = sigma
    )
}


my_sim_umi <- function(seed) {
    symsimtrue <- my_sim_true(seed)
    symsimumi <- mysymsim$sim_symsim_obs("UMI", symsimtrue)
    suffix <- str_glue("{ncell}_{ngene}_{npop}_{minmodn}_{seed}.rds")
    saveRDS(
        symsimtrue,
        here("src", "simu", "symsimdata", str_glue("true_{suffix}"))
    )
    saveRDS(
        symsimumi,
        here("src", "simu", "symsimdata", str_glue("umi_{suffix}"))
    )
}

## * get repeats of simulations
# lapply(seq_len(10), mysymsim)

## * add individual effects
symsimtrue <- my_sim_true()
symsimumi <- mysymsim$sim_symsim_obs("UMI", symsimtrue)
batchids <- mysymsim$assign_batch_for_cells(symsimumi, nbatch)
