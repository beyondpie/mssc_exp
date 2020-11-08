## Check MSSC Hierarhical Bayesian Model

## Instead of using SymSim to simulate the dataset
## Here we firstly use the data simulated from our model.
## Hyper parameters are learned from the real dataset: PBMC.

## * set R environment
source("hbnb_set_r_lib_env_01.R")

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## * load mssc hbnb script
hbnbmssc <- modules::import("hbnb_mssc_03")

## * configs
## number of conditions
k <- 10
j <- 2

## ** mssc data settings
g <- 200

## * functions
simulate_params <- function(){}
simulate_data <- function(){}

## TODO
get_ground_truth_params <- function(data, varnm) {}

## * run hbnb
default_hi_params <- hbnbmssc$get_default_hi_params(k = k, j = j, g = g)

hi_params <- hbnbmssc$set_hi_params(default_hi_params,
                                    k = k,
                                    j = j,
                                    g = g,
                                    cnt = cnt,
                                    s = s,
                                    cond = cond,
                                    ind = ind,
                                    scale = 1.96^2,
                       murnm = c("mu0", "r0"),
                       mucondnm = str_glue_vec("mu_cond", seq_len(2)),
                       muindnm = str_glue_vec("mu_ind", seq_len(k))
                                    )
init_params <- hbnbmssc$calibrate_init_params(hi_params$ip, data)

hbnbvifit <- hbnbmssc$run_hbnb_vi(data, init_params)

## * analyze results
