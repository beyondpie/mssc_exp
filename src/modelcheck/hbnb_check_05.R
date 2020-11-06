## Check MSSC Hierarhical Bayesian Model

## Instead of using SymSim to simulate the dataset
## Here we firstly use the data simulated from our model.
## Hyper parameters are learned from the real dataset: PBMC.

## * set R environment
source("set_r_lib_env.R")

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## * configs
## number of conditions
k <- 10
j <- 2

## ** mssc data settings
g <- 200


## * functions
