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

## ** mssc hbnb params default settings
## *** default hyper params
dhp <- list(mu0 = rep(0.0, g),
            hp_varofmu = c(1.0, 1.0),
            hp_r = c(1.0, 1.0),
            hp_alpha_varofind = c(1.0, 1.0),
            hp_beta_varofind = c(1.0, 1.0),
            hp_alpha_varofcond = c(1.0, 1.0),
            hp_beta_varofcond = c(1.0, 1.0))
## *** default init params
dip <- list(nb_r = rep(10.0, g),
            varofmu = 25.0,
            raw_mu = rep(0.0, g),
            raw_mu_cond = array(0.0, dim = c(j, g)),
            raw_mu_ind = array(0.0, dim = c(k, g)),
            hp_varofcond = rep(4.0, g),
            varofind = rep(1.0, g))

## * load stan models
snbm <- cmdstan_model(
  here::here("src", "stan", "scale_nb.stan")
)

hbnbm <- cmdstan_model(
  here::here("src", "dirty_stan", "hbnb_rndeff.stan"),
  compile = T, quiet = FALSE
)

## * functions
