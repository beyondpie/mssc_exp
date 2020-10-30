## Check MSSC Hierarhical Bayesian Model

## Use SymSim to simulate the dataset

## * set R environment
import::from(here, here)
suppressPackageStartupMessages(library(tidyverse))
library(MCMCpack)

## develop version of cmdstanr
## devtools::install_github("stan-dev/cmdstanr")
library(cmdstanr)
## use the code below to reload the library
## detach("package:cmdstanr", unload = TRUE)
## library(cmdstanr)

library(grid)
library(gtable)
library(gridExtra)
library(bayesplot)
## for bayesplot plotting
## color_scheme_set("brewer-Spectral")

library(posterior)
library(bbmle)
library(sads)
## install_github("olafmersmann/truncnorm")
library(truncnorm)

## local modules
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
## use the code below to reload the modules
## modules::reload(myfit)
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## warnings/errors traceback settings

## options(error = traceback)
options(warn = 0)
options(mc.cores = 3)

## * configs
## for simulation
num_of_cell_per_ind <- 280
num_of_ind <- 10
num_of_ind_per_cond <- 5
num_of_cond <- 2
sgn <- "NFKB1"

