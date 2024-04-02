## Simulate data from the model
## - check the model fitting

## model: high, hbnb with indeff across genes
## previous, hide, hbnb with indeff across individuals.

## * set R env
suppressPackageStartupMessages(library(tidyverse))
library(MCMCpack)
library(cmdstanr)
library(grid)
library(gtable)
library(gridExtra)
library(bayesplot)
library(posterior)
library(bbmle)
library(sads)
library(truncnorm)

## warnings/errors traceback settings
options(error = traceback)
options(warn = 1)
options(mc.cores = 3)

highm <- modules::import("mssc_hbnb_indeff_across_gene")
options("import.path" = here::here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")


## * load pbmc demo data and vifit
