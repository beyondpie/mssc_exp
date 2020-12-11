## mssc 2-1 model

## This inherits from high2
## But model the gene-specific variances of conditions,
## which is actually used in the first-version of mssc.


## * load depenences
library(posterior)
library(cmdstanr)
library(loo)
library(R6)

## * settings
options(error = traceback)
options(warn = 1)
options(mc.cores = 2)

## * common functions

