options(error = traceback)
options(warn = -1)

library(Seurat)
library(rstan)
library(ggplot2)
library(bayesplot)
library(data.table)

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(import::from(here::here))
import::from(optparse, make_option, OptionParser, parse_args)
import::from(stringr, str_glue)

options("import.path" = here("rutils"))
myt <- modules::import("transform")

## * load genes information

## * utils for gettign rstan results
load_stan_vi <- function(path) {
  rstan::read_stan_csv(path)
}

load_stan_mc <- function(dirpath, modelnm) {
  csvfiles <- dir(path = dirpath, pattern = paste0(modelnm, "[0-9].csv"),
                  full.names = T)
  rstan::read_stan_csv(csvfiles)
}

## TODO: add explanation about which is case and control
get_ctrlmnscase_par <- function(mystanfit, par="MuCond") {
  mus <- rstan::extract(mystanfit, pars = par)[[par]]
  ## TODO: add colnames of delta: gsymbols
  delta <- as.data.frame(mus[, , 1] - mus[, , 2])
  return(delta)
}

## point relative to regions
mypntrela2rgn <- function(myintvals, myprobs = c(0.025, 0.975), pnt=0.0) {
  q <- as.data.frame(apply(myintvals, 2, quantile, probs=myprobs))
  apply(q, 2, myt$pntrela2rgn, pnt=pnt)
}






