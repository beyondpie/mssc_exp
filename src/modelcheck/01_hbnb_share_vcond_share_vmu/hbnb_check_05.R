## Check MSSC Hierarhical Bayesian Model

## Instead of using SymSim to simulate the dataset
## Here we firstly use the data simulated from our model.
## Hyper parameters are learned from the real dataset: PBMC.

## * set R environment
source("hbnb_set_r_lib_env_01.R")

hbnbm <- modules::import("hbnb_mssc_03")

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## * load pbmc demo data and vifit.
refdir <- here::here("src", "modelcheck", "01_hbnb_share_vcond_share_vmu")
demo_pbmc <- readRDS(file.path(refdir, "snb_pool_ref_pbmc.rds"))
demo_vifit <- readRDS(file.path(refdir, "est_hbnb_params_ref_pbmc.rds"))

## * functions
simulate_mu_and_r <- function(ngene = 40, mu, r) {
  ## mu: n by g matrix from vifit draw
  ## r: n by g matrix from vifit draw

  ## simulate mu and r: g by 2 matrix

  cols <- sample(seq_len(ncol(mu)),
    size = ngene,
    replace = T
  )
  invisible(t(vapply(seq_len(ngene), function(i) {
    col <- cols[i]
    c(
      sample(mu[, col], size = 1),
      sample(r[, col], size = 1)
    )
  }, FUN.VALUE = c(0.0, 0.0))))
}

simulate_muind <- function(nind = 10, ngene = 40, muind) {
  ## muind: n by g by nind* from vifit draw

  ## simulate muind: g by nind

  gcols <- sample(seq_len(dim(muind)[2]),
    size = ngene,
    replace = T
  )
  invisible(t(vapply(gcols, function(g) {
    icols <- sample(seq_len(dim(muind)[3]),
      size = nind,
      replace = T
    )
    invisible(vapply(icols, function(i) {
      sample(muind[, g, i], size = 1)
    }, FUN.VALUE = 0.0))
  }, FUN.VALUE = rep(0.0, nind))))
}

simulate_mucond <- function(ngene = 40, ncond = 2, mucond,
                            sd_noise = 0.05) {
  ## mucond: n by g by ncond (2)
  ndiff <- ngene / 2
  gcols <- sample(seq_len(dim(mucond)[2]),
    size = ndiff,
    replace = T
  )
  diff <- t(vapply(gcols, function(g) {
    row <- sample(dim(mucond)[1], size = 1)
    invisible(mucond[row, g, ])
  }, FUN.VALUE = c(0.0, 0.0)))

  undiff <- t(vapply(seq_len(ngene - ndiff, function(i) {
    invisible(rnorm(ncond, mean = 0.0, sd = sd_noise))
  }, FUN.VALUE = c(0.0, 0.0))))

  invisible(rbind(diff, undiff))
}

simulate_s <- function(ncell=100, nind = 10, s) {
}

simulate_params <- function(nind = 10, ngene = 40,
                            ncond = 2) {
  ## simulate the params for data generation.
  ## parameters are mainly learned from pbmc dataset.
}
simulate_data <- function() {}

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
