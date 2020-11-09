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
      ## sample(mu[, col], size = 1),
      ## sample(r[, col], size = 1)
      mean(mu[, col]),
      min(10.0, median(r[, col]))
    )
  }, FUN.VALUE = c(0.0, 0.0))))
}

simulate_muind <- function(tnind = 10, ngene = 40, muind) {
  ## tnind: total number of individual in the dataset
  ## muind: n by g by tnind from vifit draw

  ## simulate muind: ngene by tnind

  gcols <- sample(seq_len(dim(muind)[2]),
    size = ngene,
    replace = T
  )
  invisible(t(vapply(gcols, FUN = function(g) {
    icols <- sample(seq_len(dim(muind)[3]),
      size = tnind,
      replace = T
    )
    invisible(vapply(icols, FUN = function(i) {
      ## sample(muind[, g, i], size = 1)
      mean(muind[, g, i])
    }, FUN.VALUE = 0.0))
  }, FUN.VALUE = rep(0.0, tnind) ) ))
}

simulate_mucond <- function(ngene = 40, ncond = 2, mu_cond,
                            sd_noise = 0.05) {
  ## mu_cond: n by g by ncond (2)

  ## simulate mu_cond: ngene by ncond
  ndiff <- ngene / 2
  gcols <- sample(seq_len(dim(mu_cond)[2]),
    size = ndiff,
    replace = T
  )
  diff <- t(vapply(gcols, function(g) {
    row <- sample(dim(mu_cond)[1], size = 1)
    invisible(mu_cond[row, g, ])
  }, FUN.VALUE = c(0.0, 0.0)))

  undiff <- t(vapply(seq_len(ngene - ndiff), FUN = function(i) {
    invisible(rnorm(ncond, mean = 0.0, sd = sd_noise))
  }, FUN.VALUE = c(0.0, 0.0)))

  invisible(rbind(diff, undiff))
}

simulate_data <- function(nind = 5, ncell = 100, ngene = 40,
                          s, params_vifit) {
  ## nind is per condition
  ## ncell is per individual
  ## ncond equals to 2
  ncond <- 2
  n <- nind * ncell * ncond
  tnind <- nind * ncond
  s <- sample(s, size = n, replace = T)
  logs <- log(s)
  cond <- rep(seq_len(ncond), each = nind * ncell)
  ind <- rep(seq_len(ncond * nind), each = ncell)

  mur <- simulate_mu_and_r(ngene,
    mu = params_vifit$mu,
    r = params_vifit$nb_r
  )
  mu_ind <- simulate_muind(tnind * ncond,
    ngene = ngene,
    params_vifit$mu_ind
  )
  mu_cond <- simulate_mucond(
    ngene = ngene, ncond = ncond,
    mu_cond = params_vifit$mu_cond,
    sd_noise = 0.05
  )
  params <- list(
    mu = mur[, 1], r = mur[, 2],
    mu_ind = mu_ind, mu_cond = mu_cond
  )
  c2y <- vapply(seq_len(ngene), function(g) {
    loglambda <- logs + mur[g, 1] + mu_cond[g, cond] + mu_ind[g, ind]
    invisible(rnbinom(n = n, mu = exp(loglambda), size = mur[g, 2]))
  }, FUN.VALUE = rep(0.0, n))
  data <- list(s = s, cond = cond, ind = ind, y2c = t(c2y))
  return(invisible(list(params = params, data = data)))
}

## * main
## could be setting as function argument
nind <- 5
ncell <- 100
ngene <- 100
ncond <- 2

## simulate the system
hbnbsim <- simulate_data(nind = nind, ncell = ncell, ngene = ngene,
                        s = demo_pbmc$s, params_vifit = demo_vifit)
## run hbnb

get_vifit <- function(hbnbsim, calibrate_with_opt=FALSE) {
  ind <- hbnbsim$data$ind
  cond <- hbnbsim$data$cond
  s <- hbnbsim$data$s
  y2c <- hbnbsim$data$y2c
  hi_params <- hbnbm$set_hi_params(k = max(ind),
                                   j = max(cond),
                                   g = nrow(y2c),
                                   cnt = y2c,
                                   s = s,
                                   cond = cond, ind = ind,
                                   scale = 1.96^2)
  data <- hbnbm$to_hbnb_data(y2c, ind, cond, s, hi_params$hp)
  if (calibrate_with_opt) {
    ip <- hbnbm$calibrate_init_params(hi_params$ip, data = data)
  } else {
    ip <- hi_params$ip
  }
  vifit <- hbnbm$run_hbnb_vi(data = data, ip = ip)
  est_params <- lapply(hbnbm$nm_params, function(nm) {
    hbnbm$extract_vifit(vifit, data, nm)
  })
  names(est_params) <- hbnbm$nm_params
  return(invisible(list(vifit = vifit, data = data,
                        ip = ip, hip = hi_params,
                        est_params = est_params)))
}

hbnb_vifit <- get_vifit(hbnbsim)

## save intermidiate result
saveRDS(object = list(hbnb_vifit=hbnb_vifit, hbnbsim = hbnbsim),
        file = str_glue("hbnb_vifit_{nind}_{ncell}_{ngene}.rds"))

## run pseudobulk
get_pseudobulk <- function(data) {}
