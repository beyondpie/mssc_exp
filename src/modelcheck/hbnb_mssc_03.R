## hbnb model

## param names in model
## 1) names defined in  stan script (by us)
## 2) names returned by stan sampler/opt
##    (stan adds [i,j] or [i] for matrix/vec)
## 3) names from the estimation directly from the data (only part of params)
##    see fig_mg_sng for details.

## * load R env
source("hbnb_set_r_lib_env_01.R")

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## * load param fitting functions
pf <- modules::import("hbnb_param_fitting_02")

## * load stan model
nm_params <- c("nb_r", "hp_r", "varofmu", "mu",
               "varofcond", "mu_cond",
               "hp_varofind", "varofind", "mu_ind",
               "raw_mu", "raw_mu_cond", "raw_mu_ind")

snbm_for_mur <- cmdstanr::cmdstan_model(
  here::here("src", "stan", "scale_nb.stan"),
  compile = T, quiet = F
)

snbm_for_mucond <- cmdstanr::cmdstan_model(
  here::here("src", "stan", "scale_nb_fixed_r.stan"),
  compile = T, quiet = F)

hbnbm <- cmdstan_model(
  here::here("src", "dirty_stan", "mssc_hbnb_vi.stan"),
  compile = T, quiet = F
)

## * constants
## ** default gamma param
dftgamma <- c(1.0, 1.0)
## default invgamma param
dftinvg <- c(1.0, 10.0)

## ** stan training
num_iter <- 5000
fresh <- 5000
eval_elbo <- 100

## * functions
get_hbnb_param_nms <- function(k, j, g) {

  index_mu_cond <- matrix(rep(seq_len(g),2), nrow = 2, byrow = TRUE)
  index_mu_ind <- matrix(rep(seq_len(g), k), nrow = k, byrow = TRUE)

  hbnb_param_nms <- list(
    hp_r = pf$str_glue_vec(nm = "hp_r", a = seq_len(2)),
    nb_r = pf$str_glue_vec(nm = "nb_r", a = seq_len(g)),
    varofmu = "varofmu",
    raw_mu = pf$str_glue_vec(nm = "raw_mu", a = seq_len(g)),
    mu = pf$str_glue_vec(nm = "mu", a = seq_len(g)),
    varofcond = "varofcond",
    raw_mu_cond = pf$str_glue_mat(nm = "raw_mu_cond", a = index_mu_cond),
    mu_cond = pf$str_glue_mat(nm = "mu_cond", a = index_mu_cond),
    raw_mu_ind = pf$str_glue_mat(nm = "raw_mu_ind", a = index_mu_ind),
    mu_ind = pf$str_glue_mat(nm = "mu_ind", a = index_mu_ind),
    hp_varofind = pf$str_glue_vec(nm = "hp_varofind", a = seq_len(2)),
    varofind = pf$str_glue_vec(nm = "varofind", a = seq_len(g))
  )
  return(invisible(hbnb_param_nms))
}

get_default_hi_params <- function(k, j, g) {
  ## hi: hyper and initial
  ## default hyper params
  dhp <- list(mu0 = rep(0.0, g),
    hp_varofmu = dftinvg,
    hp_alpha_r = dftgamma,
    hp_beta_r = dftgamma,
    hp_alpha_varofind = dftgamma,
    hp_beta_varofind = dftgamma,
    hp_alpha_varofcond = dftgamma,
    hp_beta_varofcond = dftgamma)
  ## default init params
  dip <- list(
    hp_r = dftgamma,
    nb_r = rep(10.0, g),
    varofmu = 25.0,
    mu = rep(0.0, g),
    raw_mu = rep(0.0, g),
    ## TODO: check stan dimension settings
    mu_cond = array(0.0, dim = c(j, g)),
    raw_mu_cond = array(0.0, dim=c(j,g)),
    mu_ind = array(0.0, dim = c(k, g)),
    raw_mu_ind = array(0.0, dim = c(k, g)),
    ## TODO: more reasonable settings.
    hp_varofcond = rep(4.0, g),
    varofind = rep(1.0, g))
  return(invisible(list(hp = dhp, ip = dip)))
}

set_hi_params <- function(k, j, g,
                          cnt, s, cond, ind,
                          scale = 1.96^2
                          ) {

  murnm = c("mu0", "r0")
  mucondnm = pf$str_glue_vec("mu_cond", seq_len(2))
  muindnm = pf$str_glue_vec("mu_ind", seq_len(k))

  dhip <- get_default_hi_params(k, j, g)
  mat <- pf$fit_mg_snb(cnt = cnt, s = s, cond = cond, ind = ind,
                       murnm = murnm, mucondnm = mucondnm,
                       muindnm = muindnm)
  r <- pf$init_hbnb_params(mat, murnm = murnm, mucondnm = mucondnm,
                           muindnm = muindnm, scale = scale)

  ## update the gene scale log mean expression estimation
  dhip$hp$mu0 <- r$hp$mu0

  for (n in nm_params) {
    dhip$ip[[n]]  <- r$init[[n]]
  }
  return(invisible(list(hp = dhip$hp, ip = dhip$ip)))
}

calibrate_init_params <- function(ip, data, seed = 1L) {
  opt <- hbnbm$optimize(
                 data = data,
                 init = list(ip),
                 seed = seed,
                 refresh = refresh,
                 iter = num_iter,
                 algorithm = "lbfgs")
  if (pf$is_vi_or_opt_success()) {
    map <- opt$mle()
    k <- data$k
    j <- data$j
    g <- data$g

    param_nms <- get_hbnb_param_nms(k = k, j = j, g = g)
    ip$raw_mu <- map[param_nms$raw_mu]
    ip$mu <- map[param_nms$mu]

    ip$raw_mu_ind <- matix(map[param_nms$raw_mu_ind], nrow = k, byrow=TRUE)
    ip$mu_ind <- matrix(map[param_nms$mu_ind], nrow = k, byrow = TRUE)

    ip$raw_mu_cond <- matix(map[param_nms$raw_mu_cond], nrow = 2, byrow = TRUE)
    ip$mu_cond <- matrix(map[param_nms$mu_cond], nrow = 2, byrow = TRUE)
  }
  return(invisible(ip))
}

run_hbnb_vi <- function(data, ip, seed = 1L) {
  invisible(
    hbnbm$variational(
            data = data,
            init = list(ip),
            seed = seed,
            refresh = refresh,
            iter = num_iter,
            eval_elbo = eval_elbo,
            adapt_engaged = TRUE
          )
  )
}

