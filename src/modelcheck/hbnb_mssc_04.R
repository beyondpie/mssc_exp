## hbnb model

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

nm_params <- c("nb_r", "varofmu", "mu",
               "varofcond", "mu_cond",
               "hp_varofind", "varofind", "mu_ind")

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


## * functions
get_default_hi_params <- function(k, j, g) {
  ## hi: hyper and initial
  ## default hyper params
  dhp <- list(mu0 = rep(0.0, g),
    hp_varofmu = c(1.0, 1.0),
    hp_r = c(1.0, 1.0),
    hp_alpha_varofind = c(1.0, 1.0),
    hp_beta_varofind = c(1.0, 1.0),
    hp_alpha_varofcond = c(1.0, 1.0),
    hp_beta_varofcond = c(1.0, 1.0))
  ## default init params
  dip <- list(nb_r = rep(10.0, g),
    varofmu = 25.0,
    mu = rep(0.0, g),
    mu_cond = array(0.0, dim = c(j, g)),
    mu_ind = array(0.0, dim = c(k, g)),
    hp_varofcond = rep(4.0, g),
    varofind = rep(1.0, g))
  return(invisible(list(hp = dhp, ip = dip)))
}

set_hi_params <- function(hip, k, j, g,
                          cnt, s, cond, ind,
                          scale = 1.96^2,
                       murnm = c("mu0", "r0"),
                       mucondnm = str_glue_vec("mu_cond", seq_len(2)),
                       muindnm = str_glue_vec("mu_ind", seq_len(k))
                       ) {
  mat <- pf$fit_mg_snb(cnt = cnt, s = s, cond = cond, ind = ind,
                       murnm = murnm, mucondnm = mucondnm,
                       muindnm = muindnm)
  r <- pf$init_hbnb_params(mat, murnm = murnm, mucondnm = mucondnm,
                           muindnm = muindnm, scale = scale)

  hp <- hip$hp
  hp$mu0 <- r$hp$mu0

  ip <- hip$ip
  for (n in nm_params) {
    ip[[n]]  <- r$init[[n]]
  }
  return(invisible(list(hp = hp, ip = ip)))
}

calibrate_init_params <- function(ip, data, seed = 1L) {
  opt <- hbnbm$optimize(
                 data = data,
                 init = list(ip),
                 seed = seed,
                 refresh = 5000,
                 iter = 5000,
                 algorithm = "lbfgs")
  if (pf$is_vi_or_opt_success()) {
    map <- opt$mle()
    ## TODO
    update(ip, map)
  }
  return(invisible(ip))
}

run_hbnb_vi <- function(data, ip, seed = 1L) {
  invisible(
    hbnbm$variational(
            data = data,
            init = list(ip),
            seed = seed,
            ## TODO: set common config
            refresh = 5000,
            iter = 5000,
            eval_elbo = 100,
            adapt_engaged = TRUE
          )
  )
}

