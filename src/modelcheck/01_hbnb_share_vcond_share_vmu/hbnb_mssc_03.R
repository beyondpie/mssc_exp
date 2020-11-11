## hbnb model

## param names in model
## 1) names defined in  stan script (by us)
## 2) names returned by stan sampler/opt
##    (stan adds [i,j] or [i] for matrix/vec)
## 3) names from the estimation directly from the data (only part of params)
##    see fig_mg_sng for details.

## hbnb row-major order considers the situation:
## [g1,1] [g1,2], ..., [g1,second_dim], [g2,1], ...
## - this is hbnb param nms order
## - cmdstan will named the varibles column-major order
##   like: [g1,1], [g2,1], [g3, 1], ..., [g1, 2], [g2, 2], ...
## As long as we use variable names to extract the parameters,
## it's fine.


## * load R env
source("hbnb_set_r_lib_env_01.R")

options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## * load param fitting functions
## modules::reload(pf)
pf <- modules::import("hbnb_param_fitting_02")

## * load stan model
nm_params <- c(
  "nb_r", "hp_r", "varofmu", "mu",
  "varofcond", "mu_cond",
  "hp_varofind", "varofind", "mu_ind",
  "raw_mu", "raw_mu_cond", "raw_mu_ind"
)

snbm_for_mur <- cmdstanr::cmdstan_model(
  here::here("stanutils", "scale_nb.stan"),
  compile = T, quiet = F
)

snbm_for_mucond <- cmdstanr::cmdstan_model(
  here::here("stanutils", "scale_nb_fixed_r.stan"),
  compile = T, quiet = F
)

## hbnbm <- cmdstan_model(
##   here::here("src", "dirty_stan", "mssc_hbnb.stan"),
##   compile = T, quiet = F
## )

hbnbm <- cmdstan_model(
  here::here("src", "dirty_stan", "mssc_hbnb_mat.stan"),
  compile = T, quiet = F
)

## * constants
## ** stan training
num_iter <- 5000
## set refresh as zero to reduce the output log from stan
## https://github.com/stan-dev/cmdstanr/issues/341
refresh <- 0
eval_elbo <- 100
## vi_algorithm <- "fullrank"
vi_algorithm <- "meanfield"
output_samples <- 1500

## * functions
get_hbnb_param_nms <- function(k, j, g) {
  hbnb_param_nms <- list(
    hp_r = pf$str_glue_vec(nm = "hp_r", a = seq_len(2)),
    nb_r = pf$str_glue_vec(nm = "nb_r", a = seq_len(g)),
    varofmu = "varofmu",
    raw_mu = pf$str_glue_vec(nm = "raw_mu", a = seq_len(g)),
    mu = pf$str_glue_vec(nm = "mu", a = seq_len(g)),
    varofcond = "varofcond",
    raw_mu_cond = pf$str_glue_mat(nm = "raw_mu_cond", nr = g, nc = j),
    mu_cond = pf$str_glue_mat(nm = "mu_cond", nr = g, nc = j),
    raw_mu_ind = pf$str_glue_mat(nm = "raw_mu_ind", nr = g, nc = k),
    mu_ind = pf$str_glue_mat(nm = "mu_ind", nr = g, nc = k),
    hp_varofind = pf$str_glue_vec(nm = "hp_varofind", a = seq_len(2)),
    varofind = pf$str_glue_vec(nm = "varofind", a = seq_len(g))
  )
  return(invisible(hbnb_param_nms))
}

get_default_hi_params <- function(k, j, g) {
  ## hi: hyper and initial
  ## default hyper params
  dhp <- list(
    mu0 = rep(0.0, g),
    hp_varofmu = pf$hpinvg_default,
    hp_alpha_r = pf$hpgamma_default,
    hp_beta_r = pf$hpgamma_default,
    hp_alpha_varofind = pf$hpgamma_default,
    hp_beta_varofind = pf$hpgamma_default,
    hp_varofcond = pf$hpinvg_default
  )
  ## default init params
  dip <- list(
    hp_r = pf$hpgamma_default,
    nb_r = rep(pf$r_default, g),
    varofmu = pf$varofmu_default,
    mu = rep(0.0, g),
    raw_mu = rep(0.0, g),
    mu_cond = array(0.0, dim = c(g, j)),
    raw_mu_cond = array(0.0, dim = c(g, j)),
    varofcond = pf$varofcond_default,
    mu_ind = array(0.0, dim = c(g, k)),
    raw_mu_ind = array(0.0, dim = c(g, k)),
    varofind = rep(1.0, g)
  )
  return(invisible(list(hp = dhp, ip = dip)))
}

set_hi_params <- function(k, j, g,
                          cnt, s, cond, ind,
                          scale = 1.96^2) {
  murnm <- c("mu0", "r0")
  mucondnm <- pf$str_glue_vec("mu_cond", seq_len(2))
  muindnm <- pf$str_glue_vec("mu_ind", seq_len(k))

  dhip <- get_default_hi_params(k, j, g)
  mat <- pf$fit_mg_snb(
    cnt = cnt, s = s, cond = cond, ind = ind,
    snbm = snbm_for_mur,
    snbm_for_mucond = snbm_for_mucond,
    murnm = murnm, mucondnm = mucondnm,
    muindnm = muindnm
  )
  r <- pf$init_hbnb_params(mat,
    murnm = murnm, mucondnm = mucondnm,
    muindnm = muindnm, scale = scale
  )

  ## update the gene scale log mean expression estimation
  dhip$hp$mu0 <- r$hp$mu0

  for (n in nm_params) {
    dhip$ip[[n]] <- r$init[[n]]
  }
  return(invisible(list(hp = dhip$hp, ip = dhip$ip)))
}

to_hbnb_data <- function(cnt, ind, cond, s, hp) {
  ## given the basic data, we translate it into
  ## what hbnb needs.

  return(invisible(c(list(
    n = ncol(cnt),
    k = max(ind),
    j = max(cond),
    g = nrow(cnt),
    s = s,
    cond = cond,
    ind = ind,
    y = t(cnt)
  ), hp)))
}

mu_transform_from_raw <- function(raw_mu, mu0, varofmu) {
  return(invisible(raw_mu * sqrt(varofmu) + mu0))
}

mucond_transfrom_from_raw <- function(raw_mu_cond, varofcond) {
  return(invisible(raw_mu_cond * sqrt(varofcond)))
}

muind_transform_from_raw <- function(raw_mu_ind, varofind) {
  return(invisible(diag(sqrt(varofind)) %*% raw_mu_ind))
}

vi_mu_transform_from_raw <- function(vi_raw_mu, mu0,
                                     vi_varofmu, genenms = NULL) {
  ## get mu: n by g

  ## vi_raw_* are got by vi_fit$draw()

  ## element-wise multiply, when second is column
  ## it will column by column.
  res <- vi_raw_mu * sqrt(as.numeric(vi_varofmu))

  ## n <- nrow(vi_raw_mu)
  ## d <- ncol(vi_raw_mu)
  ## t_res <- vapply(seq_len(n), function(i) {
  ##   mu_transform_from_raw(vi_raw_mu[i, ], mu0, vi_varofmu[i, ])
  ## }, FUN.VALUE = rep(0.0, d))
  ## res <- t(t_res)
  if (!is.null(genenms)) {
    colnames(res) <- genenms
  }
  return(invisible(res))
}

split_matrix_col <- function(mat, second_dim, second_dim_nms = NULL) {
  ## split matrix col into a matrix row-major order

  t <- vapply(seq_len(nrow(mat)), function(i) {
    matrix(mat[i, ], nrow = second_dim, byrow = TRUE)
  }, FUN.VALUE = matrix(mat[1, ], nrow = second_dim, byrow = TRUE))
  r <- aperm(t, c(3, 1, 2))
  if (!is.null(second_dim_nms)) {
    dimnames(r)[[2]] <- second_dim_nms
  }
  return(invisible(r))
}

vi_mucond_transform_from_raw <- function(vi_raw_mu_cond, vi_varofcond,
                                         g, genenms = NULL) {
  ## get mucond: n by g by j (num_of_cond, default is 2)

  ## vi_raw_mu_cond shape: n by g*j
  ## order of the names of vi_raw_mu_cond: by row

  ## vi_varofcond shape: n by 1 if we use cmdstanr draw matrix
  ## then let it be a vector length of n
  t <- vi_raw_mu_cond * sqrt(as.numeric(vi_varofcond))
  r <- split_matrix_col(t, second_dim = g, second_dim_nms = genenms)
  return(invisible(r))
}

vi_muind_transform_from_raw <- function(vi_raw_mu_ind, vi_varofind, g, k,
                                        genenms = NULL) {
  ## vi_raw_mu_ind: n by g * k
  ## vi_varofind: n by g

  ## repeat_var_per_gene: n by g * k:
  ## [g1, 1]. [g1, 2], ..., [g1,k], [g2, 1] ...
  repeat_var_per_gene <- t(apply(vi_varofind, 1, rep, each = k))
  ## element-wise multiplification
  t <- vi_raw_mu_ind * sqrt(repeat_var_per_gene)
  r <- split_matrix_col(t, second_dim = g, second_dim_nms = genenms)
  return(r)
}

calibrate_init_params <- function(ip, data, scale = 1.96^2, seed = 1L) {
  ## use optimization to update the mean related parameters
  ## only raw_mu, raw_mu_ind, raw_mu_cond will be updated.

  ## R will not modify the init_params (ip) even it's a list.

  opt <- hbnbm$optimize(
    data = data,
    init = list(ip),
    seed = seed,
    refresh = refresh,
    iter = num_iter,
    algorithm = "lbfgs"
  )
  if (pf$is_vi_or_opt_success(opt)) {
    map <- opt$mle()
    k <- data$k
    j <- data$j
    g <- data$g

    param_nms <- get_hbnb_param_nms(k = k, j = j, g = g)

    ip$raw_mu <- map[param_nms$raw_mu]
    ip$mu <- mu_transform_from_raw(
      ip$raw_mu, data$mu0,
      map[param_nms$varofmu]
    )

    ip$raw_mu_ind <- matrix(map[param_nms$raw_mu_ind], nrow = g, byrow = TRUE)
    ip$mu_ind <- muind_transform_from_raw(
      ip$raw_mu_ind,
      map[param_nms$varofind]
    )
    r3 <- pf$est_muind(ip$mu_ind, scale)
    ip$varofind <- r3$varofind
    ip$hp_varofind <- r3$hp_varofind

    ip$raw_mu_cond <- matrix(map[param_nms$raw_mu_cond], nrow = g, byrow = TRUE)
    ip$mu_cond <- mucond_transfrom_from_raw(
      ip$raw_mu_cond,
      map[param_nms$varofcond]
    )
    ip$varofcond <- pf$est_mucond(ip$mu_cond, scale)
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
      adapt_engaged = TRUE,
      algorithm = vi_algorithm,
      output_samples = output_samples
    )
  )
}

extract_vifit <- function(vifit, data, param) {
  ## get the draw matrix for the param
  ## if not exist, return NaN

  if (!(param %in% nm_params)) {
    message(stringr::str_glue("{param} is not recognized."))
    return(NaN)
  }
  if (param %in% c("varofmu", "varofcond")) {
    return(invisible(vifit$draws(param)))
  }
  if (param %in% c("hp_r", "hp_varofind")) {
    return(invisible(vifit$draws(pf$str_glue_vec(param, seq_len(2)))))
  }

  gbc <- t(data$y)
  if (!is.null(rownames(gbc))) {
    genenms <- rownames(gbc)
  } else {
    genenms <- seq_len(data$g)
  }

  if (param %in% c("nb_r", "varofind")) {
    r <- vifit$draws(pf$str_glue_vec(param, seq_len(data$g)))
    colnames(r) <- genenms
    return(invisible(r))
  }

  if (param %in% c("raw_mu", "mu")) {
    r <- vifit$draws(pf$str_glue_vec("raw_mu", seq_len(data$g)))
    if (param == "mu") {
      varofmu <- vifit$draws("varofmu")
      mu0 <- data$mu0
      return(invisible(vi_mu_transform_from_raw(r, mu0, varofmu, genenms)))
    }
    colnames(r) <- genenms
    return(invisible(r))
  }

  if (param %in% c("raw_mu_cond", "mu_cond")) {
    t1 <- vifit$draws(pf$str_glue_mat("raw_mu_cond", nr = data$g, nc = data$j))
    if (param == "mu_cond") {
      t2 <- vifit$draws("varofcond")
      r <- vi_mucond_transform_from_raw(t1, t2, data$g, genenms)
      return(invisible(r))
    }
    return(invisible(split_matrix_col(t1, data$g, genenms)))
  }

  if (param %in% c("raw_mu_ind", "mu_ind")) {
    t1 <- vifit$draws(pf$str_glue_mat("raw_mu_ind", nr = data$g, nc = data$k))
    if (param == "mu_ind") {
      t2 <- vifit$draws(pf$str_glue_vec("varofind", seq_len(data$g)))
      r <- vi_muind_transform_from_raw(t1, t2,
        g = data$g, k = data$k,
        genenms
      )
      return(invisible(r))
    }
    return(invisible(split_matrix_col(t1, data$g, genenms)))
  }
  message(stringr::str_glue("{param} is missed."))
  return(NaN)
}

get_rank_statistics <- function(mu_cond, c1 = 1, c2 = 2,
                                epsilon = 0.1) {
  ## mu_cond: n by ngene by ncond
  ## c1, c2 correspond to different conditions

  delta <- as.matrix(abs(mu_cond[, , c1] - mu_cond[, , c2]))
  if (!is.null(dimnames(mu_cond)[[2]])) {
    colnames(delta) <- dimnames(mu_cond)[[2]]
  }
  abs_delta <- abs(delta)
  ## z score
  z <- colMeans(abs_delta) / matrixStats::colSds(abs_delta)
  ## probability larger than a given epision
  p <- colSums(abs_delta >= epsilon) / nrow(delta)
  return(invisible(list(z = z, p = p, delta = delta)))
}

get_auc <- function(rank_stats, diffg, ndiffg) {
  true_class <- c(rep(TRUE, length(diffg)), rep(FALSE, length(ndiffg)))
  z <- rank_stats$z[c(diffg, ndiffg)]
  p <- rank_stats$p[c(diffg, ndiffg)]
  auc_z <- caTools::colAUC(z, true_class)
  auc_p <- caTools::colAUC(p, true_class)
  return(invisible(list(auc_z = auc_z, auc_p = auc_p)))
}

## * test hbnb_mssc here
test <- function() {
  pbmc <- readRDS(here::here(
    "src", "modelcheck",
    "snb_pool_ref_pbmc.rds"
  ))
  k <- max(pbmc$ind)
  j <- 2
  g <- nrow(pbmc$y2c)

  hi_params <- set_hi_params(
    k = k, j = j, g = g,
    cnt = pbmc$y2c, s = pbmc$s, cond = pbmc$cond,
    ind = pbmc$ind, scale = 1.96^2
  )

  data <- to_hbnb_data(pbmc$y2c,
    ind = pbmc$ind, cond = pbmc$cond,
    s = pbmc$s, hp = hi_params$hp
  )

  ## debug optimization
  ip <- hi_params$ip
  scale <- 1.96^2
  opt <- hbnbm$optimize(
    data = data,
    init = list(ip),
    seed = 1,
    refresh = refresh,
    iter = num_iter,
    algorithm = "lbfgs"
  )

  map <- opt$mle()
  k <- data$k
  j <- data$j
  g <- data$g

  param_nms <- get_hbnb_param_nms(k = k, j = j, g = g)

  ip$raw_mu <- map[param_nms$raw_mu]
  ip$mu <- mu_transform_from_raw(
    ip$raw_mu, data$mu0,
    map[param_nms$varofmu]
  )

  ip$raw_mu_ind <- matrix(map[param_nms$raw_mu_ind], nrow = g, byrow = TRUE)
  ip$mu_ind <- muind_transform_from_raw(
    ip$raw_mu_ind,
    map[param_nms$varofind]
  )
  r3 <- pf$est_muind(ip$mu_ind, scale)
  ip$varofind <- r3$varofind
  ip$hp_varofind <- r3$hp_varofind

  ip$raw_mu_cond <- matrix(map[param_nms$raw_mu_cond], nrow = g, byrow = TRUE)
  ip$mu_cond <- mucond_transfrom_from_raw(
    ip$raw_mu_cond,
    map[param_nms$varofcond]
  )
  ip$varofcond <- pf$est_mucond(ip$mu_cond, scale)

  ## use simpler ip
  vi_sampler <- run_hbnb_vi(data = data, ip = hi_params$ip)
  ## use ip calibrated by optimization
  ip <- calibrate_init_params(hi_params$ip, data = data)
  vi_sampler <- run_hbnb_vi(data = data, ip = ip)

  est_params <- lapply(nm_params, function(nm) {
    extract_vifit(vi_sampler, data, nm)
  })
  names(est_params) <- nm_params
  return(invisible(est_params))
}
