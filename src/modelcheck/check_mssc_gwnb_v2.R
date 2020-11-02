## Check MSSC variational inference model.

## Use simulation based calibration.
## - simulation from prior.
## - [hyper] parameters are setted based on a real dataset:
##   - PBMC dataset

## in the second version, we use stan to support the
## negative binomial fitting.

## * set R environment
import::from(here, here)
suppressPackageStartupMessages(library(tidyverse))
library(optparse)
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
opts <- list(
  make_option(c("--ncell", action = "store", type = "integer", default = 200)),
  make_option(c("--nind", action = "store", type = "integer", default = 10)),
  make_option(c("--celltype",
    action = "store", type = "character",
    default = "Naive CD4+ T"
  ))
)

args <- parse_args(OptionParser(option_list = opts))
num_of_cell_per_ind <- args$ncell
num_of_ind_per_cond <- args$nind
num_of_cond <- 2
num_of_ind <- num_of_ind_per_cond * num_of_cond
celltype <- args$celltype
num_top_gene <- 200

## for debug
num_of_cell_per_ind <- 200
num_of_ind_per_cond <- 5
num_of_cond <- 2
num_of_ind <- num_of_cond * num_of_ind_per_cond
celltype <- "Naive CD4+ T"
num_top_gene <- 200
sgn <- "NFKB1"

## * load stan models
## if compiling not work, try rebuild_cmdstan()
snbm <- cmdstan_model(
  here::here("stanutils", "scale_nb.stan"),
  compile = T, quiet = T
)

snb_fixr_m <- cmdstan_model(
  here::here("stanutils", "scale_nb_fixed_r.stan"),
  compile = T, quiet = T
)

mssc_gwnb_muind_model <- cmdstan_model(here::here(
  "src", "dirty_stan",
  "gwnb_simu_from_prior.stan"
), compile = T, quiet = T)

mssc_gwnb_rndeff_model <- cmdstan_model(
  here::here("src", "dirty_stan", "gwnb_rndeffect.stan"),
  compile = T, quiet = T
)

## * load the dataset
pbmc_seurat <- mypbmc$load_pbmc_seurat() %>%
  mypbmc$extract_from_seurat(pbmc_seurat = .)

subscdata <- mypbmc$get_celltype_specific_scdata(
  pbmc_seurat$cnt,
  pbmc_seurat$resp,
  pbmc_seurat$inds,
  pbmc_seurat$ct,
  ## limit to the cell type
  celltype
)

cnt <- subscdata$cnt
inds <- subscdata$inds
resp <- subscdata$resp
sumcnt <- colSums(cnt)

## functions
fit_singlegene_nb <- function(gn = "NFKB1", cnt,
                              inds,
                              resp,
                              sumcnt){
  ## get mu0, mucond, r as in mssc gwnb model
  ## from the real dataset
  ## use stan fit

  d1g_cnt <- cnt[gn, ]
  outliers <- myfit$is_outlier(d1g_cnt)
  d1g_cnt <- d1g_cnt[!outliers]
  d1g_inds <- inds[!outliers]
  d1g_resp <- resp[!outliers]
  d1g_sumcnt <- sumcnt[!outliers]

  ## fit NB dist
  fit_mur <- myfit$stanfit_gwsnb_till_cond_level(
    s = d1g_sumcnt,
    y = d1g_cnt,
    s_control = d1g_sumcnt[d1g_resp == 0],
    y_control = d1g_cnt[d1g_resp == 0],
    s_case = d1g_sumcnt[d1g_resp == 1],
    y_case = d1g_cnt[d1g_resp == 1],
    scale_nb_model = snbm,
    scale_nb_fixed_r_model = snb_fixr_m,
  )
  invisible(list(
    mur = fit_mur,
    cnt = d1g_cnt,
    inds = d1g_inds,
    resp = d1g_resp,
    sumcnt = d1g_sumcnt
  ))
}

set_gwnb_hyper_params <- function(sg_mur, sigmaG0 = 4.0,
                                  default_control_value = 1.0,
                                  default_case_value = -2.0) {
  ## set hyper params based on the fitted sg_mur
  ## from a real dataset

  ## inv-gamma mode estimate.
  ## instead of using var(sg_mur$mu_cond), which might be
  ## small when both of them are in the same sign:
  ## - we assume the mean is around zero

  ## we use the max absolute value, and cover it by 1.5 fold
  mu_control <- sg_mur$mur$mu_cond[1]
  if (is.nan(sg_mur$mur$mu_cond[1])) {
    mu_control <- default_control_value
  }
  mu_case <- sg_mur$mur$mu_cond[2]
  if (is.nan(sg_mur$mur$mu_cond[2])) {
    mu_case <- default_case_value
  }

  alphaTau2G <- 1.0
  betaTau2G <- (alphaTau2G + 1) * max(abs(c(mu_control, mu_case))) * 1.5


  ## when phi2G ~ inv-gamma (alpha, beta)
  ## alphaPhi2G <- 1.0
  ## use mode
  ## betaPhi2G <- min(5.0, (alphaPhi2G + 1) * sg_mur$r0)

  ## when phi2G ~ gamma (alpha, beta)
  betaPhi2G <- 1.0
  ## use mean
  alphaPhi2G <- betaPhi2G * min(10.0, sg_mur$r0)

  invisible(list(
    alphaKappa2G = 1.0,
    betaKappa2G = 1.0,
    alphaTau2G = alphaTau2G,
    betaTau2G = betaTau2G,
    ## control the mean not too small in simulaiton.
    muG0 = max(-7.0, sg_mur$mu0),
    sigmaG0 = sigmaG0,
    alphaPhi2G = alphaPhi2G,
    betaPhi2G = min(5, betaPhi2G)
  ))
}

get_vec_of_repeat_int <- function(to_int, repeatimes, from_int = 1) {
  invisible(c(vapply(
    seq(from_int, to_int),
    FUN = function(i) {
      rep(i, repeatimes)
    }, FUN.VALUE=rep(0L, repeatimes)
  )))
}

simulate_from_gwnb_prior <- function(hp, num_of_cond = 2) {
  ## sampling the parameters from the prior
  ## defined by the hyper parameters.

  ## limit tau size
  tau2g <- min(10.0, rinvgamma(1,
    shape = hp$alphaTau2G,
    scale = hp$betaTau2G
  ))
  ## in reality, mucond should be in different signs.
  mucond <- c(
    truncnorm::rtruncnorm(1, a = -Inf, b = 0.0, mean = 0.0, sd = sqrt(tau2g)),
    truncnorm::rtruncnorm(1, a = 0.0, b = Inf, mean = 0.0, sd = sqrt(tau2g))
  )

  ## mucond <- rnorm(num_of_cond, 0.0, sqrt(tau2g))

  invisible(list(
    mug = max(
      -6.0,
      rnorm(1, hp$muG0, hp$sigmaG0)
    ),
    kappa2g = min(10.0, rinvgamma(1,
      shape = hp$alphaKappa2G,
      scale = hp$betaKappa2G
    )),
    tau2g = tau2g,
    phi2g = min(10.0, rinvgamma(1,
      shape = hp$alphaPhi2G,
      scale = hp$betaPhi2G
    )),
    mucond = mucond
  ))
}

generate_gwnc_y <- function(mug, mucond, kappa2g, phi2g, vec_sumcnt) {
  ## geneerate data based on the gwnb model
  ## given the parameters.
  sample_sumcnt <- sample(vec_sumcnt, size = num_of_cell_per_ind * num_of_ind,
                          replace = TRUE)
  half_n <- num_of_ind_per_cond * num_of_cell_per_ind
  n <- 2 * half_n

  ## cond in simulation start from 1 following stan.
  vec_of_cond <- c(rep(1, half_n), rep(2, half_n))

  index_of_ind_per_cond <- get_vec_of_repeat_int(
    num_of_ind_per_cond,
    num_of_cell_per_ind
  )
  vec_of_ind_under_cond <- c(
    paste0("NR", index_of_ind_per_cond),
    paste0("R", index_of_ind_per_cond)
  )
  vec_of_ind <- get_vec_of_repeat_int(num_of_ind, num_of_cell_per_ind)
  y <- vapply(
    seq_len(n),
    function(i) {
      logmu <- rnorm(1,
        mean = mug + mucond[vec_of_cond[i]],
        sd = sqrt(kappa2g)
      )
      rnbinom(
        1,
        mu = sample_sumcnt[i] * exp(logmu), size = phi2g
      )
    }, 0.0
  )

  invisible(list(
    n = n,
    vec_of_cond = vec_of_cond,
    vec_ind_under_cond = vec_of_ind_under_cond,
    vec_of_ind = vec_of_ind,
    s = sample_sumcnt,
    y = y
  ))
}

set_init_params <- function(simu_data, hp, default_control_value = 1.0,
                            default_case_value = -2.0) {
  sy <- simu_data$y
  ss <- simu_data$s
  resp <- simu_data$vec_of_cond

  ## cond in simulation start from 1 following stan.
  index_of_control  <- (resp == 1)
  index_of_case <- (resp != 1)
  kappa2g <- 1.0
  fit_mur <- myfit$stanfit_gwsnb_till_cond_level(
                     s = ss,
                     y = sy,
                     s_control = ss[index_of_control],
                     y_control = sy[index_of_control],
                     s_case = ss[index_of_case],
                     y_case = sy[index_of_case],
                     scale_nb_model = snbm,
                     scale_nb_model = snb_fixr_m
  )

  if (is.nan(fit_mur$mu_cond[1])) {
    fit_mur$mucond[1] <- default_control_value
  }

  if (is.nan(fit_mur$mu_cond[2])) {
    fit_mur$mucond[2] <- default_case_value
  }

  ## we use the max absolute value, and cover it by 1.5 fold
  tau2g <- max(abs(fit_mur$mu_cond)) * 1.5
  invisible(list(
    MuG = fit_mur$mu0,
    MuIndRaw = rnorm(simu_data$n, 0.0, sqrt(kappa2g)),
    MuCond = fit_mur$mu_cond,
    MuCondRaw = fit_mur$mu_cond / sqrt(tau2g),
    Kappa2G = 1.0,
    Tau2G = tau2g,
    Phi2G = fit_mur$r0
  ))
}

set_gwnb_env <- function(sgn, cnt, inds, resp, sumcnt) {
  ## given a single gene name, and the corresponding real data
  ## - estimate the hyper params
  ## - sampling params from the prior
  ## - generate the data
  ## - init the params

  ## use pbmc to estimate and set the parameters.
  pbmc_sg_fit <- fit_singlegene_nb(
    sgn, cnt,
    inds, resp, sumcnt
  )
  ## par(mfrow = c(1, 2))
  ## hist(pbmc_sg_fit$cnt[pbmc_sg_fit$resp == 0])
  ## hist(pbmc_sg_fit$cnt[pbmc_sg_fit$resp == 1])


  ## simulation from prior
  gwnb_hyperparmas <- set_gwnb_hyper_params(pbmc_sg_fit, sigmaG0 = 4.0)
  params_from_prior <- simulate_from_gwnb_prior(gwnb_hyperparmas)
  gwnb_genrt_data <- generate_gwnc_y(params_from_prior$mug,
    params_from_prior$mucond,
    params_from_prior$kappa2g,
    params_from_prior$phi2g,
    vec_sumcnt = sumcnt
  )

  ## set data for mssc.
  input_data <- c(
    list(
      N = gwnb_genrt_data$n,
      K = num_of_ind,
      J = num_of_cond,
      S = gwnb_genrt_data$s,
      Ind = gwnb_genrt_data$vec_of_ind,
      Cond = gwnb_genrt_data$vec_of_cond,
      y = gwnb_genrt_data$y
    ),
    gwnb_hyperparmas
  )

  ## set initial params
  init_params <- set_init_params(
    gwnb_genrt_data,
    gwnb_hyperparmas
  )
  invisible(list(
    sg_nbfit = pbmc_sg_fit,
    hp = gwnb_hyperparmas,
    params_from_prior = params_from_prior,
    data = input_data,
    init_params = init_params
  ))
}

run_stan_vi <- function(model, model_env,
                        iter = 5000) {
  ## run vi and MAP (maximum a posterior mode) method
  ## in addition, we also run the corresponding methods without
  ## initial values.

  vi <- model$variational(
    data = model_env$data,
    init = list(model_env$init_params),
    seed = 355113,
    refresh = 100,
    iter = iter,
    eval_elbo = 100,
    adapt_engaged = TRUE
  )
  noi_vi <- model$variational(
    data = model_env$data,
    seed = 355113,
    refresh = 100,
    iter = iter,
    eval_elbo = 100,
    adapt_engaged = TRUE
  )

  opt <- model$optimize(
    data = model_env$data,
    init = list(model_env$init_params),
    seed = 355113,
    refresh = 100,
    iter = iter,
    algorithm = "lbfgs"
  )
  noi_opt <- model$optimize(
    data = model_env$data,
    seed = 355113,
    refresh = 100,
    iter = iter,
    algorithm = "lbfgs"
  )
  invisible(list(
    vi = vi,
    noi_vi = noi_vi,
    opt = opt,
    noi_opt = noi_opt
  ))
}

hist_vi_opt <- function(vi, opt, gwnb_env, varnm,
                        lsize = 1.2 * c(1.5, 0.7, 1.2),
                        lalpha = c(0.75, 0.75, 0.75),
                        lcolor = c("red", "gray20", "blue"),
                        ltype = c(1, 2, 2), show_opt = FALSE) {
  if (varnm == "MuG") {
    value_from_prior <- gwnb_env$params_from_prior$mug
    init_value <- gwnb_env$init_params[[varnm]]
  }
  if (varnm == "MuCond[1]") {
    value_from_prior <- gwnb_env$params_from_prior$mucond[1]
    init_value <- gwnb_env$init_params[["MuCond"]][1]
  }
  if (varnm == "MuCond[2]") {
    value_from_prior <- gwnb_env$params_from_prior$mucond[2]
    init_value <- gwnb_env$init_params[["MuCond"]][2]
  }
  if (varnm == "Kappa2G") {
    value_from_prior <- gwnb_env$params_from_prior$kappa2g
    init_value <- gwnb_env$init_params[[varnm]]
  }
  if (varnm == "Tau2G") {
    value_from_prior <- gwnb_env$params_from_prior$tau2g
    init_value <- gwnb_env$init_params[[varnm]]
  }
  if (varnm == "Phi2G") {
    value_from_prior <- gwnb_env$params_from_prior$phi2g
    init_value <- gwnb_env$init_params[[varnm]]
  }

  ## add hist graph with two lines:
  ## simulation truth and designed init params.
  p <- bayesplot::mcmc_hist(vi$draws(varnm)) +
    bayesplot::vline_at(value_from_prior,
      size = lsize[1], alpha = lalpha[1],
      color = lcolor[1], linetype = ltype[1]
    ) +
    bayesplot::vline_at(init_value,
      size = lsize[2], alpha = lalpha[2],
      color = lcolor[2], linetype = ltype[2]
    )

  if (show_opt) {
    p <- p +
      bayesplot::vline_at(opt$mle(varnm),
        size = lsize[3],
        alpha = lalpha[3], color = lcolor[3], linetype = ltype[3]
      )
  }
  p <- p + bayesplot::facet_text(size = 15, hjust = 0.5)
  invisible(p)
}

visualize_vi_init <- function(vi, gwnb_env,
                              varnms = c(
                                "MuG", "MuCond[1]", "MuCond[2]",
                                "Kappa2G", "Tau2G", "Phi2G"
                              ),
                              init_show_opt = TRUE) {
  l <- lapply(varnms, function(varnm) {
    p <- bayesplot::bayesplot_grid(
      plots = list(
        without_init_vi = hist_vi_opt(vi$noi_vi, vi$noi_opt, gwnb_env,
          varnm,
          show_opt = init_show_opt
        ),
        with_init_vi = hist_vi_opt(vi$vi, vi$opt, gwnb_env,
          varnm,
          show_opt = init_show_opt
        )
      ),
      grid_args = list(nrow = 1)
    )
    invisible(p)
  })
  invisible(l)
}

vis_compare_vi <- function(vi1, vi2, gwnb_env, vi_nms,
                           varnms = c(
                             "MuG", "MuCond[1]", "MuCond[2]",
                             "Kappa2G", "Tau2G", "Phi2G"
                           )) {
  l <- lapply(varnms, function(varnm) {
    plots <- list(
      vi1 = hist_vi_opt(vi1$vi, vi1$opt, gwnb_env,
        varnm,
        show_opt = FALSE
      ),
      vi2 = hist_vi_opt(vi2$vi, vi2$opt, gwnb_env,
        varnm,
        show_opt = FALSE
      )
    )
    names(plots) <- vi_nms

    p <- bayesplot::bayesplot_grid(
      plots = plots,
      grid_args = list(nrow = 1)
    )
    invisible(p)
  })
  invisible(l)
}

run_and_view_variational <- function(sgn = "NFKB1", cnt = cnt, inds = inds,
                                     resp = resp, sumcnt = sumcnt,
                                     tag = 1) {

  ## set hyper param, data, and init values for gwnb mssc
  ## may also have failed fitting
  gwnb_env <- set_gwnb_env(
    sgn = sgn, cnt = cnt,
    inds = inds,
    resp = resp,
    sumcnt = sumcnt
  )

  ## run vi
  muind_vi <- run_stan_vi(mssc_gwnb_muind_model, gwnb_env)
  rndeff_vi <- run_stan_vi(mssc_gwnb_rndeff_model, gwnb_env)

  ## summary the results of vi
  ## initilization plays important roles

  ## init with showing opt, model muind.
  visualize_muind_vi_init <- do.call(
    rbind,
    visualize_vi_init(muind_vi, gwnb_env,
      init_show_opt = TRUE
    )
  )

  comp_init_col_annot <- paste0(
    "Left: without designed init params. ",
    "Right: with designed init params."
  )
  comp_vi_col_annot <- paste0(
    "Left: model individual effect on mean. ",
    "Right: model inidividual effect as variance."
  )
  lines_annot <- paste0(
    "Red line: simulation truth. ",
    "Gray line: designed init params. ",
    "Blue line: MAP estimation."
  )

  p_muind_vi_init <- gridExtra::grid.arrange(visualize_muind_vi_init,
    top = grid::textGrob(paste0("MSSC: model individual mean"),
      gp = grid::gpar(fontsize = 20, font = 3)
    ),
    bottom = grid::textGrob(paste(comp_init_col_annot,
      lines_annot,
      sep = "\n"
    ),
    gp = grid::gpar(fontsize = 15)
    )
  )

  ## init without showing opt, rand eff
  visualize_rndeff_vi_init_nopt <- do.call(
    rbind,
    visualize_vi_init(rndeff_vi, gwnb_env, init_show_opt = FALSE)
  )
  p_rndeff_vi_init_nopt <- gridExtra::grid.arrange(
    visualize_rndeff_vi_init_nopt,
    top = grid::textGrob(paste0("MSSC: model rand effect"),
      gp = grid::gpar(fontsize = 20, font = 3)
    ),
    bottom = grid::textGrob(paste(comp_init_col_annot,
      lines_annot,
      sep = "\n"
    ),
    gp = grid::gpar(fontsize = 15)
    )
  )

  ## compare different vi models
  ## muind: mle estimation is more robust

  visualize_comp_vi <- do.call(
    rbind,
    vis_compare_vi(
      muind_vi, rndeff_vi, gwnb_env,
      c("Model individual mean", "Random effect model")
    )
  )

  p_compare_vi <- gridExtra::grid.arrange(visualize_comp_vi,
    top = grid::textGrob(paste0("MSSC"),
      gp = grid::gpar(fontsize = 20, font = 3)
    ),
    bottom = grid::textGrob(paste(
      comp_vi_col_annot,
      lines_annot,
      sep = "\n"
    ), gp = grid::gpar(fontsize = 15))
  )

  outdir <- here::here(
    "src", "modelcheck",
    stringr::str_glue("gwnb_figures")
  )
  dir.create(outdir, showWarnings = FALSE)
  pdf(stringr::str_glue("{outdir}/{tag}_{nind}_{ncell}.pdf"),
    width = 12, height = 10
  )

  grid::grid.newpage()
  grid::grid.draw(p_muind_vi_init)
  grid::grid.newpage()
  grid::grid.draw(p_rndeff_vi_init_nopt)
  grid::grid.newpage()
  grid::grid.draw(p_compare_vi)
  dev.off()
}

## Classical genes as DE
##   - SNHG16, OASL, NAMPT, NFKB1, BCL2L11, TRAF4, ICAM1, XCL2, XCL1, CCL3L3, CCL3L1

## Strong individual effect genes
##   - HBA1, HBA2, HBD

## possible non-DE
##   - TOX, YIPF5, CCL3, KDM6A, HDDC2


lapply(seq_len(5), function(i) {
  result <- tryCatch(
    {
      run_and_view_variational(
        sgn = "NFKB1", cnt = cnt, inds = inds,
        resp = resp, sumcnt = sumcnt,
        tag = i
      )
      0L
    },
    error = function(e) {
      message(e)
      1L
    }
  )
})
