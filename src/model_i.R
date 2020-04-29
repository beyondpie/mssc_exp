import::from(
  rstan, stan, extract, read_rdump, read_stan_csv,
  check_hmc_diagnostics, sampling
)

data <- readRDS("scHBB.Rds")

N <- nrow(data)
I <- ncol(data) - 2 - 2
scale <- 10000
K <- 2
di <- data[, (ncol(data) - 1):ncol(data)]
ic <- data[, 3:(I + 2)]
x_ <- data[[2]]
x_cg <- data[[1]]

modeli_dat <- list(
  N = N, I = I, K = K, scale = scale, di = di,
  ic = ic, x_ = x_, x_cg = x_cg
)

modeli_fit <- stan(
  file = "model_i.stan",
  model_name = "model_i",
  data = modeli_dat,
  chains = 4, warmup = 1000, iter = 2000,
  cores = 2
)
