import::from(
  rstan, stan, extract, read_rdump, read_stan_csv,
  check_hmc_diagnostics, sampling, stan_rdump
)

data <- readRDS("scHBB.Rds")

N <- nrow(data)
I <- ncol(data) - 2 - 2
scale <- 10000
K <- 2
di <- as.matrix(data[, (ncol(data) - 1):ncol(data)])
ic <- as.matrix(data[, 3:(I + 2)])
x_ <- data[[2]]
x_cg <- data[[1]]


stan_rdump(c("N","I","K","scale","di","ic","x_","x_cg"), file="./model_i.rdump")
# In terminal, use stan make to compile and then sample
# stancm `pwd`/model_i.stan
# ./model_i methpd=sample adopt delta=0.9 data file=model_i.rdump

model_i_res  <- read_stan_csv("output.csv")


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


