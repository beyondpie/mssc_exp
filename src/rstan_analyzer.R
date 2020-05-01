import::from(
  rstan, stan, extract, read_rdump, read_stan_csv,
  check_hmc_diagnostics, sampling, stan_rdump, plot,
  get_sampler_params,traceplot
)

model <- "model_i_3"
model <- "model_i_5"
gene <- "HBB"
gene  <- "LYZ"
model_fnm <- paste0("./result/", model, "_", gene, ".csv")

fit <- read_stan_csv(model_fnm)

sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

param <- "mu_g_ic"
param <- "mu_g_di"

## par(mfrow=c(2,2))
plot(fit, pars = c(param))
plot(fit, pars = c(param), plot.fun = "hist")
traceplot(fit, pars=c(param), inc_warmup=TRUE, nrow=2)
## pairs(fit, pars = c(param), las = 1)
