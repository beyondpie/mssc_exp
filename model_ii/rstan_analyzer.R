import::from(
  rstan, stan, extract, read_rdump, read_stan_csv,
  check_hmc_diagnostics, sampling, stan_rdump, plot,
  get_sampler_params,traceplot
)
model <- "model_ii_1"

# model 2, 4 not that good.
gene <- "HBB"
gene  <- "LYZ"
model_fnm <- paste0("./", model, "_", gene, ".csv")
fit <- read_stan_csv(model_fnm)

sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

param <- "mu_g_ic"
param <- "mu_g_di"

## par(mfrow=c(2,2))
plot(fit, pars = c(param))
plot(fit, pars = c(param), plotfun = "hist")
traceplot(fit, pars=c(param), inc_warmup=TRUE, nrow=2)
## pairs(fit, pars = c(param), las = 1)

dfr <- as.data.frame(fit)

diff  <- dfr[["mu_g_di[1]"]]  - dfr[["mu_g_di[2]"]]
hist(diff, breaks=100)
quantile(diff, 0.025)
quantile(diff, 0.9725)
ecdf(diff)(0.0)

# foldchange, influenced by small values quite a lot.
# fc <- exp(diff)
# hist(fc, breaks=100)

