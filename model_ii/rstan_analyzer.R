import::from(
  rstan, stan, extract, read_rdump, read_stan_csv,
  check_hmc_diagnostics, sampling, stan_rdump, plot,
  get_sampler_params, traceplot
)

draw_invgamma <- function(alpha, rate = 1 / alpha, from = 0.0001, to = 20) {
  library(invgamma)
  library(ggplot2)
  theme_set(theme_bw())
  x <- seq(from, to, 0.01)
  qplot(x, dinvgamma(x, alpha, rate), geom = "line")
}

plot_posterior <- function(modelnm, genenm, param = c("mu_g_ic", "mu", "mu_g_di")) {
  library(ggplot2)
  library(ggpubr)

  model_fnm <- paste0("./", modelnm, "_", genenm, ".csv")
  fit <- read_stan_csv(model_fnm)
  p1 <- plot(fit, pars = c(param))

  dfr <- as.data.frame(fit)
  diff_di <- data.frame(diff = dfr[["mu_g_di[1]"]] - dfr[["mu_g_di[2]"]])

  hist <- ggplot(diff_di, aes(x = diff)) +
    geom_histogram(color = "blue", fill = "blue")

  p <- ggarrange(p1, hist, ncol = 1, nrow = 2)
  print(p)

  print(quantile(diff, 0.025))
  print(quantile(diff, 0.975))
}

plot_posterior("model_ii_1", "HBB")
plot_posterior("model_ii_2", "HBB")
plot_posterior("model_ii_3", "HBB")

plot_posterior("model_ii_1", "LYZ")
plot_posterior("model_ii_2", "LYZ")
plot_posterior("model_ii_3", "LYZ")

## TEST
## model <- "model_ii_1"
## model <- "model_ii_2"
## model <- "model_ii_3"

## # model 2, 4 not that good. gene <- "HBB"
## gene <- "HBA2"
## gene <- "LYZ"
## model_fnm <- paste0("./", model, "_", gene, ".csv")
## fit <- read_stan_csv(model_fnm)
## sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
## summary(do.call(rbind, sampler_params), digits = 2)
## param <- "mu_g_ic"
## param <- "mu_g_di"
## param <- "mu"
## param <- c("mu_g_ic", "mu", "mu_g_di")
## ## par(mfrow=c(2,2))
## plot(fit, pars = c(param))
## plot(fit, pars = c(param), plotfun = "hist")
## traceplot(fit, pars = c(param), inc_warmup = TRUE, nrow = 2)
## ## ecdf(diff)(0.0)

# foldchange, influenced by small values quite a lot.
# fc <- exp(diff)
# hist(fc, breaks=100)
