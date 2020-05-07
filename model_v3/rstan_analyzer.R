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

plot_posterior <- function(modelnm, genenm,
                           param = c("mu_g_ic", "mu", "mu_g_di")) {
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
  p <- annotate_figure(p,
                       top = text_grob(paste(modelnm, genenm),
                                          color = "red",
                                          face = "bold", size = 14))
  print(p)

  print(quantile(diff, 0.025))
  print(quantile(diff, 0.975))
}
