library(R6)


Hi <- R6Class(
  "Hi",
  public = list(
    a = 1,
    initialize = function(a = NULL){
      if (!is.null(a)) {
        self$a = a
      }
      self$hi()
    },
    hi = function() {
      print("Just say hi")
    }
  ))

a  <- Hi$new()
b <- Hi$new(a = 1000)


## test nb fitting
library(cmdstanr)
library(loo)
library(posterior)
snb <- cmdstanr::cmdstan_model(
  stan_file = here::here("src", "modelcheck", "high2", "stan",
                         "snb.stan")
)

n <- 100
mu <- 0.1
size <- 20
s <- rnorm(n = n, mean = 1.0, sd = 0.1)

y <- rnbinom(n = n, mu = s*exp(mu), size = 10)
hpg <- c(0.05, 0.05)

mu1 <- log(mean(y/s))
mu2 <- log(mean(y))- log(median(s))
r1 <- mu1^2 / (var(y) - mu1)
r2 <- mu2^2 / (var(y) - mu2)

opt1 <- snb$optimize(
  data = list(n = n, s = s, y = y, hpg = hpg),
  seed = 1L,
  iter = 1000,
  init = list(list(mu = mu1, r = r1)),
  algorithm = "lbfgs"
)
mu1
r1
opt1$mle()

vi1 <- snb$variational(
  data = list(n = n, s = s, y = y, hpg = hpg),
  seed = 1L,
  init = list(list(mu = mu1, r = r1)),
  ## init = NULL,
  output_samples = 2000,
  tol_rel_obj = 0.001
)
a <- vi1$draws("mu")
b <- vi1$draws("r")
c <- vi1$draws("r_u")
par(mfrow = c(3,1))
hist(a)
hist(b)
hist(c)
mean(a)
mean(b)
mean(c)

log_ratios <- vi1$lp() - vi1$lp_approx()
mypsis <- loo::psis(log_ratios = log_ratios,
                    r_eff = NA)
mypsis$diagnostic
lw <- weights(mypsis, log = FALSE, normalize = TRUE)

aa <- posterior::resample_draws(x = a, weights = lw,
                                method = "stratified")
mean(a)
mean(aa)

bb <- posterior::resample_draws(x=b, weights = lw,
                                method = "stratified")
cc <- posterior::resample_draws(x = c, weights = lw,
                                method = "stratified")
mean(b)
mean(bb)

mean(c)
mean(cc)


hpg <- c(1, 1)
opt2 <- snb$optimize(
  data = list(n = n, s = s, y = y, hpg = hpg),
  seed = 1L,
  iter = 1000,
  init = list(list(mu = mu2, r = r2)),
  algorithm = "lbfgs"
)
mu2
r2
opt2$mle()


opt3 <- snb$optimize(
  data = list(n = n, s = s, y = y, hpg = hpg),
  seed = 1L,
  iter = 1000,
  init = NULL,
  algorithm = "lbfgs"
)

opt3$mle()


## test model loading and run.
## * test
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

  vi_sampler <- run_hbnb_vi(data = data, ip = hi_params$ip)
  est_params <- lapply(nm_params, function(nm) {
    extract_vifit(vi_sampler, data, nm)
  })
}
