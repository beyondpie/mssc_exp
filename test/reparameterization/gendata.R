library(rstan)

## * generate train datas

mu <- 0.0
sigmu <- 10

sigma <- 1

N <- 100
y <- rnorm(N, mean=mu, sd = sigma)

rstan::stan_rdump(c("N", "y", "sigmu"), file="data.rdump")

## * analyze the results
