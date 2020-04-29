data {
    int<lower=0> N;
    int<lower=0> y[N];
}
parameters {
    real<lower=0> lambda;
}

transformed parameters {
    real ln_lambda;
    ln_lambda = log(lambda);
}
model {
    ln_lambda ~ normal(0,1); // uniform prior on interval 0,1
    y ~ poisson(lambda);
}
