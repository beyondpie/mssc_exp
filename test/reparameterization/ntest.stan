data {
    int N;
    vector[N] y;
    real<lower=0> sigmu;
}

parameters {
    real mu;
    real<lower=0> sigma;
}

model {
    /* mu ~ normal(0.0, sigmu); */
    mu ~ normal(0.0, sigma);
    sigma ~ inv_gamma(1,1);
    print("sigma is: ", sigma);
    y ~ normal(mu, sigma);
}
