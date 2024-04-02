data {
    int N;
    vector[N] y;
    real<lower=0> sigmu;
}

parameters {
    real muraw;
    real<lower=0> sigma;
}

transformed parameters {
    /* real mu = sigmu * muraw; */
    real mu = sigma * muraw;
}

model {
    muraw ~ std_normal();
    sigma ~ inv_gamma(1,1);
    print("sigma is: ", sigma);
    y ~ normal(mu, sigma);
}
