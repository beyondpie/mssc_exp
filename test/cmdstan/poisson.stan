data {
    int<lower=0> N;
    int<lower=0> y[N];
}

transformed data {
    matrix[N, 1] x;
    for (i in 1:N) {
        x[i,1] = 1;
    }
}
parameters {
    real alpha;
    vector[1] beta;
    real<lower=0> l1;
    real<lower=0> l2;
    real<lower=0> l10;
    real<lower=0> l11;
}

transformed parameters {
    real ln_lambda;
}
model {
    l10 ~ gamma(1, 1);
    l11 ~ gamma(1,1);
    l1 ~ gamma(l10, l11);
    l2 ~ gamma(1, 1);
    alpha ~ normal(0,1);
    beta[1] ~ gamma(l1, l2);
    y ~ poisson_log_glm(x, alpha, beta);
    /* for (i in 1:N) { */
        /* target += poisson_log_glm_lpmf(y[i] | x[i], alpha, beta); */
    /* } */
}
