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
}

transformed parameters {
    real ln_lambda;
}
model {
    y ~ poisson_log_glm(x, alpha, beta);
    /* for (i in 1:N) { */
        /* target += poisson_log_glm_lpmf(y[i] | x[i], alpha, beta); */
    /* } */
}
