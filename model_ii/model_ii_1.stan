data {
    int<lower=0> N;  // total sample size
    /* int<lower=0> J;   // number of covariats. */
    int<lower=0> K; // number of individuals

    /* int<lower=0> T; // number of cell type */

    int<lower=10> scale; // scale factor

    int J; // number of condistions

    matrix[N, J] di; // conditions
    matrix[N, K] ic; // individual indicator
    vector[N] x_; // total UMK counts
    int x_cg[N]; // read counts for gene g in different cells.
}

// define hyper parameters here.
transformed data {
    real<lower=0> alpha_0;
    real<lower=0> beta_0;
    alpha_0 = 1;
    beta_0 = 1;

    matrix[N,1] ones;
    ones = rep_matrix(1, N, 1);

    matrix[N, 1+J+K] x_dic;
    x_dic = append_col(ones, append_col(di, ic));
}

parameters {
    real mu_g_ic[K];
    real mu_g_di[J];
    /* real mu_0; */
    real<lower=0> LambdaK;
    real<lower=0> LabmdaJ;

    real mu[1];
}

transformed parameters {
    vector[J+K] betas;
    betas = to_vector(append_array(mu, append_array(mu_g_di, mu_g_ic)));
}

model {
    /* choose variance as var estimated from samples, cover the estimated mean,
       and multiplied by constant */
    mu[1] ~ normal(0.0, 400);
    LabmdaJ ~ inv_gamma(0.001, 0.001);
    for (i in 1:J) {
        mu_g_di[i] ~ normal(0.0, LabmdaJ);
    }
    LambdaK ~ inv_gamma(0.001, 0.001);
    for (i in 1:K) {
        mu_g_ic[i] ~ normal(0.0, LambdaK);
    }
    target += poisson_log_glm_lpmf(x_cg | x_dic, log(x_), betas);
}
