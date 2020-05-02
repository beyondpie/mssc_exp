data {
    int<lower=0> N;  // total sample size
    /* int<lower=0> K;   // number of covariats. */
    int<lower=0> I; // number of individuals

    /* int<lower=0> T; // number of cell type */

    int<lower=10> scale; // scale factor

    int K; // number of condistions

    matrix[N, K] di; // conditions
    matrix[N, I] ic; // individual indicator
    vector[N] x_; // total UMI counts
    int x_cg[N]; // read counts for gene g in different cells.
}

// define hyper parameters here.
transformed data {
    real<lower=0> alpha_0;
    real<lower=0> beta_0;
    alpha_0 = 1;
    beta_0 = 1;

    matrix[N, K+I] x_dic;
    x_dic = append_col(di, ic);
}

parameters {
    real mu_g_ic[I];
    real mu_g_di[K];
    /* real mu_0; */
    real mu_g;
    real<lower=0> Lambda_0;
    real<lower=0> Lambda_g;

}

transformed parameters {
    vector[K+I] betas;
    betas = to_vector(append_array(mu_g_di, mu_g_ic));
}

model {
    Lambda_g ~ inv_gamma(0.001, 0.001);
    mu_g ~ normal(0.0, 400);
    for (i in 1:K) {
        mu_g_di[i] ~ normal(mu_g, Lambda_g);
    }
    Lambda_0 ~ inv_gamma(0.001, 0.001);
    for (i in 1:I) {
        mu_g_ic[i] ~ normal(0.0, Lambda_0);
    }
    target += poisson_log_glm_lpmf(x_cg | x_dic, log(x_), betas);
}
