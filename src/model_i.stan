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
}

parameters {
    real alpha;
    real beta;

    real<lower=0> lambda_cg;

    real<lower=0> Lambda_cg;

    vector[I] mu_g_ic;
    vector[K] mu_g_di;

    vector[I] mu_0;

    /* cholesky_factor_cov[2] Lambda_0; */
    /* cov_matrix[I] Lambda_0; */
    vector<lower=0>[I] Lambda_0;

    vector[K] mu_g;

    /* cholesky_factor_cov[2] Lambda_g; */
    /* cov_matrix[2] Lambda_g; */
    vector<lower=0>[K] Lambda_g;

}

transformed parameters {
    real ln_xcg;
    ln_xcg = log(lambda_cg);
}

model {
    mu_g_di ~ normal(mu_g, Lambda_g);
    mu_g_ic ~ normal(mu_0, Lambda_0);
    ln_xcg ~ normal(ic * mu_g_ic + di * mu_g_di, Lambda_cg);
    lambda_cg ~ inv_gamma(alpha, beta);
    x_cg ~ poisson(x_ * lambda_cg);
}
