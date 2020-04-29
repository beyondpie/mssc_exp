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
    real<lower=0> alpha;
    real<lower=0> beta;

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
    for (i in 1:K) {
        mu_g_di[i] ~ normal(mu_g[i], Lambda_g[i]);
    }

    for (i in 1:I) {
        mu_g_ic[i] ~ normal(mu_0[i], Lambda_0[i]);
    }

    Lambda_cg ~ inv_gamma(alpha, beta);

    lambda_cg ~ lognormal(ic * mu_g_ic + di * mu_g_di, Lambda_cg);
    x_cg ~ poisson(x_ * lambda_cg);
}
