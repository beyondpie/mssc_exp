data {
    int<lower=0> N;  // total sample size
    /* int<lower=0> K;   // number of covariats. */
    int<lower=0> I; // number of individuals
    int<lower=10> scale; // scale factor
    matrix[N, 2] ds; // conditions
    matrix[N, I] ic; // individual indicator
    vector<N> x_; // total UMI counts
    vector<N> x_cg; // in
}


parameters {
    real alpha;
    real beta;

    real<lower=0> lambda_cg;
    vector[I] mu_g_ic;
    vector[2] mu_g_di;

    vector[I] mu_0;
    vector[I] Lambda_0;
    vector[2] mu_g;
    vector[2] Lambda_g;

}

transformed parameters {
    vector[N] ln_xcg;
    ln_xcg <- ln(lambda_cg);
}



model {
    lambda_cg ~ invgamma_lpdf(alpha, beta)
    mu_g_di ~ normal_lpdf(mu_g, Lambda_g)
    mu_g_ic ~ normal_lpdf(mu_0, Lambda_0)
    ln_xcg ~ normal_lpdf(ic * mu_g_ic + ds * mu_g_di, lambda_cg)
    x_cg ~ poisson_lpdf(x_ * lambda_cg)
}
