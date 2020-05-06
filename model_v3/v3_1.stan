data {
    int<lower=10> N;  // number of cells
    int<lower=2> K; // number of individuals
    int<lower=10> scale; // scale factor
    int<lower=2> J; // number of condistions
    int<lower=1> G; // number of genes

    int x_cg[N, G]; //counts of gene(col) in cell(row)
    vector<lower=1>[N] x_; // total umk counts

    matrix[N, J] di; // conditions
    matrix[N, K] ic; // individual indicator
    matrix<lower=0>[K, K] nnegcor; // correlation between individuals.
}

transformed data {
    // hyperparam for variances of ind and cond.
    vector<lower=0>[K] alpha_ind = rep_vector(0.1, K);
    vector<lower=0>[K] beta_ind = rep_vector(0.1, K);
    real<lower=0> alpha_cond = 0.5;
    real<lower=0> beta_cond = 0.5;

    vector[G] gzeros = rep_vector(0.0, G);
    vector[K] izeros = rep_vector(0.0, K);
    real<lower=0> Lambda = 400;
    // add 1 for intercept.
    matrix[N,1] ones = rep_matrix(1, N, 1);
    matrix[N,1+J+K] x_dic = append_col(ones, append_col(di, ic));
}

parameters {
    real<lower=0> LambdaInd;
    vector<lower=0>[G] LambdaCond;
    vector[K] mu_g_ic;
    matrix[G, J] mu_g_di;
    real mu;
}

transformed parameters {
    matrix[1 + J+K, G] gbetas;
    for (g in 1:G) {
        gbetas[, g] = append_row(mu,append_row(mu_g_di[g]', mu_g_ic));
    }
}

model {

    mu ~ normal(0.0, Lambda);

    LambdaInd ~ inv_gamma(alpha_ind, beta_ind);
    mu_g_ic ~ multi_normal(izeros, LambdaInd * nnegcor);

    LambdaCond  ~ inv_gamma(alpha_cond, beta_cond);

    for (j in 1:J){
        mu_g_di[ ,j] ~ normal(gzeros, LambdaCond);
    }

    for (g in 1:G) {
        target += poisson_log_glm_lpmf(x_cg[, g] | x_dic, log(x_), gbetas[,g]);
    }
}
