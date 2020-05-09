// exactly model II
// global mean shared among genes
// individual effect shared among genes.
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
    vector<lower=0>[G] alpha_cond = rep_vector(0.5, G);
    vector<lower=0>[G] beta_cond = rep_vector(0.5, G);
    real<lower=0> alpha_ind = 0.1;
    real<lower=0> beta_ind = 0.1;

    vector[K] izeros = rep_vector(0.0, K);
    real Lambda = 20;
    // add 1 for intercept.
    matrix[N,1] ones = rep_matrix(1, N, 1);
    matrix[N,1+J+K] x_dic = append_col(ones, append_col(di, ic));

    // reparameterization for MVN
    matrix[K, K] iL = cholesky_decompose(nnegcor);
}

parameters {
    real<lower=0> LambdaInd;
    vector<lower=0>[G] LambdaCond;
    vector[K] mu_ic_raw;
    matrix[G, J] mu_di_raw;
    real mu;
}

transformed parameters {
    matrix[G, J] mu_g_di = diag_pre_multiply(LambdaCond, mu_di_raw);
    vector[K] mu_g_ic = LambdaInd * ( iL * mu_ic_raw);
    matrix[1 + J+K, G] gbetas;
    for (g in 1:G) {
        gbetas[, g] = append_row(mu,append_row(mu_g_di[g]', mu_g_ic));
    }
}

model {

    mu ~ normal(0.0, Lambda);

    LambdaInd ~ inv_gamma(alpha_ind, beta_ind);
    LambdaCond  ~ inv_gamma(alpha_cond, beta_cond);

    // opotimize this using reparameterization.
    mu_ic_raw ~ std_normal();
    for (j in 1:J){
        mu_di_raw[,j] ~ std_normal();
    }

    // change to map-reduce for parallel.
    for (g in 1:G) {
        target += poisson_log_glm_lpmf(x_cg[, g] | x_dic, log(x_), gbetas[,g]);
    }
}

// add block for diff of mu_g_di under different conditions.
