// Model batch effect on gene module level.
// Gene module is directly estimated from PCA results of genes
//     on the cells from a cluster.
// Test the scalability on 1,000 genes together.


data {
    // number of cells
    int<lower=10> N;
    // number of individuals
    int<lower=2> K;
    // scRNAseq scale factor
    int<lower=1000> scale;
    // number of conditions
    int<lower=2> J;
    // number of genes;
    int<lower=1> G;
    // dimension of gene module matrix
    int<lower=1> P;

    // gene module matrix
    matrix[G, P] B;

    // count matrix of cell by gene
    int<lower=0> Xcg[N,G];

    // total counts per cell
    vector<lower=0>[N] S;

    // condition matrix of cell by cond
    matrix[N, J] XCond;
    // individual index matrix of cell by indi
    matrix[N, K] XInd;
}

transformed data{
    // fixed parameter for hLambdaG
    real<lower=0> alphahLambdaG = 5;
    // fixed parameter for hLambdaF
    real<lower=0> alphahLambdaF = 1;
    // fixed parameter for hLambdaCond
    real<lower=0> alphahLambdaCond = 1;

    // add ones for the total counts per cell;
    matrix[N, 1] Ones = rep_matrix(1, N, 1);
    // append glm data matrix
    matrix[N, 1 + J + K] X = append_col(Ones, append_col(XCond, XInd));

    // log total counts per cell;
    vector[N] logS = log(S);
}

parameters {
    // population-wise gene expressions: gene by 1
    // use raw for reparameterization
    vector[G] MuRaw;
    // variances per gene: gene by 1
    vector<lower=0>[G] LambdaG;
    // LambdaG share the same hyper prior
    real<lower=0.0001> hLambdaG;

    // individual effects on gene module: gene_module by inds;
    // use raw for reparameterization
    matrix[P, K] MuFRaw;
    // variances for individual effects: gene_module by inds;
    matrix<lower=0>[P, K] LambdaF;
    // shared variances for inds for a given gene module:
    // gene module by 1;
    vector<lower=0.0001>[P] hLambdaF;

    // gene expression level per condition: gene by cond
    // use raw for reparameterization
    matrix[G, J] MuCondRaw;
    // shared variances for genes under condition: gene by 1.
    vector<lower=0>[G] LambdaCond;
    // LambdaCond share the same hyper prior.
    real<lower=0.0001> hLambdaCond;
}

transformed parameters{
    // population-wise gene expression
    vector[G] Mu = LambdaG .* MuRaw;
    // individual effects on gene module
    matrix[P,K] MuF = LambdaF .* MuFRaw;
    // gene expression level per cond
    matrix[G,J] MuCond = diag_pre_multiply(LambdaCond, MuCondRaw);

}

model {
    // variances hyper priors
    hLambdaG ~ inv_gamma(alphahLambdaG, alphahLambdaG);
    hLambdaF ~ inv_gamma(alphahLambdaF, alphahLambdaF);
    hLambdaCond ~ inv_gamma(alphahLambdaCond, alphahLambdaCond);

    // variances prior
    LambdaG ~ inv_gamma(hLambdaG, hLambdaG);
    for (k in 1:K) {
        LambdaF[ ,k] ~ inv_gamma(hLambdaF, hLambdaF);
    }
    LambdaCond ~ inv_gamma(hLambdaCond, hLambdaCond);

    // due to reparameterization, std normal per element
    MuRaw ~ std_normal();
    for (k in 1:K) {
        MuFRaw[, k] ~ std_normal();
    }
    for (j in 1:J) {
        MuCondRaw[, j] ~ std_normal();
    }

    for (g in 1:G) {
        target += poisson_log_glm_lpmf(Xcg[, g] | X, logS,
                                       append_row(Mu[g],
                                                  append_row( MuCond[g]',
                                                              (B*MuF)[g]')));
    }
}
