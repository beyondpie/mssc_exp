// Model batch effect on gene module level.
// Gene module is directly estimated from PCA results of genes
//     on the cells from a cluster.
// Remove hyper prior for variances, directly use predefined.
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
    real<lower=0> SigmaG = 20;
    // fixed parameter for hLambdaF
    real<lower=0> alphahLambdaF = 2;
    // fixed parameter for hLambdaCond
    real<lower=0> alphahLambdaCond = 2;

    // add ones for the total counts per cell;
    matrix[N, 1] Ones = rep_matrix(1, N, 1);
    // append glm data matrix
    matrix[N, 1 + J + K] X = append_col(Ones, append_col(XCond, XInd));

    // log total counts per cell;
    vector[N] logS = log(S);
}

parameters {
    // population-wise gene expressions: gene by 1
    // variances per gene: gene by 1
    vector[G] Mu;

    // individual effects on gene module: gene_module by inds;
    // use raw for reparameterization
    matrix[P, K] MuFRaw;
    // variances for individual effects: gene_module by 1;
    vector<lower=0>[P] LambdaF;

    // gene expression level per condition: gene by cond
    // use raw for reparameterization
    matrix[G, J] MuCondRaw;
    // shared variances for genes under condition: gene by 1.
    vector<lower=0>[G] LambdaCond;
}

transformed parameters{
    // individual effects on gene module
    matrix[P,K] MuF = diag_pre_multiply(sqrt(LambdaF), MuFRaw);
    // gene expression level per cond
    matrix[G,J] MuCond = diag_pre_multiply(sqrt(LambdaCond), MuCondRaw);

}

model {

    // variances prior
    LambdaF ~ inv_gamma(alphahLambdaF, alphahLambdaF);
    LambdaCond ~ inv_gamma(alphahLambdaCond, alphahLambdaCond);

    Mu ~ normal(0.0, SigmaG);
    for (k in 1:K) {
        MuFRaw[, k] ~ std_normal();
    }
    for (j in 1:J) {
        MuCondRaw[, j] ~ std_normal();
    }
    matrix[G, K] MuInd = B * MuF;
    for (g in 1:G) {
        vector[1 + J + K] weight = append_row(Mu[g], append_row(MuCond[g]', MuInd[g]'));
        /* print("The weight is ", weight); */
        target += poisson_log_glm_lpmf(Xcg[, g] | X, logS,weight);
    }
}
