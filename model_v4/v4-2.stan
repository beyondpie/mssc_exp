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

    // count matrix of gene by cell
    int<lower=0> Xgc[G, N];

    // total counts per cell
    vector<lower=0>[N] S;

    // condition matrix of cell by cond
    matrix[N, J] XCond;
    // individual index matrix of cell by indi
    matrix[N, K] XInd;
}

transformed data{
    // fixed parameter for SigmaG, this is std.
    real<lower=0> SigmaG = 20;
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
    matrix[N, G] logSs = rep_matrix(logS, G);
}

parameters {
    // gene-wise population mean expressions: gene by 1
    // mean per gene: gene by 1
    // use raw for reparameterization
    vector[G] MuRaw;

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
    // gene-wise population mean expression: gene by 1
    vector[G] Mu = SigmaG * MuRaw;
    // individual effects on gene module
    matrix[P,K] MuF = diag_pre_multiply(sqrt(LambdaF), MuFRaw);
    // gene expression level per cond
    matrix[G,J] MuCond = diag_pre_multiply(sqrt(LambdaCond), MuCondRaw);

}

model {
    // variances prior
    LambdaF ~ inv_gamma(alphahLambdaF, alphahLambdaF);
    LambdaCond ~ inv_gamma(alphahLambdaCond, alphahLambdaCond);

    MuRaw ~ std_normal();
    to_vector(MuFRaw) ~ std_normal();
    to_vector(MuCondRaw) ~ std_normal();

    matrix[G, K] MuInd = B * MuF;
    matrix[G, 1 + J + K] W = append_col(Mu, append_col(MuCond, MuInd));
    matrix[N,G] Lcg = X * W' + logSs;

    // NOTE: to_array_1d for int array is row major order.
    //       while to_vector for matrix is column major order.
    to_array_1d(Xgc) ~ poisson_log(to_vector(Lcg));

    // TODO: test
    // This version seems to be faster
    /* for (g in 1:G) { */
        /* vector[1 + J + K] weight = append_row(Mu[g], */
                                              /* append_row(MuCond[g]', MuInd[g]')); */
        /* print("The weight is ", weight); */
        /* target += poisson_log_glm_lpmf(Xcg[, g] | X, logS,weight); */
    /* } */
}
