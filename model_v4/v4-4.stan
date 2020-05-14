// Model batch effect on gene module level.
// Gene module is directly estimated from PCA results of genes
//     on the cells from a cluster.
// Remove hyper prior for variances, directly use predefined.
// Use indicators for Xind, Xcond, instead of one-hot matrix production.
// Test the scalability on 1,000 genes together.


data {
    // number of cells
    int N;
    // number of individuals
    int K;
    // scRNAseq scale factor
    int scale;
    // number of conditions
    int J;
    // number of genes;
    int G;
    // dimension of gene module matrix
    int P;
    // gene module matrix
    matrix[G, P] B;
    // int count 2d array of cell by gene
    int IXcg[N, G];
    // condition indicator: control 1, case 2
    int IXCond[N];
    // individual indicator: from 1 to 10
    int IXInd[N];
    // total counts per cell
    vector[N] S;
}

transformed data{
    // fixed parameter for SigmaG, this is std.
    real<lower=0> SigmaG = 20;
    // fixed parameter for hLambdaF
    real<lower=0> alphahLambdaF = 2;
    // fixed parameter for hLambdaCond
    real<lower=0> alphahLambdaCond = 2;

    // log total counts per cell;
    vector[N] LogS = log(S);
    matrix[G,N] LogSGN = rep_matrix(LogS, G)';

    vector[N*G] myzeros = rep_vector(0.0, N*G);
    vector[N*G] myones = rep_vector(1.0, N*G);
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
    matrix[G, N] MuIndGN = MuInd[:, IXInd];
    matrix[G, N] MuCondGN = MuCond[:, IXCond];
    matrix[G, N] Lgc = LogSGN + MuIndGN + MuCondGN;
    matrix[G*N, 1] X = rep_matrix(to_vector(Lgc),1);

    // NOTE: to_array_1d for int array is row major order.
    //       while to_vector for matrix is column major order.
    target += poisson_log_glm_lpmf(to_array_1d(IXcg) | X, 0.0, [1.0]');
}
