data {
    int N; // number of cells
    int K; // number of individuals
    int J; // number of conditions
    int G; // number of genes
    int P; // number of gene module

    int Xcg[N, G]; // counts: cell by gene
    vector<lower=100>[N] S; // total count per cell

    matrix[N, J] XCond; // one-hot condition repr
    matrix[N, K] XInd; // one-hot individual repr

    matrix[G,P] B; // gene-moduel matrix

    int scale; // for scRNAseq

    real<lower=0> alphaSigma2G;
    real<lower=0> betaSigma2G;
    real<lower=0> sigmaMu;

    real<lower=0> alphaLambda2Ind;
    real<lower=0> betaLambda2Ind;

    real<lower=0> sigmaMuCond;

    int<lower=0, upper=1> GRAINSIZE;
}

transformed data {
    vector[N] logS = log(S);
    matrix[N, 1 + J + K] X = append_col(rep_matrix(1,N,1),
                                        append_col(XCond, XInd));
}
