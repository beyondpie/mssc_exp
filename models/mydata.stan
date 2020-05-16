data {
    int N; // number of cells
    int K; // number of individuals
    int J; // number of conditions
    int G; // number of genes

    int Xcg[N, G]; // counts: cell by gene
    vector<lower=100>[N] S; // total count per cell

    matrix[N, J] XCond; // one-hot condition repr
    matrix[N, K] XInd; // one-hot individual repr

    int scale; // for scRNAseq

    real<lower=0> alphaSigma2G;
    real<lower=0> betaSigma2G;

    real<lower=0> sigmaMu;
}

