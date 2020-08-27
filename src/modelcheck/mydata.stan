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
}

transformed data {
    int scale = 10000; // for scRNAseq

    vector<lower=0>[G] sigmaG0 = rep_vector(10.0, G);

    real<lower=0> alphaKappa2G = 1.0;
    real<lower=0> betaKappa2G = 1.0;

    real<lower=0> alphaTauG = 1.0;
    real<lower=0> betaTauG = 1.0;

    real<lower=0> alphaAlphaKappa = 1.0;
    real<lower=0> betaAlphaKappa = 1.0;
    real<lower=0> alphaBetaKappa = 1.0;
    real<lower=0> betaBetaKappa = 1.0;

    real<lower=0> alphaAlphaTau = 1.0;
    real<lower=0> betaAlphaTau = 1.0;
    real<lower=0> alphaBetaTau = 1.0;
    real<lower=0> betaBetaTau = 1.0;

    real<lower=0> alphaKappaP = 1.0;
    real<lower=0> betaKappaP = 1.0;

  	real<lower=0> alphaPhi2G = 1.0;
	  real<lower=0> betaPhi2G = 1.0;

	  real<lower=0> alphaSigma2G = 1.0;
	  real<lower=0> betaSigma2G = 1.0;

    int<lower=0, upper=1> GRAINSIZE = 1;

    vector[N] logS = log(S);
    matrix[N, 1 + J + K] X = append_col(rep_matrix(1,N,1),
                                        append_col(XCond, XInd));
}
