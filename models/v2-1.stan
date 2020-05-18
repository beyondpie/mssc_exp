// Model batch effect on gene module level

#include mydata.stan

parameters {
    real<lower=0> Sigma2G;
    real Mu;
    vector[G] MuGRaw;
    vector[P] MuIndMean;

    matrix[P, K] MuIndRaw;

    matrix[G, J] MuCondRaw;
    vector<lower=0>[P] Lambda2Ind;
}

transformed parameters {
    real SigmaG = sqrt(Sigma2G);
    vector[G] MuG = Mu + SigmaG * MuGRaw;
    vector[P] LambdaInd = sqrt(Lambda2Ind);
    matrix[P, K] MuInd;
    for (k in 1:K) {
        MuInd[,k] = dot_product(LambdaInd, MuIndRaw[, k]) + MuIndMean;
    }
    matrix[G, J] MuCond = sigmaMuCond * MuCondRaw;
}


model {
    Sigma2G ~ inv_gamma(alphaSigma2G, betaSigma2G);
    Lambda2Ind ~ inv_gamma(alphaLambda2Ind, betaLambda2Ind);

    Mu ~ normal(0.0, sigmaMu);
    MuIndMean ~ normal(0.0, sigmaMuIndMean);
    MuGRaw ~ std_normal(); // implicit MuG ~ normal(Mu, sqrt(Sigma2G))
    for (k in 1:K) {
        // implicit MuInd ~ normal(MuIndMean, diag(LambdaInd))
        MuIndRaw[, k] ~ std_normal();
    }
    for (j in 1:J) {
        MuCondRaw[, j] ~ std_normal(); // implicit MuCond ~ normal(0.0, sigmaMuCond)
    }
    matrix[G, K] MuIndonG = B * MuInd;
    vector[G] scores;
    vector[1 + J + K] W;
    for (g in 1:G) {
        W[1] = MuG[g];
        for (j in 1:J) {
            W[j+1] = MuCond[g, j];
        }
        for (k in 1: K) {
            W[k+J+1] = MuIndonG[g,k];
        }
        scores[g] = poisson_log_glm_lpmf(Xcg[,g] | X, logS, W);
    }
    target += sum(scores);
}
