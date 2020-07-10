// Model batch effect gene-wise.
// This version model each gene indepdendently.

#include mydata.stan

parameters {
    vector[G] MuG;
    matrix[G, K] MuIndRaw;
    matrix[G, J] MuCondRaw;
    vector<lower=0>[G] Kappa2G;
    vector<lower=0>[G] Tau2G;
}

transformed parameters {
    vector[G] KappaG = sqrt(Kappa2G);
    matrix[G, K] MuInd;
    for (k in 1:K) {
        MuInd[, k] = KappaG .* MuIndRaw[, k];
    }
    vector[G] TauG = sqrt(Tau2G);
    matrix[G, J] MuCond;
    for (j in 1:J) {
        MuCond[ ,j] = TauG .* MuCondRaw[, j];
    }
}


model {
    Kappa2G ~ inv_gamma(alphaKappa2G, betaKappa2G);
    Tau2G ~ inv_gamma(alphaTauG, betaTauG);
    MuG ~ normal(0.0, sigmaG0);

    for (k in 1:K) {
        MuIndRaw[, k] ~ std_normal();//implicit MuInd[,k] ~ normal(0, diag(KappaG))
    }

    for (j in 1:J) {
        MuCondRaw[, j] ~ std_normal();//implicit MuCond[,j] ~ normal(0.0, diag(TauG))
    }

    vector[G] scores;
    vector[1 + J + K] W;
    for (g in 1:G) {
        W[1] = MuG[g];
        for (j in 1:J) {
            W[j+1] = MuCond[g, j];
        }
        for (k in 1: K) {
            W[k+J+1] = MuInd[g, k];
        }
        scores[g] = poisson_log_glm_lpmf(Xcg[,g] | X, logS, W);
    }
    target += sum(scores);
}
