// Model batch effect on gene module level

#include mydata.stan

parameters {
    vector[G] MuG;
    matrix[P, K] MuIndRaw;
    matrix[G, J] MuCondRaw;
    vector<lower=0>[P] Kappa2P;
    vector<lower=0>[G] Tau2G;
}

transformed parameters {
    vector[P] KappaP = sqrt(Kappa2P);
    matrix[P, K] MuInd;
    for (k in 1:K) {
        MuInd[, k] = KappaP .* MuIndRaw[, k];
    }
    vector[G] TauG = sqrt(Tau2G);
    matrix[G, J] MuCond;
    for (j in 1:J) {
        MuCond[, j]  = TauG .* MuCondRaw[, j];
    }
}


model {
    Kappa2P ~ inv_gamma(alphaKappaP, betaKappaP);
    Tau2G ~ inv_gamma(alphaTauG, betaTauG);
    MuG ~ normal(0.0, sigmaG0);
    for (k in 1:K) {
        MuIndRaw[, k] ~ std_normal();//implicit MuInd[,k] ~ normal(0, diag(KappaP))
    }
    for (j in 1:J) {
        MuCondRaw[, j] ~ std_normal();//implicit MuCond[,j] ~ normal(0.0, diag(TauG))
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
