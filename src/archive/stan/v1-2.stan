// Model batch effect gene-wise
// Add hyper prior for the kappa and tau

#include mydata.stan

parameters {
	vector[G] MuG;
	matrix[G, K] MuIndRaw;
	matrix[G, J] MuCondRaw;
	vector<lower=0>[G] Kappa2G;
	vector<lower=0>[G] Tau2G;
	real<lower=0> AlphaKappa;
	real<lower=0> BetaKappa;
	real<lower=0> AlphaTau;
	real<lower=0> BetaTau;
}

transformed parameters {
	vector[G] KappaG = sqrt(Kappa2G);
	matrix[G, K] MuInd;
	vector[G] TauG = sqrt(Tau2G);
	matrix[G, J] MuCond;
	for (k in 1:K) {
		MuInd[,k] = KappaG .* MuIndRaw[, k];
	}
	for (j in 1:J) {
		MuCond[ ,j] = TauG .* MuCondRaw[, j];
	}
}


model {
	vector[G] scores;
	vector[1 + J + K] W;
	AlphaKappa ~ gamma(alphaAlphaKappa, betaAlphaKappa);
	BetaKappa ~ gamma(alphaBetaKappa, betaBetaKappa);
	Kappa2G ~ gamma(AlphaKappa, BetaKappa);

	AlphaTau ~ gamma(alphaAlphaTau, betaAlphaTau);
	BetaTau ~ gamma(alphaBetaTau, betaBetaTau);
	Tau2G ~ inv_gamma(AlphaTau, BetaTau);

	MuG ~ normal(0.0, sigmaG0);
	for (k in 1:K) {
		MuIndRaw[, k] ~ std_normal();//implicit MuInd ~ normal(0, diag(KappaG))
	}
	for (j in 1:J) {
		MuCondRaw[, j] ~ std_normal();//implicit MuCond ~ normal(0, sigmaMuCond)
	}
	for (g in 1:G) {
		W[1] = MuG[g];
		for (j in 1:J) {
			W[j+1] = MuCond[g, j];
		}
		for (k in 1: K) {
			W[k+J+1] = MuInd[g,k];
		}
		scores[g] = poisson_log_glm_lpmf(Xcg[,g] | X, logS, W);
	}
	target += sum(scores);
}
