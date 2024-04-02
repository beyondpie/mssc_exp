// Model batch effect on gene module level
// Add hyper prior for the kappa and tau

#include mydata.stan

parameters {
	vector[G] MuG;
	matrix[P, K] MuIndRaw;
	matrix[G, J] MuCondRaw;
	vector<lower=0>[P] Kappa2P;
	vector<lower=0>[G] Tau2G;
	real<lower=0> AlphaKappa;
	real<lower=0> BetaKappa;
	real<lower=0> AlphaTau;
	real<lower=0> BetaTau;
}

transformed parameters {
	vector[P] KappaP = sqrt(Kappa2P);
	matrix[P, K] MuInd;
	vector[G] TauG = sqrt(Tau2G);
	matrix[G, J] MuCond;
	for (k in 1:K) {
		MuInd[, k] = KappaP .* MuIndRaw[, k];
	}
	for (j in 1:J) {
		MuCond[, j]  = TauG .* MuCondRaw[, j];
	}
}


model {
	matrix[G, K] MuIndonG = B * MuInd;
	vector[G] scores;
	vector[1 + J + K] W;

	AlphaKappa ~ gamma(alphaAlphaKappa, betaAlphaKappa);
	BetaKappa ~ gamma(alphaBetaKappa, betaBetaKappa);
	Kappa2P ~ inv_gamma(AlphaKappa, BetaKappa);

	AlphaTau ~ gamma(alphaAlphaTau, betaAlphaTau);
	BetaTau ~ gamma(alphaBetaTau, betaBetaTau);
	Tau2G ~ inv_gamma(AlphaTau, BetaTau);

	MuG ~ normal(0.0, sigmaG0);
	for (k in 1:K) {
		MuIndRaw[, k] ~ std_normal();//implicit MuInd[,k] ~ normal(0, diag(KappaP))
	}
	for (j in 1:J) {
		MuCondRaw[, j] ~ std_normal();//implicit MuCond[,j] ~ normal(0.0, diag(TauG))
	}
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
