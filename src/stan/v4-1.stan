// Model batch effect gene-wise
// This version model each gene independently.
// Use Negative-binomial distirbution.
// log scaled mean models as log-normal distribution.

#include mydata.stan

parameters {
	vector[G] MuG;
	matrix[G, K] MuInd;
	matrix[G, J] MuCond;

	matrix<lower=0>[N, G] ExpMu;

	vector<lower=0>[G] Phi2G;
	vector<lower=0>[G] Sigma2G;
}

transformed parameters {
	vector[G] SigmaG = sqrt(Sigma2G);
}

model {
	vector[G] scores;
	vector[1 + J + K] W;

	Sigma2G ~ inv_gamma(alphaSigma2G, betaSigma2G);
	Phi2G ~ inv_gamma(alphaPhi2G, betaPhi2G);
	// Beta (MuG, MuIndRaw, MuCondRaw) follows
	// an non-informative improper prior: Uniform
	for(g in 1:G) {
		W[1] = MuG[g];
		for (j in 1:J) {
			W[j+1] = MuCond[g, j];
		}
		for (k in 1:K) {
			W[k + J+1] = MuInd[g, k];
		}
		ExpMu[, g] ~ lognormal(X * W, SigmaG[g]);
		scores[g] =  neg_binomial_2_lpmf(Xcg[ ,g] | dot_product(logS, ExpMu[, g]), Phi2G[g]);
	}
	target += sum(scores);
}
