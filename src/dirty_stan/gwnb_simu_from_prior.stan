// simulate data based on the prior.
// model: gene-wise negative binomial distribution.

data {
	int N; // num_of_cell
	int K; // num_of_ind
	int J; // num_of_cond
	vector<lower=100>[N] S; // cells' sum_of_cnt
	int<lower=1, upper=K> Ind[N];
  int<lower=1, upper=J> Cond[N]; // cond from 1

	// simulate model parameters from prior
	// real<lower=0> kappa2g;
	// real<lower=0> tau2g;
	// real<lower=0> phi2g;
	// real mug;
	// vector[K] muind;
	// vector[J] mucond;
	int y[N];


	// hyper parameters
	real muG0;
	real<lower=0> sigmaG0;
	real<lower=0> alphaKappa2G;
	real<lower=0> betaKappa2G;

	real<lower=0> alphaTau2G;
	real<lower=0> betaTau2G;

	real<lower=0> alphaPhi2G;
	real<lower=0> betaPhi2G;
}

transformed data {
	matrix[N, 1 + J + K] X = rep_matrix(0.0, N, 1+J+K);
	for (i in 1:N) {
		X[i, 1] = 1.0;
		X[i, 1 + Ind[i]] = 1.0;
		X[i, 1 + J + Cond[i]] = 1.0;
	}

	vector[N] logS = log(S);
}

parameters{
	real MuG;
	vector[K] MuIndRaw;
	vector[J] MuCondRaw;
	real<lower=0> Kappa2G;
	real<lower=0> Tau2G;
	real<lower=0> Phi2G;
}

transformed parameters {
	vector[K] MuInd;
	vector[J] MuCond;
	MuInd = sqrt(Kappa2G) * MuIndRaw;
	MuCond = sqrt(Tau2G) * MuCondRaw;
}

model {
	Kappa2G ~ inv_gamma(alphaKappa2G, betaKappa2G);
	Tau2G ~ inv_gamma(alphaTau2G, betaTau2G);
	Phi2G ~ inv_gamma(alphaPhi2G, betaPhi2G);

	MuG ~ normal(muG0, sigmaG0);
	MuIndRaw ~ std_normal(); //implicit MuInd[,k] ~ normal(0, diag(KappaG))
	MuCondRaw ~ std_normal(); //implicit MuCond[,j] ~ normal(0.0, diag(TauG))

	real nb_mu[N];
	for (i in 1:N) {
		nb_mu[i] = logS[i] + MuG + MuInd[Ind[i]] + MuCond[Cond[i]];
	}
	y ~ neg_binomial_2_log(nb_mu, Phi2G);
}

// generated quantities {
// 	int<lower=0, upper=1> lt_sim[4 + K + J];
// 	lt_sim[1] = Kappa2G < kappa2g;
// 	lt_sim[2] = Tau2G < tau2g;
// 	lt_sim[3] = Phi2G < phi2g;
// 	lt_sim[4] = MuG < mug;
// 	for (k in 1: K) {
// 		lt_sim[4 + k]  = (MuInd[k] < muind[k]);
// 	}
// 	for (j in 1:J) {
// 		lt_sim[4+K+j] = (MuCond[j] < mucond[j]);
// 	}
// }
