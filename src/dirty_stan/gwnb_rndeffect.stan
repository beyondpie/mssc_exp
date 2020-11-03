// gene-wise negative binomial distribution.
// integrate out mu_ind

data {
	int N; // num_of_cell
	int K; // num_of_ind
	int J; // num_of_cond
	vector<lower=100>[N] S; // cells' sum_of_cnt
	int<lower=1, upper=J> Cond[N]; // cond from 1
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
	vector[N] logS = log(S);
}

parameters {
	real MuG;
	vector[J] MuCondRaw;
	vector[N] MuRaw;
	real<lower=0> Kappa2G;
	real<lower=0> Tau2G;
	real<lower=0> Phi2G;
}

transformed parameters {
	vector[J] MuCond;
	MuCond = sqrt(Tau2G) * MuCondRaw;
	vector[N] Mu;
	for (i in 1:N)  {
		// this might be wrong.
		Mu[i] = MuG + MuCond[Cond[i]] + sqrt(Kappa2G) * MuRaw[i];
	}
}

model {
	Kappa2G ~ inv_gamma(alphaKappa2G, betaKappa2G);
	Tau2G ~ inv_gamma(alphaTau2G, betaTau2G);
	// Phi2G ~ inv_gamma(alphaPhi2G, betaPhi2G);
	Phi2G ~ gamma(alphaPhi2G, betaPhi2G);

	MuG ~ normal(muG0, sigmaG0);
	MuCondRaw ~ std_normal(); // implicit MuCond ~ normal(0, diag(TauG))
	MuRaw ~ std_normal(); // implicit Mu ~ normal(MuG + MuCond, diag(Kappa2G))

	vector[N] nb_mu = logS + Mu;
	y ~ neg_binomial_2_log(nb_mu, Phi2G);
}



