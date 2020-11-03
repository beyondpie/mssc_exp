// Hierarchial Bayesian based negative binomial model.
// integrate out mu_ind

data {
	int n; // num of cell
	int k; // num of individual
	int j; // num of cond
	int g; // num of gene

	vector<lower=0>[n] s; // sum of count in cells
	int<lower=1, upper=j> cond[n]; // experiment index
	int<lower=0> y[g, n];
	// int<lower=0> y[g, n]; // observations

	// hyper parameters
	vector<upper=0>[g] mu0; // mean of scaled log mean of exp;
	vector<lower=0>[2] hp_varofmu;

	vector<lower=0>[2] hp_r;

	vector<lower=0>[2] hp_alpha_varofind;
	vector<lower=0>[2] hp_beta_varofind;

	vector<lower=0>[2] hp_alpha_varofcond;
	vector<lower=0>[2] hp_beta_varofcond;
}

transformed data {
	vector[n] logs = log(s);
}

parameters{
	vector[g] nb_r;

	real<lower=0> varofmu;
	vector[g] raw_mu;

	vector<lower=0>[g] hp_varofcond[2];
	vector<lower=0>[g] varofcond;
	vector[g] raw_mu_cond[j];

	vector<lower=0>[g] hp_varofind[2];
	vector<lower=0>[g] varofind;

	vector<lower=0>[n] raw_sumofmu[g];
}

transformed parameters {
	vector[g] mu = raw_mu * sqrt(varofmu) + mu0;
	vector[g] mu_cond[j];
	for (i in 1:j) {
		mu_cond[i] = raw_mu_cond[i] .* sqrt(varofcond);
	}
	vector<lower=0>[n] sumofmu[g];
	for (i in 1:g) {
			sumofmu[i] = mu[i] + mu_cond[i][cond] + raw_sumofmu[g] * sqrt(varofind[i]);
	}
}

model{
	// gamma or log-normal (DESeq2 use log-normal)
	nb_r ~ gamma(hp_r[1], hp_r[2]);
	varofmu ~ inv_gamma(hp_varofmu[1], hp_varofmu[2]);
	raw_mu ~ std_normal(); // implicit mu ~ normal(mu0, sqrt(varofmu));

	hp_varofind[1] ~ gamma(hp_alpha_varofind[1], hp_alpha_varofind[2]);
	hp_varofind[2] ~ gamma(hp_beta_varofind[1], hp_beta_varofind[2]);
	varofind ~ inv_gamma(hp_varofind[1], hp_varofind[2]);

	hp_varofcond[1] ~ gamma(hp_alpha_varofcond[1], hp_alpha_varofcond[2]);
	hp_varofcond[2] ~ gamma(hp_beta_varofcond[1], hp_beta_varofcond[2]);
	varofcond ~ inv_gamma(hp_varofcond[1], hp_varofcond[2]);

	for (i in 1:j) {
		raw_mu_cond[j] ~ std_normal(); // implicit mu_cond[j] ~ normal(0.0, sqrt(varofcond))
	}

	for (i in 1:g) {
		// implicit sumofmu[g] ~ normal(mu + mucond[cond], sqrt(varofind[i]))
		raw_sumofmu[g] ~ std_normal();
	}

	for (i in 1:g) {
		y[i] ~ neg_binomial_2_log(logs + sumofmu[i], nb_r[i]);
	}
}
