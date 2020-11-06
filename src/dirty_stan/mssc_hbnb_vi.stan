// Hierarchial Bayesian based negative binomial model.

data {
	int n; // num of cell
	int k; // num of individual
	int j; // num of cond
	int g; // num of gene

	vector<lower=0>[n] s; // sum of count in cells
	int<lower=1, upper=j> cond[n]; // experiment index
	int<lower=1, upper=k> ind[n]; // individual index
	int<lower=0> y[g, n];
	// int<lower=0> y[g, n]; // observations

	// hyper parameters
	vector<upper=0>[g] mu0; // mean of scaled log mean of exp;
	vector<lower=0>[2] hp_varofmu;

	vector<lower=0>[2] hp_r;

	vector<lower=0>[2] hp_alpha_varofind;
	vector<lower=0>[2] hp_beta_varofind;

	vector<lower=0>[2] hp_varofcond;
}

transformed data {
	vector[n] logs = log(s);
}

parameters{
	vector[g] nb_r;

	real<lower=0> varofmu;
	vector[g] mu;

	real<lower=0> varofcond;
	vector[g] mu_cond[j];

	vector<lower=0>[g] hp_varofind[2];
	vector<lower=0>[g] varofind;
	vector[g] mu_ind[k];
}


model{
	// gamma or log-normal (DESeq2 use log-normal)
	nb_r ~ gamma(hp_r[1], hp_r[2]);
	varofmu ~ inv_gamma(hp_varofmu[1], hp_varofmu[2]);
	mu ~ normal(mu0, sqrt(varofmu));

	hp_varofind[1] ~ gamma(hp_alpha_varofind[1], hp_alpha_varofind[2]);
	hp_varofind[2] ~ gamma(hp_beta_varofind[1], hp_beta_varofind[2]);
	varofind ~ inv_gamma(hp_varofind[1], hp_varofind[2]);

	varofcond ~ inv_gamma(hp_varofcond[1], hp_varofcond[2]);

	for (i in 1:j) {
		mu_cond[i] ~ normal(0.0, sqrt(varofcond));
	}
	for (i in 1:k) {
		mu_ind[[i] ~ normal(0.0, sqrt(varofind));
	}
	vector[g] sumofmu[n] = rep_array(mu, n) + mu_cond[cond] + mu_ind[ind];
	for (i in 1:g) {
		y[i] ~ neg_binomial_2_log(logs + sumofmu[i], nb_r[i]);
	}
}
