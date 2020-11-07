// Hierarchical Bayesian based negative binomial model

data {
	int n;
	int k;
	int j;
	int g;

	row_vector<lower=0>[n] s;
	int<lower=1, upper=j> cond[n];
	int<lower=1, upper=k> ind[n];
	int<lower=0> y[g, n];

	vector[g] mu0;
	vector<lower=0>[2] hp_varofmu;
	vector<lower=0>[2] hp_alpha_r;
	vector<lower=0>[2] hp_beta_r;

	vector<lower=0>[2] hp_alpha_varofind;
	vector<lower=0>[2] hp_beta_varofind;
	vector<lower=0>[2] hp_varofcond;
}

transformed data {
	matrix[g, n] logs = rep_matrix(log(s), g);
	// row-major order
	int y1d[g*n] = to_array_1d(y);
}

parameters {
	vector<lower=0>[2] hp_r;
	vector<lower=0>[g] nb_r;

	real<lower=0> varofmu;
	vector[g] raw_mu;

	real<lower=0> varofcond;
	matrix[g, j] raw_mu_cond;

	vector<lower=0>[2] hp_varofind;
	vector<lower=0>[g] varofind;

	matrix[g, k] raw_mu_ind;
}

// transformed parameters {
// 	vector[g] mu = raw_mu * sqrt(varofmu)+ mu0;
// 	matrix[g, j] mu_cond = raw_mu_cond * sqrt(varofcond);
// 	matrix[g, k] mu_ind = diag_matrix(sqrt(varofind)) * raw_mu_ind;
// }

model {
	hp_r[1] ~ gamma(hp_alpha_r[1], hp_alpha_r[2]);
	hp_r[2] ~ gamma(hp_beta_r[1], hp_beta_r[2]);
	nb_r ~ gamma(hp_r[1], hp_r[2]);

	varofmu ~ inv_gamma(hp_varofmu[1], hp_varofmu[2]);
	raw_mu ~ std_normal();

	hp_varofind[1] ~ gamma(hp_alpha_varofind[1], hp_alpha_varofind[2]);
	hp_varofind[2] ~ gamma(hp_beta_varofind[1], hp_beta_varofind[2]);
	varofind ~ inv_gamma(hp_varofind[1], hp_varofind[2]);

	varofcond ~ inv_gamma(hp_varofcond[1], hp_varofcond[2]);

	to_vector(raw_mu_cond) ~ std_normal();
	to_vector(raw_mu_ind) ~ std_normal();

	// need R program to calculate mu, mu_cond, and mu_ind.
	vector[g] mu = raw_mu * sqrt(varofmu)+ mu0;
	matrix[g, j] mu_cond = raw_mu_cond * sqrt(varofcond);
	matrix[g, k] mu_ind = diag_matrix(sqrt(varofind)) * raw_mu_ind;

	matrix[g, n] lambda = (logs + rep_matrix(mu, n) + mu_cond[ , cond] + mu_ind[ , ind]);
	matrix[g, n] nb_rr = rep_matrix(nb_r, n);
	// to_vector is column-major order.
	y1d ~ neg_binomial_2_log(to_vector(lambda), to_vector(nb_rr));
}
