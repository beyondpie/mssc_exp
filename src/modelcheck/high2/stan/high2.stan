// Hierarchical Bayesian based negative binomial model

// - model the individual effects in gene-wise way.
// - model mu in hierarchical structure
//   - where the mean of mu follows a non-informative prior
// - model mucond:
//   - two condition, but only one delta is estimated
// - no more inv-gamma
// - all the gamma of gamma-related prior are replaced with uniform
//   distribution

data {
	int ncell;
	int nind;
	int ngene;

	row_vector<lower=0>[ncell]  s;
	int<lower=1, upper=2> cond[ncell];
	int<lower=1, upper=nind> ind[ncell];

	// shape: cell by gene for an two-dim array
	int<lower=0> y[ncell, ngene];

	vector<lower=0>[2] hp_varofmu;

	vector<lower=0>[2] hp_alpha_r;
	vector<lower=0>[2] hp_beta_r;
	vector<lower=0>[2] hp_uniform_r;

	vector<lower=0>[2] hp_alpha_varofind;
	vector<lower=0>[2] hp_beta_varofind;
	vector<lower=0>[2] hp_uniform_varofind;

	vector<lower=0>[2] hp_varofcond;

}

transformed data {
	// rep_matrix of row_vector: in row-wise
	matrix[ngene, ncell] logs = rep_matrix(log(s), ngene);
	// row-major order:
	// [g1_c1, g2_c1, ..., g1_c2, g2_c2, ..., g1_cn, g2_cn]
	int y1d[ngene * ncell] = to_array_1d(y);
}

parameters {
	// vector<lower=0> hp_r;
	vector<lower=hp_uniform_r[1], upper=hp_uniform_r[2]>[2] hp_r;
	vector<lower=0>[ngene] nb_r;

	real<lower=0> varofmu;
	real centerofmu;
	vector[ngene] raw_mu;

	real<lower=0> varofcond;
	// gaussian distributuon
	vector[ngene] raw_mu_cond;

	// vector<lower=0> hp_varofind;
	vector<lower=hp_uniform_varofind[1], upper=hp_uniform_varofind[2]>[2] hp_varofind;
	vector<lower=0>[nind] varofind;
	matrix[g,k] raw_mu_ind;
}

transformed parameters {
	// cmdstanr accelarates the draws
	vector[g] mu = raw_mu * sqrt(varofmu) + centerofmu;
	vector[g] mu_cond = raw_mu_cond * sqrt(varofcond);
	matrix[g, k] mu_ind = raw_mu_ind * sqrt(varofind);
}

model {
	// hp_r[1] ~ gamma(hp_alpha_r[1], hp_alpha_r[2]);
	// hp_r[2] ~ gamma(hp_beta_r[1], hp_beta_r[2]);
	hp_r ~ uniform(hp_uniform_r[1], hp_uniform_r[2]);
	nb_r ~ gamma(hp_r[1], hp_r[2]);

	// centerofmu follows non-informative prior
	varofmu ~ gamma(hp_varofmu[1], hp_varofmu[2]);
	// implicit mu ~ N(centerofmu, sqrt(varofmu))
	raw_mu ~ std_normal();

	// hp_varofind[1] ~ gamma(hp_alpha_varofind[1], hp_alpha_varofind[2]);
	// hp_varofind[2] ~ gamma(hp_beta_varofind[1], hp_beta_varofind[2]);
	hp_varofind ~ uniform(hp_uniform_varofind[1], hp_uniform_varofind[2]);
	varofind ~ gamma(hp_varofind, hp_varofind);
	// implicit mu_ind ~ N(0.0, sqrt(varofind))
	to_vector(raw_mu_ind) ~ std_normal();

  varofcond ~ gamma(hp_varofcond[1], hp_varofcond[2]);
	// implicit mu_cond ~ N(0.0, sqrt(varofcond))
	to_vector(raw_mu_cond) ~ std_normal();
	matrix[ngene, 2] mu_cond_c2 = append_col(mu_cond, -1 * mu_cond);

	matrix[ngene, ncell] lambda = logs + rep_matrix(mu, ncell)
		+ mu_cond_c2[, cond] + mu_ind[, ind];
	matrix[ngene, ncell] nb_rr = rep_matrix(n_r, ncell);

	// to_vector is column-major order
	// lambda wil in order: g1_c1, g2_c1, ..., g1_c2, g2_c2, ..., g1_cn, g2_cn
	y1d ~ neg_binomial_2_log(to_vector(lambda), to_vector(nb_rr));
}
