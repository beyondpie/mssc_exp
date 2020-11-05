// fit scaled negative binomial dist.
// one gene at a time
// fixed r and mu, estimate the individual effect.

data {
	int n; // sample size
	int k; // num of ind
	vector<lower=0>[n] s; // scale factor
	int<lower=0> y[n]; // cnt
	int<lower=1, upper=k> ind[n]; // individual index
	real mu; // fixed mean
	real<lower=0> r; // fixed dispersion
}

transformed data {
	vector[n] logs = log(s);
}

parameters {
	vector[k] mu_ind;
}

model {
	// assume mu_ind follow non informative prior.
	y ~ neg_binomial_2_log(logs + mu + mu_ind[ind], r);
}
