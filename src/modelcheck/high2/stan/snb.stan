// fit scaled negative binomial distribution
// one gene at a time

data {
	int n;
	vector<lower=0>[n] s;
	int<lower=0> y[n];
}

transformed data {
	vector[n] logs = log(s);
}

parameters {
	real mu;
	// r should not be too large
	real<lower=0> r;
}


model {
	// TODO: maybe a better prior.
	r ~ gamma(1.0, 1.0);
	mu ~ normal(0.0, 20);
	y ~ neg_binomial_2_log(logs + mu, r);
}
