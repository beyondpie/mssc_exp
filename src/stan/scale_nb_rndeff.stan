// fit scaled negative binomial dist.
// consider the random effect
// fit one gene at a time

data {
	int n;
	vector<lower=100>[n] s; // scale factor
	int<lower=0> y[n]; // count
}

transformed data {
	vector[n] logs = log(s);
}

parameters {
	real raw_mu;
	real<lower=0> r;
	real<lower=0> sigma2;
}

transformed parameters {
	real mu = raw_mu * sqrt(sigma2);
}

model {
	raw_mu ~ std_normal(); // implicit mu ~ normal(0.0, sqrt(sigma2));
	// r and sigma2 follows non informative prior.
	y ~ neg_binomial_2_log(logs + mu, r);
}
