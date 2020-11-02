// fit scaled negative binomial dist.
// consider the random effect
// fit one gene at a time

// may have numerical problem.

data {
	int n;
	vector<lower=100>[n] s; // scale factor
	int<lower=0> y[n]; // count
}

transformed data {
	vector[n] logs = log(s);
}

parameters {
	real mu;
	real lambda;
	real<lower=0> r;
	real<lower=0> sigma2;
}

model {
	// r and sigma2 follows non informative prior.
	lambda ~ normal(logs + mu, sqrt(sigma2));
	y ~ neg_binomial_2_log(lambda, r);
}
