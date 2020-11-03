// fit scaled negative binomial dist.
// consider the random effect
// fit one gene at a time

// may have numerical problem.

data {
	int n;
	vector<lower=0>[n] s; // scale factor
	int<lower=0> y[n]; // count
	real<lower=0> r; // dispersion/size in nb.
}

transformed data {
	vector[n] logs = log(s);
}

parameters {
	real mu;
}

model {
	// sigma2 follows non informative prior.
	y ~ neg_binomial_2_log(logs + mu, r);
}
