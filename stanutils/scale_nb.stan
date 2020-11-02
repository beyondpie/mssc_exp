// fit scaled negative binomial dist.
// fit one gene at a time.

data{
	int n;
	vector<lower=100>[n] s; // scale factor
	int<lower=0> y[n]; // count
}

transformed data {
	vector[n] logs = log(s);
}

parameters {
	real mu;
	real<lower=0> r;
}

model {
	// assume non-informative prior for both mu and r;
	y ~ neg_binomial_2_log(logs + mu, r);
}
