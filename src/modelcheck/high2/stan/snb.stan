// fit scaled negative binomial distribution
// one gene at a time

data {
	int n;
	vector<lower=0>[n] s;
	int<lower=0> y[n];
	vector<lower=0>[2] hpg;
}

transformed data {
	vector[n] logs = log(s);
}

parameters {
	real mu;
	// real<lower=0, upper=100> r;
	real<lower=0> r;
	real mu_u;
	real<lower=0, upper=100> r_u;
}

model {
	// assume non-informative prior
	r ~ gamma(hpg[1], hpg[2]);
	y ~ neg_binomial_2_log(logs + mu, r);
	y ~ neg_binomial_2_log(logs + mu_u, r_u);
}
