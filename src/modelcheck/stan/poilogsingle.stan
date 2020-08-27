data{
	int<lower=0> N;
	int y[N];
	int t[N];
}

parameters {
	real mu;
	real<lower=0> sigma;
	real log_lambda;
}

model {
	sigma ~ cauchy(0, 5);
	mu ~ normal(0, 10);
	log_lambda ~ normal(mu, sigma);
	target += poisson_log_lpmf(y | t * log_lambda);
}
