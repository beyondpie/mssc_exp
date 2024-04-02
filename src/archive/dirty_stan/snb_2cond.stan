// fit scaled negative binomial distribution
// one gene at a time

// fix the average mu and r

data {
  int n;
  vector<lower=0>[n] s;
  int<lower=0> y[n];
  int<lower=1, upper=2> cond[n];
  real<lower=0> r;
  real mu;
  // hp of sd for mu_cond
  // recommend using 5 or 10
  real<lower=0> hp_sigma;
}

transformed data {
  vector[n] logs = log(s);
}

parameters {
  real mucond;
}

transformed parameters {
  vector[2] mu_cond_2 = [mucond, -1 * mucond]';
}

model {
  // a relatively smooth prior
  mucond ~ normal(0.0, hp_sigma);
  vector[n] lambda = logs + mu + mu_cond_2[cond];
  y ~ neg_binomial_2_log(lambda, r);
}
