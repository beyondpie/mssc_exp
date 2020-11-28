// fit scaled negative binomial distribution
// one gene at a time

data {
  int n;
  vector<lower=0>[n] s;
  int<lower=0> y[n];
  // hp for r
  // suggest: 0.05, 0.05
  vector<lower=0>[2] hpg;
}

transformed data {
  vector[n] logs = log(s);
}

parameters {
  real mu;
  real<lower=0> r;
}

model {
  // assume non-informative prior for mu
  r ~ gamma(hpg[1], hpg[2]);
  y ~ neg_binomial_2_log(logs + mu, r);
}
