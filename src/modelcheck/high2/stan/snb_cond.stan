// fit scaled negative binomial distribution
// one gene at a time
// fix r
// no need to consider the condiitons here

data{
  int n;
  vector<lower=0>[n] s; //scale factor
  int<lower=0> y[n];
  real<lower=0> r; // dispersion / size in nb (stan/r)
}

transformed data {
  vector[n] logs = log(s);
}

parameters {
  real mu;
}

model {
  // mu follows non informative prior
  y ~ neg_binomial_2_log(logs + mu, r);
}

