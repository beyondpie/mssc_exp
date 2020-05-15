#include  myinclude.stan
data {
    real<lower=0> sigma;
}

transformed data {
    real<lower=0> s = sigma + ss;
}

parameters {
    real mu;
}

model {
    y ~ normal(mu, s);
}
