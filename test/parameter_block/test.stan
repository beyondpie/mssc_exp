transformed data {
    vector[3] y = [1, -1, 0]';
}

parameters {
    real<lower=0> sigma;
}

model {
    print("sampling sigma is ", sigma);
    y ~ normal(0, sigma);
}
