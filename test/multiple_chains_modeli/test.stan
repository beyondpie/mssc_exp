transformed data {
    vector[3] y = [-1, 0, 1]';
}

parameters {
    real mu;
}

model {
    mu ~ std_normal();
    y ~ normal(mu, 1);
}
