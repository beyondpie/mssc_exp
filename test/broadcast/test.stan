transformed data {
    int G = 10;
    int N = 5;
}

parameters {
    vector[G] mu;
}

transformed parameters{
    /* vector[N] sigma = dot_product(mu[1:N], mu[1:N]); */
    vector[N] sigma = mu[1:N] .* mu[1:N];
}

model {
    mu ~ normal(0,2);
    print("mu value are: ", mu);
    print("sigma is:", sigma);
}
