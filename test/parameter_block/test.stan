transformed data {
    vector[3] y = [1, -1, 0]';
    int N = 10000;
}

parameters {
    real<lower=0> sigma[N];
}

model {
    for (i in 1:N) {
        if(sigma[i] < 0) {
            print("sampling sigma is ", sigma[i]);
        }
    }
    y ~ normal(0, sigma[1]);
}
