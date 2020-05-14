transformed data {
    vector[3] y = [1, -1, 0]';
    int N = 10000;
}

parameters {
    /* real sigma[N]; */
    /* real sigma; */
    real<lower=0> sigma;
}

transformed parameters {
    real<lower=0> p = sigma - 10;
}

model {
    sigma ~ inv_gamma(1,1);
    /* for (i in 1:n) { */
    /*     if(sigma[i] < 0) { */
    /*         print("sampling sigma is ", sigma[i]); */
    /*     } */
    /* } */
    if (sigma <0 ) {
        print("sampling sigma is ", sigma);
    }

    if (p < 0 ) {
        print("p is out of control, ", p);
    }

    y ~ normal(0, sigma);
}
