transformed data {
    vector[3] y = [1, -1, 0]';
    int N = 10000;
}

parameters {
    /* real sigma[N]; */
    real sigma;
}

model {
    sigma ~ inv_gamma(1,1);
    /* for (i in 1:N) { */
    /*     if(sigma[i] < 0) { */
    /*         print("sampling sigma is ", sigma[i]); */
    /*     } */
    /* } */
    if (sigma <0 ) {
        print("sampling sigma is ", sigma);
    }

    y ~ normal(0, sigma[1]);
}
