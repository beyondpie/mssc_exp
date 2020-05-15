transformed data {
    int mydata[10] = {0, 1, 2, 0, 1, 0, 1,2,3,2};
}

parameters {
    /* this will work */
    vector[10] y;
    vector<lower=0>[10] sigma;

    /* real y[2]; */
}

transformed parameters {
    vector[10] Mu = sqrt(sigma) .* y;
}

model {
   // below will error
    /* vector[2] z; */

    // this will work.
    real z[2];
    print("no define for z, z is :", z);

    print("Before sampling y, Mu is:", Mu);
    print("Before sampling y, y is:", y);

    /* print("log density before = ", target()); */
    y ~ normal(0,1);
    print("sampling y: ", y);
    print("After sampling y, Mu is:", Mu);

    print("before sampling, sigma is: ", sigma);
    sigma ~ inv_gamma(1,1);
    print("after sampling inv_gamma, sigma is: ", sigma);

    print("After sampling sigma and y, Mu is:",Mu);
    /* print("log density after = ", target()); */

    mydata ~ poisson_log(Mu);

    print("log density after = ", target());

}
