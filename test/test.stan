model {
    // below will error
    /* vector[2] y; */

    // this will work.
    real y[2];

    print("log density before = ", target());
    y ~ normal(0,1);
    print("log density after = ", target());
    print(y);
}
