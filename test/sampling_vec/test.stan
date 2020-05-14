parameters {
    /* this will work */
    vector[2] y;

    /* real y[2]; */
}

model {
   // below will error
    /* vector[2] z; */

    // this will work.
    real z[2];
    print(z);

    print("log density before = ", target());
    y ~ normal(0,1);
    z ~ normal(0,1);
    print("log density after = ", target());
    print(z);
}
