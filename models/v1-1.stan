// Model batch effect shared for all the genes

#include mydata.stan

transformed data {
}

parameters {
    real Sigma2G;
    real Mu;
    vector[G] MuGRaw;
    matrix[G, K] MuIndRaw;
    matrix[G, J] MuCondRaw;
}

transformed parameters {
    vector[G] MuG = Mu + sqrt(Sigma2G) * MuGRaw;
    matrix[]
}


model {
    Mu ~ normal(0.0, sigmaMu);
    Sigma2G ~ inv_gamma(alphaSigma2G, betaSigma2G);
    MuGRaw ~ std_normal(); // implicit MuG ~ normal(Mu, sqrt(Sigma2G))
}
