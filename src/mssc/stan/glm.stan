// Generalised linear model for differential expression analysis.
// Model the gene-wise individual effect.

data {
  int ncell;
  int nind;
  int ngene;
  int ncond;

  row_vector<lower=0>[ncell] s; // scaling factor
  int<lower=1, upper=ncond> cond[ncell];
  int<lower=1, upper=nind> ind[ncell];
  // note the dim
  int<lower=0> y[ncell, ngene];
}


transformed data {
  //rep_matrix of row_vector: in row-wise
  // [s1, s2, ..., ; s1, s2,... ]
  matrix[ngene, ncell] logs = rep_matrix(log(s), ngene);
  
  // row-major order
  // [g1_c1, g2_c1, ...., g1_c2, g2_c2, ...]
  int y1d[ngene * ncell] = to_array_1d(y);
}

parameters {
  vector[ngene] mu;
  vector[ngene] r;
  matrix[ngene, ncond] mucond;
  matrix[ngene, nind] muind;
}

transformed parameters {
  vector[ngene] nb_r = exp(r);
}

model {
  matrix[ngene, ncell] lambda = logs + rep_matrix(mu, ncell) + mucond[, cond] + muind[, ind];
  matrix[ngene, ncell] nb_rr = rep_matrix(nb_r, ncell);
  // to_vector is column-major order
  // to_vector(lambda) in order: g1_c1, g2_c1. g3_c1, ..., g1_c2, g2_c2, ...
  y1d ~ neg_binomial_2_log(to_vector(lambda), to_vector(nb_rr));
}
