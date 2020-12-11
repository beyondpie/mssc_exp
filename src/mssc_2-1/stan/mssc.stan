// Hierarchical Bayesian based NB model

// This inherits from high2.stan,
// but different when modeling the conditional variances.

// Model the individual effect in gene-wise
// - model mu has a normal prior sharing among genes
//   - centerofmu: non informative prior
//   - varofmu: inv-gamma prior (hp is estimated by data)
// - model r has a log-normal prior sharing among genes
//   - centerofr: non informative prior
//   - varofr: inv-gamma prior (hp is estimated by data)
//
//
// DIFFERENCE WITH HIGH2:
// - model mu_cond (each gene under all the conditions)
//   has a normal prior (0.0, sqrt(varofmucond_g))
//   - varofmucond_g ~ inv-gamma prior among genes
//     hp is estimated by data
//
// 
// - model mu_ind for each gene has a normal prior
//   - centerofmuind[nind] are same given a individual
//       - it then follows: N(0.0, tau)
//           - tau follows: inv-gamma prior among individual (hp is estimated)
//   - varofmuind[g] are same given a individual
//     - it then follows: inv-gamma prior for each individual (hps are estimated)


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
  
  //hps
  vector<lower=0>[2] hp_varofmu;
  vector<lower=0>[2] hp_varofr;
  // No hyper for varofcond;
  vector<lower=0>[2] hp_tau2;
  // since we have lots of genes,
  // setting hp_varofind a little big should be
  // helpful and also no need to set a hyper prior
  // This is different with mucond
  matrix<lower=0>[nind, 2] hp_varofind;
}

transformed data {
  //rep_matrix of row_vector: in row-wise
  // [s1, s2, ..., ; s1, s2,... ]
  matrix[ngene, ncell] logs = rep_matrix(log(s), ngene);
  
  // row-major order
  // [g1_c1, g2_c1, ...., g1_c2, g2_c2, ...]
  int y1d[ncell * ngene] = to_array_1d(y);
}

parameters {
  real centerofmu;
  real<lower=0> varofmu;
  vector[ngene] raw_mu;
  
  real centerofr;
  real<lower=0> varofr;
  vector[ngene] raw_r;

  vector<lower=0>[2] hp_varofcond;
  vector<lower=0>[ngene] varofcond;
  matrix[ngene, ncond] raw_mucond;

  real<lower=0> tau2;
  row_vector[nind] raw_centerofind;
  vector<lower=0>[nind] varofind;
  matrix[ngene, nind] raw_muind;
}

transformed parameters {
  vector[ngene] mu = centerofmu + sqrt(varofmu) * raw_mu;
  vector[ngene] r = centerofr + sqrt(varofr) * raw_r;
  vector[ngene] nb_r = exp(r); // nb_r ~ lognormal(centerofr, sqrt(varofr))
  //matrix product
  // DIFFERENCE with high2
  matrix[ngene, ncond] mucond = diag_matrix(sqrt(varofcond)) * raw_mucond;
  // vector product scalar
  row_vector[nind] centerofind = raw_centerofind * sqrt(tau2);
  // matrix product and add matrix
  matrix[ngene, nind] muind = raw_muind * diag_matrix(sqrt(varofind)) +
    rep_matrix(centerofind, ngene);
}

model {
  // centerofmu ~ non-informative
  varofmu ~ inv_gamma(hp_varofmu[1], hp_varofmu[2]);
  raw_mu ~ std_normal(); // implicit mu ~ N(centerofmu, sqrt(varofmu))

  // centerofr ~ non-informative
  varofr ~ inv_gamma(hp_varofr[1], hp_varofr[2]);
  raw_r ~ std_normal(); // impliccit r ~ N(centerofr, sqrt(varofr))

  // DIFFERENCE with high2
  // hp_varofcond ~ non-informative
  //matrix are column-wise store, so this is efficient
  to_vector(varofcond) ~ inv_gamma(hp_varofcond[1], hp_varofcond[2]);

  tau2 ~ inv_gamma(hp_tau2[1], hp_tau2[2]);
  raw_centerofind ~ std_normal(); //implicit centerofind ~ N(0.0, sqrt(tau2));
  varofind ~ inv_gamma(hp_varofind[ ,1], hp_varofind[ ,2]);
  to_vector(raw_muind) ~ std_normal(); // implicit muind ~ N(centerofind, sqrt(varofind))

  matrix[ngene, ncell] lambda = logs + rep_matrix(mu, ncell) + mucond[, cond] +
    muind[, ind];
  
  matrix[ngene, ncell] nb_rr = rep_matrix(nb_r, ncell);
  // to_vector is column-major order
  // to_vector(lambda) in order: g1_c1, g2_c1, ..., g1_c2, g2_c2, ...
  y1d ~ neg_binomial_2_log(to_vector(lambda), to_vector(nb_rr));
}


