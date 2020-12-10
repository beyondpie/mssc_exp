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
  
}


