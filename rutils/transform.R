to_onehot_matrix <- function(str_vec) {
  as.matrix(mltools::one_hot(data.table::data.table(as.factor(str_vec))))
}

quickdump <- function(name, myenv = parent.frame()) {
  rstan::stan_rdump(
    c(
      "N", "K", "J", "G", "XCond", "XInd", "S","B", "P"
    ),
    file = name,
    envir = myenv,
    append = TRUE
  )
}
