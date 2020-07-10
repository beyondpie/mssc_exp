to_onehot_matrix <- function(str_vec) {
  as.matrix(mltools::one_hot(data.table::data.table(as.factor(str_vec))))
}

quickdump <- function(name, myenv = parent.frame()) {
  rstan::stan_rdump(
    c(
      "N", "K", "J", "G", "XCond", "XInd", "S", "B", "P"
    ),
    file = name,
    envir = myenv,
    append = TRUE
  )
}

print_sc <- function(nr, nc, row = "gene", plat = "scRNAseq") {
  if (row == "gene") {
    message(stringr::str_glue("{plat}: {nr} genes and {nc} cells."))
  } else {
    message(stringr::str_glue("{plat}: {nc} genes and {nr} cells."))
  }
}

rm_mt <- function(seqdata) {
  mts <- grep(pattern = "^MT-", x = rownames(seqdata), value = FALSE)
  nmts <- length(mts)
  if (nmts > 1) {
    message(stringr::str_glue("num of MT genes: {nmts}"))
    mtratios <- colSums(seqdata[mts, ]) / colSums(seqdata)
    ## num of high ratio of mitochotria reads
    nhmt <- length(which(mtratios > 0.2))
    if (nhmt > 0) {
      message(stringr::str_glue(
        "num of cells with at least 20% mt reads: {nhmt}"
      ))
      ## if remove mt: seqdata <- seqdata[, mtratios <= 0.2]
    }
    seqdata <- seqdata[-mts, ]
  }
  return(seqdata)
}
