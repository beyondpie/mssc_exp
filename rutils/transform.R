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

to_bagwiff <- function(cnt_gbc, batch, conds, totcntpcell,outf, rdump = FALSE) {
  ## bagwiff: modeling batch effects on gene-wise level
  ncells <- ncol(cnt_gbc)
  if (ncells != length(batch)) {
    error(
      stringr::str_glue(
        "num of cell not match: cnt_gbc({ncells}); batch ({length(batch)})"
      )
    )
  }
  if (ncells != length(conds)) {
    error(
      stringr::str_glue(
        "num of cell not match: cnt_gbc({ncells}); conds ({length(conds)})"
      )
    )
  }
  Xcg <- t(as.matrix(cnt_gbc))
  XInd <- to_onehot_matrix(batch)
  XCond <- to_onehot_matrix(conds)
  N <- nrow(XCond)
  J <- ncol(XCond)
  K <- ncol(XInd)
  G <- ncol(Xcg)

  ## S <- rowSums(Xcg)
  S <- totcntpcell

  ## bagmiff mdel
  ## bagmiff: modeling batch effects on gene-module level
  ## add gene module infomration.
  ## a trivial one
  P <- 1L
  B <- matrix(1:G, nrow = G, ncol = P)
  if (rdump) {
    rstan::stan_rdump(c(
      "N", "J", "K", "G", "S", "P", "B",
      "XCond", "XInd", "Xcg"
    ), file = outf)
  }
  else {
    Xcg <- as.data.frame(Xcg)
    XInd <- as.data.frame(XInd)
    XCond <- as.data.frame(XCond)
    save(N, J, K, G, S, P, B, XInd, XCond, Xcg, file = outf)
  }
}

subsampling <- function(myarray, size, replace = FALSE) {
  if (length(myarray) <= size) {
    return(myarray)
  } else {
    return(sample(myarray, size = size, replace = replace))
  }
}

stat_geneset <- function(pool, geneset) {
  ovlp <- intersect(pool, geneset)
  message(stringr::str_glue("num of geneset: {length(geneset)}"))
  message(stringr::str_glue("num in pool: {length(ovlp)}"))
  return(ovlp)
}
