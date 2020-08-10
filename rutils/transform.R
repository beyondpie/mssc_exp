
umi_scale_factor <- 10000
fullength_scale_factor <- 1e+06

myggtitle <- theme(plot.title = element_text(size = 15, hjust = 0.5))

getcpm <- function(cnt, scale_level = umi_scale_factor) {
  ## * cnt is gene by cell matrix
  return(scale(cnt,
    center = F,
    scale = colSums(cnt)
  ) * scale_level)
}

getlogcpm <- function(cpm, margin = 1) {
  return(log(cpm + margin))
}

gettsne <- function(logcpm, seed = 0L) {
  set.seed(seed)
  tsne <- Rtsne::Rtsne(t(logcpm),
    pca_scale = T, pca_center = T,
    initial_dims = 50, pca = T,
    check_duplicates = F
  )
  return(tsne)
}

saveggplot <- function(p, fnm = NULL) {
  if (!is.null(fnm)) {
    message(stringr::str_glue("saving figute to {fnm}"))
    pdf(file = fnm)
    plot(p)
    dev.off()
  }
}
plottsne <- function(tsneembed, celltypes,
                     title = "tsne", fnm = NULL) {
  ## tsneembed should be tsneres$Y
  library(ggplot2)
  cellnum <- nrow(tsneembed)
  ## avoid treating celltypes as continuous variables
  celltypes <- as.character(celltypes)
  typenum <- length(celltypes)

  if (cellnum != typenum) {
    stop(stringr::str_glue(
      "dim inconsistant: ",
      "tsne rownum ({cellnum});",
      "cell type num ({typenum})"
    ))
  }
  dat <- data.frame(cbind(tsneembed), celltypes)
  colnames(dat) <- c("TSNE1", "TSNE2", "Type")
  p <- ggplot(dat, aes(x = TSNE1, y = TSNE2, color = Type)) +
    geom_point() +
    theme_bw() +
    ggtitle(title) +
    myggtitle
  saveggplot(p, fnm)
  return(p)
}



plothist <- function(cnt, mybinwidth = 200,
                     max_xlim = 2e+05, mytitle = "histogram", fnm = NULL) {
  library(ggplot2)
  tmp <- colSums(cnt)
  mydata <- data.frame(cell_library_size = tmp)
  p <- ggplot(data = mydata, aes(x = cell_library_size)) +
    geom_histogram(
      binwidth = mybinwidth
    ) +
    xlim(0, max_xlim) +
    ggtitle(mytitle) +
    myggtitle

  saveggplot(p, fnm)
  return(p)
}

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

to_bagwiff <- function(cnt_gbc, batch, conds,
                       totcntpcell, outf, rdump = FALSE) {
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
    invisible(rstan::stan_rdump(c(
      "N", "J", "K", "G", "S", "P", "B",
      "XCond", "XInd", "Xcg"
    ), file = outf))
  }
  else {
    Xcg <- as.data.frame(Xcg)
    XInd <- as.data.frame(XInd)
    XCond <- as.data.frame(XCond)
    invisible(save(N, J, K, G, S, P, B, XInd, XCond, Xcg, file = outf))
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

draw_invgamma <- function(alpha, rate = 1 / alpha, from = 0.0001, to = 20) {
  ggplot2::theme_set(ggplot2::theme_bw())
  x <- seq(from, to, 0.01)
  ggplot2::qplot(x, invgamma::dinvgamma(x, alpha, rate), geom = "line")
}

## point relative to the region
pntrela2rgn <- function(rgn = c(0.0, 1.0), pnt = 0.0) {
  post <- 0
  if (pnt < rgn[1]) {
    post <- -1
  } else if (pnt > rgn[2]) {
    post <- 1
  } else {
    post <- 0
  }
  return(post)
}

## format float
fmtflt <- function(f, nsmall = 3) {
  return(format(round(f, nsmall), nsmall = nsmall))
}

## * utils for gettign rstan results
load_stan_vi <- function(path) {
  rstan::read_stan_csv(path)
}

load_stan_mc <- function(dirpath, modelnm) {
  csvfiles <- dir(
    path = dirpath, pattern = paste0(modelnm, "[0-9].csv"),
    full.names = T
  )
  rstan::read_stan_csv(csvfiles)
}

load_stan <- function(dirnm, modelnm = "v1-1", method = "vi",
                      vi_dir = "vi", mc_dir = "mc") {
  if (method == "vi") {
    vifnm <- paste0(modelnm, ".csv")
    mystanfit <- load_stan_vi(paste(dirnm,
      vi_dir, vifnm,
      sep = "/"
    ))
  }
  if (method == "mc") {
    mcprefix <- paste0(modelnm, "_chain_")
    mystanfit <- load_stan_mc(
      dirpath = paste(dirnm, mc_dir, sep = "/"),
      modelnm = mcprefix
    )
  }
  return(mystanfit)
}

## TODO: add explanation about which is case and control
## ctrlmnscase: control minus case
get_ctrlmnscase_par <- function(mystanfit, par = "MuCond") {
  mus <- rstan::extract(mystanfit, pars = par)[[par]]
  delta <- as.data.frame(mus[, , 1] - mus[, , 2])
  return(delta)
}

## simple t statistics
calt <- function(delta, fn = matrixStats::colMedians) {
  fnhat <- fn(as.matrix(delta))
  std_hat <- matrixStats::colSds(as.matrix(delta) + 1e-10)
  sts <- fnhat / (sqrt(nrow(delta)) * std_hat)
  names(sts) <- colnames(delta)
  return(sts)
}
