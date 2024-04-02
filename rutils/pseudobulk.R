<<<<<<< HEAD
pseudobulk_deseq2 <- function(cnt_gbc,
                              mybatches,
                              myconds) {
  library(magrittr)
  names(myconds) <- as.character(mybatches)
  ubatches <- unique(mybatches)
  uconds <- myconds[as.character(ubatches)]

  mypseudobulk <- ubatches %>%
    map(.f = function(batch) {
      rowSums(cnt_gbc[, mybatches %in% batch])
    }) %>%
    do.call(what = cbind, args = .)

  colnames(mypseudobulk) <- as.character(ubatches)
  deseqds <- DESeq2::DESeqDataSetFromMatrix(
    countData = mypseudobulk,
    colData = data.frame(uconds),
    design = ~uconds
  )
  deseqres <- DESeq2::DESeq(deseqds) %>%
    DESeq2::results() %>%
    data.frame()
  invisible(deseqres)
}

=======
## functions for pseudobulk analysis

get_pseudobulk <- function(cnt_gbc, mybatches) {
  ## given the cnt data, and colmeta batches,
  ## return the pseudobulk data with the colnames.

  ubatches <- sort(unique(mybatches))
  pseudobulk <- vapply(ubatches, function(i) {
    rowSums(cnt_gbc[ ,mybatches %in% i])
  }, FUN.VALUE = rep(0, nrow(cnt_gbc)))
  colnames(pseudobulk) <- as.character(ubatches)
  invisible(pseudobulk)
}


pseudobulk_deseq2 <- function(cnt_gbc,
                              mybatches,
                              myconds) {
  ## using deseq2 to analyze pseudobulk
  ## return the data.frame format of deseq result.

  mypseudobulk <- get_pseudobulk(cnt_gbc, mybatches)
  names(myconds) <- as.character(mybatches)
  ubatches <- colnames(mypseudobulk)
  uconds <- myconds[as.character(ubatches)]
  
  coldf <- data.frame(ubatches, uconds)
  exp_design <- ~ as.factor(uconds)

  dataset <- DESeq2::DESeqDataSetFromMatrix(
    countData = mypseudobulk,
    colData = coldf, 
    design = exp_design
  )
  r <- data.frame(DESeq2::results(DESeq2::DESeq(dataset)))
  invisible(r)
}


>>>>>>> dev
calc_auc <- function(deseq2_res, degs, ndegs,
                     scorecol = "pvalue") {
  mybackend <- c(rep(TRUE, length(degs)), rep(FALSE, length(ndegs)))
  bgnms <- c(degs, ndegs)
  names(mybackend) <- as.character(bgnms)

  scores <- deseq2_res[[scorecol]]
  names(scores) <- rownames(deseq2_res)
  scores[is.na(scores)] <- 1.0
  myauc <- caTools::colAUC(scores[bgnms], mybackend)
  invisible(list(auc = myauc, sts = scores))
}
