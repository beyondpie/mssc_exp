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
