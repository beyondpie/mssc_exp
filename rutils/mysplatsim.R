splatSimulate <- function(params = newSplatParams(),
                          method = c("single", "groups", "paths"),
                          verbose = TRUE, ...) {

    checkmate::assertClass(params, "SplatParams")

    method <- match.arg(method)

    if (verbose) {message("Getting parameters...")}
    params <- setParams(params, ...)
    params <- expandParams(params)
    validObject(params)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    nBatches <- getParam(params, "nBatches")
    batch.cells <- getParam(params, "batchCells")
    nGroups <- getParam(params, "nGroups")
    group.prob <- getParam(params, "group.prob")

    if (nGroups == 1 && method == "groups") {
        warning("nGroups is 1, switching to single mode")
        method <- "single"
    }

    if (verbose) {message("Creating simulation object...")}
    # Set up name vectors
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))
    batch.names <- paste0("Batch", seq_len(nBatches))
    if (method == "groups") {
        group.names <- paste0("Group", seq_len(nGroups))
    } else if (method == "paths") {
        group.names <- paste0("Path", seq_len(nGroups))
    }

    # Create SingleCellExperiment to store simulation
    cells <-  data.frame(Cell = cell.names)
    rownames(cells) <- cell.names
    features <- data.frame(Gene = gene.names)
    rownames(features) <- gene.names
    sim <- SingleCellExperiment(rowData = features, colData = cells,
                                metadata = list(Params = params))

    # Make batches vector which is the index of param$batchCells repeated
    # params$batchCells[index] times
    batches <- lapply(seq_len(nBatches), function(i, b) {rep(i, b[i])},
                      b = batch.cells)
    batches <- unlist(batches)
    colData(sim)$Batch <- batch.names[batches]

    if (method != "single") {
        groups <- sample(seq_len(nGroups), nCells, prob = group.prob,
                         replace = TRUE)
        colData(sim)$Group <- factor(group.names[groups], levels = group.names)
    }

    if (verbose) {message("Simulating library sizes...")}
    sim <- splatSimLibSizes(sim, params)
    if (verbose) {message("Simulating gene means...")}
    sim <- splatSimGeneMeans(sim, params)
    if (nBatches > 1) {
        if (verbose) {message("Simulating batch effects...")}
        sim <- splatSimBatchEffects(sim, params)
    }
    sim <- splatSimBatchCellMeans(sim, params)
    if (method == "single") {
        sim <- splatSimSingleCellMeans(sim, params)
    } else if (method == "groups") {
        if (verbose) {message("Simulating group DE...")}
        sim <- splatSimGroupDE(sim, params)
        if (verbose) {message("Simulating cell means...")}
        sim <- splatSimGroupCellMeans(sim, params)
    } else {
        if (verbose) {message("Simulating path endpoints...")}
        sim <- splatSimPathDE(sim, params)
        if (verbose) {message("Simulating path steps...")}
        sim <- splatSimPathCellMeans(sim, params)
    }
    if (verbose) {message("Simulating BCV...")}
    sim <- splatSimBCVMeans(sim, params)
    if (verbose) {message("Simulating counts...")}
    sim <- splatSimTrueCounts(sim, params)
    if (verbose) {message("Simulating dropout (if needed)...")}
    sim <- splatSimDropout(sim, params)

    if (verbose) {message("Done!")}
    return(sim)
}


splatSimLibSizes <- function(sim, params) {

    nCells <- getParam(params, "nCells")
    lib.loc <- getParam(params, "lib.loc")
    lib.scale <- getParam(params, "lib.scale")
    lib.norm <- getParam(params, "lib.norm")

    if (lib.norm) {
        exp.lib.sizes <- rnorm(nCells, lib.loc, lib.scale)
        min.lib <- min(exp.lib.sizes[exp.lib.sizes > 0])
        exp.lib.sizes[exp.lib.sizes < 0] <- min.lib / 2
    } else {
        exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
    }

    colData(sim)$ExpLibSize <- exp.lib.sizes

    return(sim)
}

splatSimGeneMeans <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    out.prob <- getParam(params, "out.prob")
    out.facLoc <- getParam(params, "out.facLoc")
    out.facScale <- getParam(params, "out.facScale")

    # Simulate base gene means
    base.means.gene <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)

    # Add expression outliers
    outlier.facs <- getLNormFactors(nGenes, out.prob, 0, out.facLoc,
                                    out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]

    rowData(sim)$BaseGeneMean <- base.means.gene
    rowData(sim)$OutlierFactor <- outlier.facs
    rowData(sim)$GeneMean <- means.gene

    return(sim)
}

splatSimBatchEffects <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    nBatches <- getParam(params, "nBatches")
    batch.facLoc <- getParam(params, "batch.facLoc")
    batch.facScale <- getParam(params, "batch.facScale")
    means.gene <- rowData(sim)$GeneMean

    for (idx in seq_len(nBatches)) {
        batch.facs <- getLNormFactors(nGenes, 1, 0.5, batch.facLoc[idx],
                                        batch.facScale[idx])
        batch.means.gene <- means.gene * batch.facs
        rowData(sim)[[paste0("BatchFacBatch", idx)]] <- batch.facs
    }

    return(sim)
}

splatSimBatchCellMeans <- function(sim, params) {

    nBatches <- getParam(params, "nBatches")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    gene.means <- rowData(sim)$GeneMean

    if (nBatches > 1) {
        batches <- colData(sim)$Batch
        batch.names <- unique(batches)

        batch.facs.gene <- rowData(sim)[, paste0("BatchFac", batch.names)]
        batch.facs.cell <- as.matrix(batch.facs.gene[,
                                                  as.numeric(factor(batches))])
    } else {
        nCells <- getParam(params, "nCells")
        nGenes <- getParam(params, "nGenes")

        batch.facs.cell <- matrix(1, ncol = nCells, nrow = nGenes)
    }

    batch.means.cell <- batch.facs.cell * gene.means

    colnames(batch.means.cell) <- cell.names
    rownames(batch.means.cell) <- gene.names
    assays(sim)$BatchCellMeans <- batch.means.cell

    return(sim)
}

#' @importFrom SummarizedExperiment rowData
splatSimGroupDE <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    nGroups <- getParam(params, "nGroups")
    de.prob <- getParam(params, "de.prob")
    de.downProb <- getParam(params, "de.downProb")
    de.facLoc <- getParam(params, "de.facLoc")
    de.facScale <- getParam(params, "de.facScale")
    means.gene <- rowData(sim)$GeneMean

    for (idx in seq_len(nGroups)) {
        de.facs <- getLNormFactors(nGenes, de.prob[idx], de.downProb[idx],
                                   de.facLoc[idx], de.facScale[idx])
        group.means.gene <- means.gene * de.facs
        rowData(sim)[[paste0("DEFacGroup", idx)]] <- de.facs
    }

    return(sim)
}


splatSimSingleCellMeans <- function(sim, params) {

    nCells <- getParam(params, "nCells")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    exp.lib.sizes <- colData(sim)$ExpLibSize
    batch.means.cell <- assays(sim)$BatchCellMeans

    cell.means.gene <- batch.means.cell
    cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)

    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names
    assays(sim)$BaseCellMeans <- base.means.cell

    return(sim)
}

#' @importFrom SummarizedExperiment rowData colData assays assays<-
splatSimGroupCellMeans <- function(sim, params) {

    nGroups <- getParam(params, "nGroups")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    groups <- colData(sim)$Group
    group.names <- levels(groups)
    exp.lib.sizes <- colData(sim)$ExpLibSize
    batch.means.cell <- assays(sim)$BatchCellMeans

    group.facs.gene <- rowData(sim)[, paste0("DEFac", group.names)]
    cell.facs.gene <- as.matrix(group.facs.gene[, paste0("DEFac", groups)])
    cell.means.gene <- batch.means.cell * cell.facs.gene
    cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)

    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names
    assays(sim)$BaseCellMeans <- base.means.cell

    return(sim)
}

#' Simulate BCV means
#'
#' Simulate means for each gene in each cell that are adjusted to follow a
#' mean-variance trend using Biological Coefficient of Variation taken from
#' and inverse gamma distribution.
#'
#' @param sim SingleCellExperiment to add BCV means to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated BCV means.
#'
#' @importFrom SummarizedExperiment rowData colData assays assays<-
#' @importFrom stats rchisq rgamma
splatSimBCVMeans <- function(sim, params) {

    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    bcv.common <- getParam(params, "bcv.common")
    bcv.df <- getParam(params, "bcv.df")
    base.means.cell <- assays(sim)$BaseCellMeans

    if (is.finite(bcv.df)) {
        bcv <- (bcv.common + (1 / sqrt(base.means.cell))) *
            sqrt(bcv.df / rchisq(nGenes, df = bcv.df))
    } else {
        warning("'bcv.df' is infinite. This parameter will be ignored.")
        bcv <- (bcv.common + (1 / sqrt(base.means.cell)))
    }

    means.cell <- matrix(rgamma(
        as.numeric(nGenes) * as.numeric(nCells),
        shape = 1 / (bcv ^ 2), scale = base.means.cell * (bcv ^ 2)),
    nrow = nGenes, ncol = nCells)

    colnames(means.cell) <- cell.names
    rownames(means.cell) <- gene.names

    assays(sim)$BCV <- bcv
    assays(sim)$CellMeans <- means.cell

    return(sim)
}

#' Simulate true counts
#'
#' Simulate a true counts matrix. Counts are simulated from a poisson
#' distribution where Each gene in each cell has it's own mean based on the
#' group (or path position), expected library size and BCV.
#'
#' @param sim SingleCellExperiment to add true counts to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated true counts.
#'
#' @importFrom SummarizedExperiment rowData colData assays assays<-
#' @importFrom stats rpois
splatSimTrueCounts <- function(sim, params) {

    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    cell.means <- assays(sim)$CellMeans

    true.counts <- matrix(rpois(
        as.numeric(nGenes) * as.numeric(nCells),
        lambda = cell.means),
    nrow = nGenes, ncol = nCells)

    colnames(true.counts) <- cell.names
    rownames(true.counts) <- gene.names

    assays(sim)$TrueCounts <- true.counts

    return(sim)
}

#' @importFrom stats rbinom rlnorm
getLNormFactors <- function(n.facs, sel.prob, neg.prob, fac.loc, fac.scale) {

    is.selected <- as.logical(rbinom(n.facs, 1, sel.prob))
    n.selected <- sum(is.selected)
    dir.selected <- (-1) ^ rbinom(n.selected, 1, neg.prob)
    facs.selected <- rlnorm(n.selected, fac.loc, fac.scale)
    # Reverse directions for factors that are less than one
    dir.selected[facs.selected < 1] <- -1 * dir.selected[facs.selected < 1]
    factors <- rep(1, n.facs)
    factors[is.selected] <- facs.selected ^ dir.selected

    return(factors)
}

