library(ComplexHeatmap)

# * generate single-cell gene expression demo data
set.seed(123)
nr1 <- 4
nr2 <- 8
nr3 <- 6
nr <- nr1 + nr2 + nr3
nc1 <- 6
nc2 <- 8
nc3 <- 10
nc <- nc1 + nc2 + nc3
mat <- cbind(rbind(matrix(rnorm(nr1 * nc1, mean = 1, sd = 0.5), nr = nr1),
  matrix(rnorm(nr2 * nc1, mean = 0, sd = 0.5), nr = nr2),
  matrix(rnorm(nr3 * nc1, mean = 0, sd = 0.5), nr = nr3)),
rbind(matrix(rnorm(nr1 * nc2, mean = 0, sd = 0.5), nr = nr1),
  matrix(rnorm(nr2 * nc2, mean = 1, sd = 0.5), nr = nr2),
  matrix(rnorm(nr3 * nc2, mean = 0, sd = 0.5), nr = nr3)),
rbind(matrix(rnorm(nr1 * nc3, mean = 0.5, sd = 0.5), nr = nr1),
  matrix(rnorm(nr2 * nc3, mean = 0.5, sd = 0.5), nr = nr2),
  matrix(rnorm(nr3 * nc3, mean = 1, sd = 0.5), nr = nr3))
)
dmat <- mat[sample(nr, nr), sample(nc, nc)]
rownames(dmat) <- paste0("Gene", seq_len(nr))
colnames(dmat) <- paste0("Cell", seq_len(nc))
gexp.demo <- ComplexHeatmap::Heatmap(
  dmat, rect_gp = gpar(col = "white", lwd = 2),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE)
withr::with_pdf(
  new = "out/demo.heatmap.pdf",
  width = 7, height = 7,
  code = {print(gexp.demo)}
)
