## summarize results from symsim v2
library(tidyverse)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(ggpubr)
library(here)

load_symsim_result <- function(nind_per_cond = 3, brn_len = 0.5, bimod = 1,
                               sigma = 0.4,
                               alpha = 0.2,
                               ngene = 200,
                               ncell = 300,
                               result_dir = here::here("src", "symsim", "symsim_20210528")) {
  dea_dir <- file.path(result_dir, "dea")
  nind_all <- nind_per_cond * 2
  symsim_prefix <- str_glue(
    "{ngene}gene", "{nind_all}ind",
    "{ncell}cell", "{brn_len}w", "{bimod}bimod",
    "{sigma}sigma", "{alpha}alpha", .sep = "_")
  return(invisible(readRDS(
    file = file.path(dea_dir, str_glue("dea_auc_{symsim_prefix}.rds")))))
}

f_across_rpt <- function(r, f = mean) {
  ## take f along the last dim of r, which is the repeat number.
  ## r is a three-dim array.
  ## f is designed to be mean or sd.
  ## return: matrix: dim(r)[1] by dim(r)[2]

  nrow <- dim(r)[1]
  ncol <- dim(r)[2]
  rownames <- dimnames(r)[[1]]
  colnames <- dimnames(r)[[2]]
  ## ignore the fitting error point
  r[r < 0.1] <- NA
  fr <- matrix(0.0, nrow = nrow, ncol = ncol,
    dimnames = list(rownames, colnames))
  for (i in 1:nrow) {
    for (j in 1:ncol) {
      fr[i, j]  <- f(r[i, j, ], na.rm = TRUE)
    }
  }
  invisible(fr)
}

scatter_plot <- function(r, nind = 3, expnm = NULL,
                         save_fig = FALSE,
                         save_dir = here::here("src", "symsim")) {
  ## r, a matrix, row:methods, col: cells
  ## rows: vector of string, choose the rows in r
  ## nind and use_groupshift: for title only
  ## return:
  ## - ggplot figure of geom points

  data <- as.data.frame(r)
  data$methods <- rownames(data)
  gr <- gather(data = data, key = "num_of_cell_per_ind", value = "AUC",
    matches("[0-9][0-9]"))
  gr$num_of_cell_per_ind <- as.numeric(gr$num_of_cell_per_ind)
  p <- ggplot(gr, aes(x = factor(num_of_cell_per_ind), y = AUC)) +
    geom_point(aes(colour = methods, shape = methods),
      size = 6, alpha = 0.8) +
    xlab("Number of cells / individual") +
    theme_bw(base_size = 15) +
    scale_color_brewer(palette = "Dark2") +
    ggtitle(stringr::str_glue("{expnm} \n{nind} Individuals / Condition")) +
    theme(plot.title = element_text(size = 15, hjust = 0.5))
  if (save_fig) {
    ggsave(plot = p,
      filename = file.path(save_dir,
        stringr::str_glue("{expnm}_symsim_{nind}inds.pdf")),
      device = "pdf",
      width = 10,
      height = 7
    )
  }
  invisible(p)
}



## Result
methods <- c("mssc_VI", "mssc_MAP", "GLM", "Deseq2", "Wilcox", "t")
ngene <- 200
ncell <- 300
brn_len <- 0.5
bimod <- 1

## when nind_per_cond is 3
nind_per_cond <- 3
## - 200gene_6ind_0.5w_1bimod_0.4sigma_0.2alpha
## - 200gene_6ind_0.5w_1bimod_0.4sigma_0.1alpha
## - 200gene_6ind_0.5w_1bimod_0.6sigma_0.2alpha

### when alpha is 0.2
sigma <- 0.4
alpha <- 0.2
r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
  bimod = bimod, sigma = sigma, alpha = alpha)
r_mean <- f_across_rpt(r, f = mean)
rownames(r_mean) <- methods
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = 3, expnm = expnm, save_fig = TRUE)

### when alpha is 0.1
alpha <- 0.1
r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
  bimod = bimod, sigma = sigma, alpha = alpha)
r_mean <- f_across_rpt(r, f = mean)
rownames(r_mean) <- methods
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = 3, expnm = expnm, save_fig = TRUE)

### when sigma is 0.6, alpha = 0.2
alpha <- 0.2
sigma <- 0.6
r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
  bimod = bimod, sigma = sigma, alpha = alpha)
r_mean <- f_across_rpt(r, f = mean)
rownames(r_mean) <- methods
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = 3, expnm = expnm, save_fig = TRUE)
### when sigma is 0.6, alpha = 0.1
nind_per_cond <- 3
alpha <- 0.1
sigma <- 0.6
r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
  bimod = bimod, sigma = sigma, alpha = alpha)
r_mean <- f_across_rpt(r, f = mean)
rownames(r_mean) <- methods
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = nind_per_cond, expnm = expnm, save_fig = TRUE)


## when nind_per_cond is 5
nind_per_cond <- 5
### sigma = 0.4; alpha = 0.1
alpha <- 0.1
sigma <- 0.4
r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
  bimod = bimod, sigma = sigma, alpha = alpha)
r_mean <- f_across_rpt(r, f = mean)
rownames(r_mean) <- methods
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = nind_per_cond, expnm = expnm, save_fig = TRUE)

### sigma = 0.4; alpha = 0.2
alpha <- 0.2
sigma <- 0.4
r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
  bimod = bimod, sigma = sigma, alpha = alpha)
r_mean <- f_across_rpt(r, f = mean)
rownames(r_mean) <- methods
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = nind_per_cond, expnm = expnm, save_fig = TRUE)

### sigma = 0.6; alpha = 0.1
nind_per_cond <- 5
sigma <- 0.6
alpha <- 0.1

r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
  bimod = bimod, sigma = sigma, alpha = alpha)
r_mean <- f_across_rpt(r, f = mean)
rownames(r_mean) <- methods
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = nind_per_cond, expnm = expnm, save_fig = TRUE)

### sigma = 0.6; alpha = 0.2
nind_per_cond <- 5
sigma <- 0.6
alpha <- 0.2

r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
  bimod = bimod, sigma = sigma, alpha = alpha)
r_mean <- f_across_rpt(r, f = mean)
rownames(r_mean) <- methods
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = nind_per_cond, expnm = expnm, save_fig = TRUE)

## when nind_per_cond is 10, and methods is opt
methods <- c("mssc_VI", "mssc_MAP", "GLM", "Deseq2", "Wilcox", "t")
ngene <- 200
ncell <- 300
brn_len <- 0.5
bimod <- 1
nind_per_cond <- 10
### sigma = 0.6, alpha = 0.1
sigma <- 0.6
alpha <- 0.1
r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
                        bimod = bimod, sigma = sigma, alpha = alpha,
                        result_dir = here::here("src", "symsim", "opt_symsim_20210528"))
r_mean <- f_across_rpt(r, f = mean)
r_mean <- r_mean[2:length(methods), ]
rownames(r_mean) <- methods[2:length(methods)]
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = nind_per_cond, expnm = expnm, save_fig = TRUE)
### sigma = 0.6, alpha = 0.2
sigma <- 0.6
alpha <- 0.2
r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
                        bimod = bimod, sigma = sigma, alpha = alpha,
                        result_dir = here::here("src", "symsim", "opt_symsim_20210528"))
r_mean <- f_across_rpt(r, f = mean)
r_mean <- r_mean[2:length(methods), ]
rownames(r_mean) <- methods[2:length(methods)]
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = nind_per_cond, expnm = expnm, save_fig = TRUE)

### sigma = 0.4, alpha = 0.2
sigma <- 0.4
alpha <- 0.2
r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
                        bimod = bimod, sigma = sigma, alpha = alpha,
                        result_dir = here::here("src", "symsim", "opt_symsim_20210528"))
r_mean <- f_across_rpt(r, f = mean)
r_mean <- r_mean[2:length(methods), ]
rownames(r_mean) <- methods[2:length(methods)]
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = nind_per_cond, expnm = expnm, save_fig = TRUE)

### sigma = 0.4, alpha = 0.1
sigma <- 0.4
alpha <- 0.1
r <- load_symsim_result(nind_per_cond = nind_per_cond, brn_len = brn_len,
                        bimod = bimod, sigma = sigma, alpha = alpha,
                        result_dir = here::here("src", "symsim", "opt_symsim_20210528"))
r_mean <- f_across_rpt(r, f = mean)
r_mean <- r_mean[2:length(methods), ]
rownames(r_mean) <- methods[2:length(methods)]
nind_all <- nind_per_cond * 2
expnm <- str_glue(
  "{ngene}gene", "{nind_all}ind",
  "{ncell}cell", "{brn_len}w", "{bimod}bimod",
  "{sigma}sigma", "{alpha}alpha", .sep = "_")
p <- scatter_plot(r = r_mean, nind = nind_per_cond, expnm = expnm, save_fig = TRUE)

