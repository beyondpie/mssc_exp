library(tidyverse)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(ggpubr)

load_symsim_result <- function(nind, nindeff, use_groupshift = FALSE,
                               ratio_ind2cond = 0.3,
                               scale_in_diffg = 1.0,
                               scale_in_nondiffg = 1.0,
                               ngene = 80,
                               file_path = here::here("src", "symsim", "results")) {
  ## load symsim simulation results
  ## nind is the number of individuals in on condition
  expnm <- ifelse(use_groupshift, "groupshift", "nonzeroshift")
  fnm <- stringr::str_glue(
    "{ngene}gene", "{nind*2}ind", "{nindeff}indeff", expnm,
    ratio_ind2cond, scale_in_diffg, "{scale_in_nondiffg}.rds", .sep = "_")
  invisible(readRDS(file = file.path(file_path, fnm)))
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
  
  fr <- matrix(0.0, nrow = nrow, ncol = ncol,
                   dimnames = list(rownames, colnames))
  for (i in 1:nrow)  {
    for (j in 1:ncol) {
      fr[i, j]  <- f(r[i,j, ], na.rm = TRUE)
    }
  }
  invisible(fr)
}

scatter_plot <- function(r, rows, nind, use_groupshift,
                         save_fig = TRUE,
                         save_dir = here::here("src", "symsim")) {
  ## r, a matrix, row:methods, col: cells
  ## rows: vector of string, choose the rows in r
  ## nind and use_groupshift: for title only
  ## return:
  ## - ggplot figure of geom points
  
  data <- as.data.frame(r[rows, ])
  data$methods <- rownames(data)
  gr <- gather(data = data, key = "num_of_cell_per_ind", value = "AUC",
               matches("[0-9][0-9]"))
  gr$num_of_cell_per_ind <- as.numeric(gr$num_of_cell_per_ind)
  expnm <- ifelse(use_groupshift, "Group-shift", "Nonzero-shift")
  p <- ggplot(gr, aes(x = factor(num_of_cell_per_ind), y = AUC)) +
    geom_point(aes(colour = methods, shape = methods),
               size = 6, alpha = 0.8) +
    xlab("Number of cells / individual") +
    theme_bw(base_size = 15) + 
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(15, 16, 17, 18, 12)) +
    ggtitle(stringr::str_glue("{expnm} SymSim: {nind} Individuals / Condition")) +
    theme(plot.title = element_text(size = 20, hjust = 0.5))
  if (save_fig) {
    ggsave( plot = p,
           filename = file.path(save_dir,
                                stringr::str_glue("{expnm}_symsim_{nind}inds.pdf")),
           device = "pdf",
           width = 10,
           height = 7
           )
  }
  invisible(p)
}

r_mean_gather_all <- function(nind, nindeff, rows)  {
  r_mean_gather <- function(use_groupshift) {
    r <- load_symsim_result(nind = nind, nindeff = nindeff, use_groupshift = use_groupshift)
    mean_r <- f_across_rpt(r, f = mean)
    data <- as.data.frame(mean_r[rows, ])
    data$methods <- rownames(data)
    gr <- gather(data = data, key = "num_of_cell_per_ind", value = "AUC",
                 matches("[0-9][0-9]"))
    gr$num_of_cell_per_ind <- as.numeric(gr$num_of_cell_per_ind)
    gr$ninds <- nind
    gr$exp <- ifelse(use_groupshift, "GroupShift", "NonZeroShift")
    invisible(gr)
  }
  invisible(rbind(r_mean_gather(T), r_mean_gather(F)))
}



# main

## when nind equals to 3
nind <- 3
nindeff <- 1
### use groupshift
r <- load_symsim_result(nind = nind, nindeff = nindeff, use_groupshift = TRUE)
mean_r <- f_across_rpt(r, f = mean)
rows <- c("mssc_2-0_m", "mssc_2-1_m", "mssc_2-0_pm", "mssc_2-1_pm", "pseudo")
p <- scatter_plot(mean_r, rows = rows, nind = nind, use_groupshift = TRUE)

### do not use groupshift
r <- load_symsim_result(nind = nind, nindeff = nindeff, use_groupshift = FALSE)
mean_r <- f_across_rpt(r, f = mean)
rows <- c("mssc_2-0_m", "mssc_2-1_m", "mssc_2-0_pm", "mssc_2-1_pm", "pseudo")
p <- scatter_plot(mean_r, rows = rows, nind = nind, use_groupshift = FALSE)

## when nind equals to 5
nind <- 5
nindeff <- 2
### use groupshift
r <- load_symsim_result(nind = nind, nindeff = nindeff, use_groupshift = TRUE)
mean_r <- f_across_rpt(r, f = mean)
rows <- c("mssc_2-0_m", "mssc_2-1_m", "mssc_2-0_pm", "mssc_2-1_pm", "pseudo")
p <- scatter_plot(mean_r, rows = rows, nind = nind, use_groupshift = TRUE)

### do not use groupshift
r <- load_symsim_result(nind = nind, nindeff = nindeff, use_groupshift = FALSE)
mean_r <- f_across_rpt(r, f = mean)
rows <- c("mssc_2-0_m", "mssc_2-1_m", "mssc_2-0_pm", "mssc_2-1_pm", "pseudo")
p <- scatter_plot(mean_r, rows = rows, nind = nind, use_groupshift = FALSE)


## when nind equals to 10
nind  <- 10
nindeff <- 3

### use groupshift
r <- load_symsim_result(nind = nind, nindeff = nindeff, use_groupshift = TRUE)
mean_r <- f_across_rpt(r, f = mean)
rows <- c("mssc_2-0_m", "mssc_2-1_m", "mssc_2-0_pm", "mssc_2-1_pm", "pseudo")
p <- scatter_plot(mean_r, rows = rows, nind = nind, use_groupshift = TRUE)

### do not use groupshift
r <- load_symsim_result(nind = nind, nindeff = nindeff, use_groupshift = FALSE)
mean_r <- f_across_rpt(r, f = mean)
rows <- c("mssc_2-0_m", "mssc_2-1_m", "mssc_2-0_pm", "mssc_2-1_pm", "pseudo")
p <- scatter_plot(mean_r, rows = rows, nind = nind, use_groupshift = FALSE)


## merge different indiviudals

### use groupshift
use_groupshift <- TRUE
r_3 <- load_symsim_result(nind = 3, nindeff = 1, use_groupshift = use_groupshift)
r_5 <- load_symsim_result(nind = 5, nindeff = 2, use_groupshift = use_groupshift)
r_10 <- load_symsim_result(nind = 10, nindeff = 3, use_groupshift = use_groupshift)

mean_r3 <- f_across_rpt(r_3, f = mean)
mean_r5 <- f_across_rpt(r_5, f = mean)
mean_r10 <- f_across_rpt(r_10, f = mean)


rows <- c("mssc_2-0_m", "mssc_2-1_m", "mssc_2-0_pm", "mssc_2-1_pm", "pseudo")

data_r3 <- as.data.frame(mean_r3[rows, ])
data_r3$methods <- rownames(data_r3)
gr_3 <- gather(data = data_r3, key = "num_of_cell_per_ind", value = "AUC",
             matches("[0-9][0-9]"))
gr_3$num_of_cell_per_ind <- as.numeric(gr_3$num_of_cell_per_ind)
gr_3$ninds <- 3

data_r5 <- as.data.frame(mean_r5[rows, ])
data_r5$methods <- rownames(data_r5)
gr_5 <- gather(data = data_r5, key = "num_of_cell_per_ind", value = "AUC",
             matches("[0-9][0-9]"))
gr_5$num_of_cell_per_ind <- as.numeric(gr_5$num_of_cell_per_ind)
gr_5$ninds <- 5


data_r10 <- as.data.frame(mean_r10[rows, ])
data_r10$methods <- rownames(data_r10)
gr_10 <- gather(data = data_r10, key = "num_of_cell_per_ind", value = "AUC",
             matches("[0-9][0-9]"))
gr_10$num_of_cell_per_ind <- as.numeric(gr_10$num_of_cell_per_ind)
gr_10$ninds <- 10

gr_groupshift <- do.call("rbind", list(gr_3, gr_5, gr_10))
gr_groupshift$exp <- "GroupShift"
  
expnm <- ifelse(use_groupshift, "Group-shift", "Nonzero-shift")
p <- ggplot(gr_groupshift, aes(x = factor(num_of_cell_per_ind), y = AUC)) +
  geom_point(aes(colour = methods, shape = methods),
             size = 6, alpha = 0.8) +
  xlab("Number of cells / individual") +
  theme_bw(base_size = 15) + 
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 16, 17, 18, 12)) +
  ggtitle(stringr::str_glue("{expnm} SymSim: Individuals / Condition")) +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) + facet_grid(~ninds)



## use_nonzeroshift
use_groupshift <- FALSE
r_3 <- load_symsim_result(nind = 3, nindeff = 1, use_groupshift = use_groupshift)
r_5 <- load_symsim_result(nind = 5, nindeff = 2, use_groupshift = use_groupshift)
r_10 <- load_symsim_result(nind = 10, nindeff = 3, use_groupshift = use_groupshift)

mean_r3 <- f_across_rpt(r_3, f = mean)
mean_r5 <- f_across_rpt(r_5, f = mean)
mean_r10 <- f_across_rpt(r_10, f = mean)


rows <- c("mssc_2-0_m", "mssc_2-1_m", "mssc_2-0_pm", "mssc_2-1_pm", "pseudo")

data_r3 <- as.data.frame(mean_r3[rows, ])
data_r3$methods <- rownames(data_r3)
gr_3 <- gather(data = data_r3, key = "num_of_cell_per_ind", value = "AUC",
             matches("[0-9][0-9]"))
gr_3$num_of_cell_per_ind <- as.numeric(gr_3$num_of_cell_per_ind)
gr_3$ninds <- 3

data_r5 <- as.data.frame(mean_r5[rows, ])
data_r5$methods <- rownames(data_r5)
gr_5 <- gather(data = data_r5, key = "num_of_cell_per_ind", value = "AUC",
             matches("[0-9][0-9]"))
gr_5$num_of_cell_per_ind <- as.numeric(gr_5$num_of_cell_per_ind)
gr_5$ninds <- 5


data_r10 <- as.data.frame(mean_r10[rows, ])
data_r10$methods <- rownames(data_r10)
gr_10 <- gather(data = data_r10, key = "num_of_cell_per_ind", value = "AUC",
             matches("[0-9][0-9]"))
gr_10$num_of_cell_per_ind <- as.numeric(gr_10$num_of_cell_per_ind)
gr_10$ninds <- 10

gr_nonzero <- do.call("rbind", list(gr_3, gr_5, gr_10))
gr_nonzero$exp <- "NonZeroShift"

gr <- rbind(gr_groupshift, gr_nonzero)
p <- ggplot(gr, aes(x = factor(num_of_cell_per_ind), y = AUC)) +
  geom_point(aes(colour = methods, shape = methods),
             size = 6, alpha = 0.8) +
  xlab("Number of cells / individual") +
  theme_bw(base_size = 15) + 
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 16, 17, 18, 12)) +
  ggtitle(stringr::str_glue("SymSim: Individuals/Condition versus Experiments")) +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) + facet_grid(exp ~ ninds)

save_dir <- here::here("src", "symsim")
ggsave(filename = file.path(save_dir,
                            stringr::str_glue("SymSim_Summarize.pdf")),
       plot = p,
       device = "pdf",
       width = 12,
       height = 8
       )

## compare different measures
### mssc_20
rows <- paste("mssc_2-0", c("t", "bf", "m", "pbf", "pm"), sep = "_")
gr3 <- r_mean_gather_all(nind = 3, nindeff = 1, rows)
gr5 <- r_mean_gather_all(nind = 5, nindeff = 2, rows)
gr10 <- r_mean_gather_all(nind = 10, nindeff = 3, rows)
gr <- rbind(gr3, gr5, gr10)

p <- ggplot(gr, aes(x = factor(num_of_cell_per_ind), y = AUC)) +
  geom_point(aes(colour = methods, shape = methods),
             size = 5, alpha = 0.8) +
  xlab("Number of cells / individual") +
  theme_bw(base_size = 15) + 
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 16, 17, 10, 14)) +
  ggtitle(stringr::str_glue("SymSim: Individuals/Condition versus Experiments")) +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) + facet_grid(exp ~ ninds)

save_dir <- here::here("src", "symsim")
ggsave(filename = file.path(save_dir,
                            stringr::str_glue("SymSim_mssc20_Summarize.pdf")),
       plot = p,
       device = "pdf",
       width = 12,
       height = 8
       )


### mssc_21
rows <- paste("mssc_2-1", c("t", "bf", "m", "pbf", "pm"), sep = "_")
gr3 <- r_mean_gather_all(nind = 3, nindeff = 1, rows)
gr5 <- r_mean_gather_all(nind = 5, nindeff = 2, rows)
gr10 <- r_mean_gather_all(nind = 10, nindeff = 3, rows)
gr <- rbind(gr3, gr5, gr10)

p <- ggplot(gr, aes(x = factor(num_of_cell_per_ind), y = AUC)) +
  geom_point(aes(colour = methods, shape = methods),
             size = 5, alpha = 0.8) +
  xlab("Number of cells / individual") +
  theme_bw(base_size = 15) + 
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 16, 17, 10, 14)) +
  ggtitle(stringr::str_glue("SymSim: Individuals/Condition versus Experiments")) +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) + facet_grid(exp ~ ninds)

save_dir <- here::here("src", "symsim")
ggsave(filename = file.path(save_dir,
                            stringr::str_glue("SymSim_mssc21_Summarize.pdf")),
       plot = p,
       device = "pdf",
       width = 12,
       height = 8
       )

## mssc_20 vs mssc_21
rows <- c("mssc_2-0_m", "mssc_2-0_pm", "mssc_2-1_m", "mssc_2-1_pm")
gr3 <- r_mean_gather_all(nind = 3, nindeff = 1, rows)
gr5 <- r_mean_gather_all(nind = 5, nindeff = 2, rows)
gr10 <- r_mean_gather_all(nind = 10, nindeff = 3, rows)
gr <- rbind(gr3, gr5, gr10)

p <- ggplot(gr, aes(x = factor(num_of_cell_per_ind), y = AUC)) +
  geom_point(aes(colour = methods, shape = methods),
             size = 5, alpha = 0.8) +
  xlab("Number of cells / individual") +
  theme_bw(base_size = 15) + 
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(16, 17, 10, 14)) +
  ggtitle(stringr::str_glue("SymSim: Individuals/Condition versus Experiments")) +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) + facet_grid(exp ~ ninds)

save_dir <- here::here("src", "symsim")
ggsave(filename = file.path(save_dir,
                            stringr::str_glue("SymSim_mssc-20vs21_Summarize.pdf")),
       plot = p,
       device = "pdf",
       width = 12,
       height = 8
       )

