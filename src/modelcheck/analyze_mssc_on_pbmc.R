## Analyze the mssc vi results on PBMC dataset.

suppressPackageStartupMessages(library(tidyverse))
import::from(here, here)
library(GGally)
library(hrbrthemes)
library(viridis)

## local modules
options("import.path" = here("rutils"))
myt <- modules::import("transform")
myfit <- modules::import("myfitdistr")
mypseudo <- modules::import("pseudobulk")
mypbmc <- modules::import("pbmc")

## warnings/errors traceback settings
options(error = traceback)
options(warn = 0)


## * config
cell_type <- "Naive CD4+ T"
num_top_gene <- 1000L
cmdstan_version <- "2.23.0"
mssc_stan_fnm <- "v3-1.stan"
output_stan_fnm <- "pbmc_naivecd4t_top1000_stan_v3-1_vi.rds"



## ** for plots
myggtitle <- theme(plot.title = element_text(size = 15,
  hjust = 0.5, family = "Helvetica", color = "black",
  face = "bold"),
axis.title = element_text(size = 13, family = "Helvetica", color = "black"),
axis.text = element_text(size = 13, family = "Helvetica", color = "black"),
legend.text = element_text(size = 13, family = "Helvetica", color = "black"),
legend.title = element_text(size = 13, family = "Helvetica", color = "black"))


## * load cmdstan
library(cmdstanr)
set_cmdstan_path(path = paste(Sys.getenv("HOME"), "softwares",
  paste0("cmdstan-", cmdstan_version), sep = "/"))

library(bayesplot)
library(posterior)

## * load mssc result
mssc_pbmc_res <- readRDS(here("exps", "pbmc",
  "vi", output_stan_fnm))

## * result analysis
## ** mssc rank genes
mucondt <- mssc_pbmc_res$mucondt
mssc_ranked_index <- order(abs(mucondt), descreasing = T)
mssc_ranked_genes <- mucondt[order(abs(mucondt), decreasing = T)]

pseudo_ranked_genes <- mssc_pbmc_res$pseudwilcoxp

## *** generate parallel figure
pseudo_scores <- seq(from = length(pseudo_ranked_genes), to = 1, by = -1)
names(pseudo_scores) <- names(pseudo_ranked_genes)

mssc_scores <- seq(from = length(mssc_ranked_genes), to = 1, by = -1)
names(mssc_scores) <- names(mssc_ranked_genes)


get_comp_figure <- function(ingroupnum = 50, genenum = 150,
                            printnm = T, hjust = 1.1, vjust = 1) {
  genes <- names(pseudo_scores)[seq_len(genenum)]
  rank_comp <- data.frame(pseudo = pseudo_scores[genes],
    mssc = mssc_scores[genes],
    group = c(rep("Top_pseudo", ingroupnum),
      rep("Middle_pseudo", genenum - 2 * ingroupnum),
      rep("Last_pseudo", ingroupnum)),
    genenm = genes
  )
  rownames(rank_comp) <- genes

  p <- ggparcoord(rank_comp,
    columns = 1:2,
    showPoints = TRUE,
    groupColumn = 3,
    alphaLines = 0.1,
    title = "Naive CD4+ T cell on PBMC",
  )
  if (printnm) {
    p <- p +
      geom_text(data = rank_comp %>%
        select(genenm) %>%
        mutate(x = 1,
          y = pseudo),
      aes(x = x, y = scale(pseudo), label = genenm),
      hjust = hjust,
      vjust = vjust,
      inherit.aes = FALSE)
  }

  p <- p +
    scale_color_manual(values = c("Top_pseudo" = "red",
      "Last_pseudo" = "blue",
      "Middle_pseudo" = "black")) +
    xlab("Method") +
    ylab("Inverse Rank Index") +
    myggtitle
  invisible(p)
}

p <-  get_comp_figure()
ggsave(path = here("exps", "pbmc", "vi"),
  filename = "pbmc_naivecd4t_top150.pdf",
  device = "pdf",
  plot = p,
  width = 7,
  height = 30)

