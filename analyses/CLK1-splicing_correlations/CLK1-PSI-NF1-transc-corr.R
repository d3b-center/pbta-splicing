################################################################################
# CLK1-PSI-NF1-transc-corr.R
# Plot scatter and compute R for CLK1 PSI and NF1 transcripts in HGGs
# written by Ammar S Naqvi
# Usage: Rscript CLK1-PSI-NF1-transc-corr.R 
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("ggplot2")
  library("tidyverse")
  library("ggpubr")
  library("vroom")
  library('data.table')
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
results_dir <- file.path(analysis_dir, "results")

plots_dir <- file.path(analysis_dir, "plots")
figures_dir <- file.path(root_dir, "figures")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
source(file.path(figures_dir, "theme_for_plots.R"))
source(file.path(analysis_dir, "util", "function-create-scatter-plot.R"))

## define input files
clin_file <- file.path(data_dir,"histologies-plot-group.tsv")
rsem_transc_counts <- file.path(data_dir,"rna-isoform-expression-rsem-tpm.rds")
rmats_clk1_file <- file.path(data_dir, "clk1-splice-events-rmats.tsv")
rmats_nf1_file <- file.path(results_dir, "nf1-splice-events-rmats.tsv")

## get CLK1 psi values 
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_df <- vroom(indep_file)

## read in histology file and count data
## filter histology file for all HGG, only stranded samples
histologies_df  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         experimental_strategy == "RNA-Seq",
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID,
         RNA_library == 'stranded'
  )

hgg_bs_id <- histologies_df %>%
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) 


## all transcripts
rsem_all_transc_df <- readRDS(rsem_transc_counts) %>%
  filter(grepl("^NF1-", gene_symbol))

rsem_all_transc_tv_df <- rsem_all_transc_df %>% gather(key = "sample_id", value = "TPM", -transcript_id, -gene_symbol) %>%
  arrange(transcript_id,sample_id)

rmats_exp_CLK1_NF1_df <- rsem_all_transc_tv_df %>%
  inner_join(clk1_HGG_rmats, by= "sample_id") %>%
  # need to add this so that we can use this in the plot axis
  #mutate(transcript = paste0(transcript_id)) %>%
  dplyr::filter(TPM > 1) %>% 
  mutate(logExp = log(TPM, 2))


# Convert transcript_id to factor for easier plotting
rmats_exp_CLK1_NF1_df$transcript_id <- factor(rmats_exp_CLK1_NF1_df$transcript_id)

p <- ggplot(rmats_exp_CLK1_NF1_df, aes(x = IncLevel1, y = logExp)) +
  geom_point(colour = "black") +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              colour = "red",
              fill = "pink",
              linetype="dashed") +
  labs(x = "Exon 4 CLK1 PSI",y = "NF1 TPM (log2)") + 
  stat_cor(method = "pearson",
           label.x = 0, label.y =0, size = 2) +
  xlim(c(0,1.1)) +
  #ylim(c(0,4)) +
  facet_wrap(~transcript_id, nrow = 7, scales = "free_y") + 
  theme_Publication()

# save plot
pdf(file.path(paste(plots_dir, "/NF1-transc-expr_vs_CLK1Ex4-PSI_hgg.pdf", sep = "")), width = 14, height = 10)
print(p)
dev.off()
