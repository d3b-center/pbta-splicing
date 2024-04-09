# 05-plot-tmb-sbi.R
# written by Ammar Naqvi & Jo Lynne Rokita
#
# This script plots TMB based on high vs low splicing burden 
# usage: Rscript 05-plot-tmb-sbi.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("tidyverse")
  library("ggplot2")
  library("ggpubr")
} )

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
map_dir <- file.path(root_dir, "analyses", "cohort_summary", "results")
input_dir <- file.path(analysis_dir, "input")
data_dir <- file.path(root_dir, "data")
tmb_dir <- file.path(root_dir, "analyses", "oncoprint", "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

## Call plot publication theme script 
source(file.path(root_dir, "figures", "theme_for_plots.R"))

## define output files
boxplot_sbi_vs_tmb_by_mutation_status_file  <- file.path(plots_dir, "boxplot_sbi-tmb-by-mutation-status.pdf")
boxplot_sbi_vs_tmb_by_cg_file  <- file.path(plots_dir, "boxplot_sbi-tmb-by-cg.pdf")
corplot_sbi_vs_tmb_file <- file.path(plots_dir, "corplot_sbi-tmb.pdf")
corplot_sbi_vs_tmb_by_cg_file <- file.path(plots_dir, "corplot_sbi-tmb-by-cg.pdf")

## input files
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_wgs_file <- file.path(data_dir, "independent-specimens.wgswxspanel.primary.prefer.wgs.tsv")
tmb_coding_file  <- file.path(tmb_dir,"snv-mutation-tmb-coding.tsv") # OPC v13 TMB file
sbi_coding_file  <- file.path(results_dir,"splicing_index.SE.txt")
palette_file <- file.path(map_dir,"histologies-plot-group.tsv") 

# read in files
hist_pal <- read_tsv(palette_file) %>%
  filter(!is.na(pathology_diagnosis),
         !is.na(plot_group)) %>%
  dplyr::select(Kids_First_Biospecimen_ID, match_id, cancer_group, plot_group)

indep_rna_df <- read_tsv(indep_rna_file) %>% 
  filter(cohort == "PBTA") %>%
 # dplyr::rename(Kids_First_Biospecimen_ID_RNA = Kids_First_Biospecimen_ID) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID)
indep_wgs_df <- read_tsv(indep_wgs_file) %>% 
  filter(cohort == "PBTA") %>%
  #dplyr::rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID)

## get tmb file (source: OpenPedCan v13)
tmb_coding_df  <-  read_tsv(tmb_coding_file)  %>% 
  mutate(mutation_status = case_when(tmb <10 ~ "Normal",
                                     tmb >=10 & tmb < 100 ~ "Hypermutant",
                                     tmb >=100 ~ "Ultra-hypermutant")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  left_join(hist_pal[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  dplyr::select(-experimental_strategy) %>%
  right_join(indep_wgs_df) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID)


sbi_coding_df  <-  read_tsv(sbi_coding_file) %>% 
  dplyr::rename(Kids_First_Biospecimen_ID = Sample) %>%
  left_join(hist_pal[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  dplyr::select(-Histology) %>%
  right_join(indep_rna_df) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_RNA = Kids_First_Biospecimen_ID)


## intersect tmb values with SBI tumors 
sbi_vs_tmb_innerjoin_df <- tmb_coding_df %>%
  inner_join(sbi_coding_df, by=c("Kids_First_Participant_ID", "match_id")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_RNA) %>%
  # add histology
  left_join(hist_pal) 

sbi_tmb_no_hyper <- sbi_vs_tmb_innerjoin_df %>%
  filter(mutation_status == "Normal")

# generate corplot
## create plot
df_list <- list(sbi_vs_tmb_innerjoin_df, sbi_tmb_no_hyper)
  
pdf(file.path(corplot_sbi_vs_tmb_file), width = 4.5, height = 4.5)

for (each_df in df_list) {
  p <- ggscatter(each_df, 
                         y= "SI", 
                         x= "tmb", 
                         add = "reg.line", 
                         conf.int = TRUE, 
                         cor.coef = TRUE, 
                         cor.method = "pearson",
                         add.params = list(color = "red",
                                           fill = NA),
                         ticks = TRUE,
                 alpha = 0.6) + 
  ylab("Splicing Burden Index") +
  xlab("Tumor Mutation Burden") +
  theme_Publication()
  print(p)
}

dev.off()


## identify samples by high vs low SBI tumors
quartiles_sbi <- quantile(sbi_vs_tmb_innerjoin_df$SI, probs=c(.25, .75), na.rm = TRUE)
lower_sbi <- quartiles_sbi[1]
upper_sbi <- quartiles_sbi[2]

high_vs_low_df <- sbi_vs_tmb_innerjoin_df %>%
  mutate(SBI_level = case_when(SI > upper_sbi ~ "High",
                               SI < lower_sbi ~ "Low")) %>%
  # remove samples without either TMB or SI data
  filter(!is.na(SBI_level),
         !is.na(mutation_status))

# relevel for plotting
high_vs_low_df$mutation_status <- factor(high_vs_low_df$mutation_status, levels = c("Normal", "Hypermutant", "Ultra-hypermutant"))
high_vs_low_df$SBI_level <- factor(high_vs_low_df$SBI_level, levels = c("Low", "High"))

sbi_tmb_plot <- ggplot(high_vs_low_df, aes(SBI_level, log10(tmb))) +  
  ggforce::geom_sina(aes(color = SBI_level, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  stat_compare_means(label.y = 2.9) + 
  facet_wrap("mutation_status") +
  scale_color_manual(name = "SBI_level", values = c(High = "#0C7BDC", Low = "#FFC20A")) + 
  theme_Publication() + 
  labs(y="log10 (TMB)", x="Splicing Burden Level") + 
  ylim(c(-3,3)) +
  theme(legend.position="none")

## save pdf
pdf(boxplot_sbi_vs_tmb_by_mutation_status_file, width = 7, height = 4)
sbi_tmb_plot
dev.off()


## remove hyper/ultra-mutant samples to assess histology

# Calculate upper and lower quantiles
quantile_data <- sbi_vs_tmb_innerjoin_df %>%
  # for normals only
  dplyr::filter(mutation_status == "Normal") %>%
  # remove samples without either TMB or SI data
  filter(!is.na(SI),
         !is.na(mutation_status),
         !is.na(plot_group)) %>%
  group_by(plot_group) %>%
  summarize(lower_quantile = quantile(SI, 0.25),
            upper_quantile = quantile(SI, 0.75)) %>%
  inner_join(sbi_vs_tmb_innerjoin_df, by = "plot_group")

by_hist <- quantile_data %>%
  mutate(SBI_level = case_when(SI > upper_quantile ~ "High",
                               SI < lower_quantile ~ "Low",
                               TRUE ~ "Mid")) %>%
  filter(SBI_level != "Mid") %>%
  filter(plot_group != "Other tumor",
         plot_group != "Non−neoplastic tumor")

# relevel for plotting
by_hist$SBI_level <- factor(by_hist$SBI_level, levels = c("Low", "High"))
# calculate n per group and retain those n >= 3
counts <- by_hist %>%
  group_by(plot_group) %>%
  count(SBI_level) %>%
  filter(n >=3)
  
by_hist <- by_hist %>%
  filter(plot_group %in% counts$plot_group)
  
# plot
pdf(boxplot_sbi_vs_tmb_by_cg_file, width = 12.5, height = 5.5)
ggplot(by_hist, aes(SBI_level, log10(tmb))) +  
  ggforce::geom_sina(aes(color = SBI_level, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  stat_compare_means(label.y = 2,size = 3) + 
  facet_wrap("plot_group", labeller = labeller(plot_group = label_wrap_gen(width = 18)), nrow  = 2) +
  scale_color_manual(name = "SBI_level", values = c(High = "#0C7BDC", Low = "#FFC20A")) + 
  theme_Publication() + 
  labs(y="log10 (TMB)", x="Splicing Burden Level") + 
  ylim(c(-2,2.5)) +
  theme(legend.position = "none", strip.text = element_text(size = 10))  # Adjust the size here
dev.off()

# plot cancer group corplots
sbi_tmb_no_hyper_subset <- sbi_tmb_no_hyper %>%
  filter(plot_group != "Other tumor",
         plot_group != "Non−neoplastic tumor")

pdf(corplot_sbi_vs_tmb_by_cg_file, width = 16, height = 8)
ggscatter(sbi_tmb_no_hyper, 
          y= "SI", 
          x= "tmb", 
          add = "reg.line", 
          conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.method = "pearson",
          add.params = list(color = "red",
                            fill = NA),
          ticks = TRUE,
          alpha = 0.6) + 
  facet_wrap("plot_group", labeller = labeller(plot_group = label_wrap_gen(width = 18)), 
             nrow  = 3,
             scales = "free") +
  theme_Publication() + 
  labs(x="Tumor Mutation Burden", y="Splicing Burden Index") + 
  theme(legend.position = "none", strip.text = element_text(size = 13))  # Adjust the size here
dev.off()




