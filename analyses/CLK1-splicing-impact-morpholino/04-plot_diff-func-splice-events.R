################################################################################
# 05-plot_diff-splice-events.R
# written byAmmar Naqvi and Jo Lynne Rokita
#
# usage: Rscript 03-plot_diff-splice-events.R
################################################################################

suppressPackageStartupMessages({
  library("ggplot2")
  library("tidyverse")
  library("ggrepel")
  library("vroom")
  library("ggpubr")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories and file paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")
input_dir <- file.path(analysis_dir, "input")

results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}


##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## define output files
file_dpsi_plot <- file.path(plots_dir,"dPSI-distr-func.pdf")
file_dpsi_goi_plot <- file.path(plots_dir,"dPSI-distr-func-goi.pdf")

## get and setup input
## retrieve psi values from tables
file_psi_SE_func <- file.path(results_dir,"splicing_events.morpho.SE.intersectUnip.ggplot.txt")
file_psi_RI_func <- file.path(results_dir,"splicing_events.morpho.RI.intersectUnip.ggplot.txt")
file_psi_A5SS_func <- file.path(results_dir,"splicing_events.morpho.A5SS.intersectUnip.ggplot.txt")
file_psi_A3SS_func <- file.path(results_dir,"splicing_events.morpho.A3SS.intersectUnip.ggplot.txt")

## gene categories file 
oncokb_gene_ref <- read_tsv(file.path(results_dir, "gene_categories.tsv"))

## combine all splice types together, everything relative to CLK1 exon 4 cells (untreated)
dpsi_unip_incl <- vroom(c(file_psi_SE_func, file_psi_RI_func, file_psi_A5SS_func, file_psi_A3SS_func)) %>%
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>%
  filter(dPSI<0) %>% 
  mutate(Preference='Inclusion',
         dPSI=abs(dPSI)) 

dpsi_unip_skp <- vroom(c(file_psi_SE_func, file_psi_RI_func, file_psi_A5SS_func, file_psi_A3SS_func)) %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>%
  filter(dPSI>0) %>% 
  mutate(Preference='Skipping')

psi_comb <- rbind(dpsi_unip_incl,dpsi_unip_skp) %>% 
  mutate(Uniprot = case_when(Uniprot == 'DisulfBond' ~ "Disulfide Bond",
                             Uniprot == 'LocSignal' ~ "Localization Signal",
                             Uniprot == 'Mod' ~ 'Modification',
                             .default = Uniprot),
         Uniprot_wrapped = stringr::str_wrap(Uniprot, width = 10)
  )

## ggstatplot across functional sites
set.seed(123)
counts_psi_comb <- psi_comb %>% 
  dplyr::count(Type, Uniprot_wrapped, Preference) %>%
  # add n = for first n
  mutate(n = as.character(n),
         n = case_when(Type == "A3SS" & Uniprot_wrapped == "Disulfide\nBond" & Preference == "Inclusion" ~ paste0("n = ",n),
                       TRUE ~ as.character(n)))

# Your ggplot code with adjustments
plot_dsp <- ggplot(psi_comb, aes(x = Uniprot_wrapped, y = dPSI * 100, fill = Preference)) +
  geom_boxplot(aes(color = Preference), position = position_dodge(width = 0.9), outlier.shape = NA, size = 0.5, alpha = 0.4) +
  ggforce::geom_sina(aes(color = Preference), pch = 16, size = 3, position = position_dodge(width = 0.9), alpha = 0.4) +
  facet_wrap(~Type, nrow = 1) +
  scale_color_manual(name = "Preference (CLK1 exon 4 high)", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A")) +
  scale_fill_manual(name = "Preference (CLK1 exon 4 high)", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A")) +
  theme_Publication() +
  labs(y = "Percent Spliced In (PSI)", x = "Uniprot-defined Functional Site") +
  geom_text(data = counts_psi_comb, aes(label = paste(n), x = Uniprot_wrapped, y = 10), vjust = 3, size = 3, position = position_dodge(width = 0.9)) +
  theme(
    legend.position = "top",  # Move legend to the top
    legend.direction = "horizontal",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylim(c(-15, 100))

# Save plot as PDF
pdf(file_dpsi_plot, 
    width = 12, height = 4)
print (plot_dsp)
dev.off()


psi_comb_goi <- psi_comb %>% 
  inner_join(oncokb_gene_ref, by="gene") %>% 
  unique() %>%
# write for supplemental 
write_tsv(file.path(results_dir, "differential_splice_by_goi_category.tsv"))

psi_comb_goi_subset_for_plot <- psi_comb_goi %>%
  select(gene, Preference, classification) %>%
  unique()

# relevel to plot
psi_comb_goi_subset_for_plot$classification <- factor(psi_comb_goi_subset_for_plot$classification, levels = c("Oncogene", "Tumor Suppressor", "Both"))

## plot num of hits per gene fam
plot_barplot_family <- ggplot(psi_comb_goi_subset_for_plot, aes(x = classification, fill= Preference)) +
  geom_bar(stat="count", position='dodge', color="black") + 
  facet_wrap(~classification, scales = "free_y", ncol = 1) +
  xlab("Cancer Gene Type")     + 
  ylab("Number of Genes Signficantly Mis-spliced") + 
  scale_fill_manual(name = "Preference\n(CLK1 exon 4 high)",
                    values=c("#FFC20A","#0C7BDC")) + 
  geom_text(stat='count',aes(label=after_stat(count)), 
            position = position_dodge(width = 1),
            hjust = -0.5, size = 3.5) +
  theme_Publication() +
  theme(legend.position = "top", legend.direction = "horizontal") +
  ylim(c(0,85))+
  coord_flip() 

# Save plot as PDF
pdf(file_dpsi_goi_plot, height = 5.5, width = 6.6, useDingbats = FALSE) 
plot_barplot_family
dev.off()


