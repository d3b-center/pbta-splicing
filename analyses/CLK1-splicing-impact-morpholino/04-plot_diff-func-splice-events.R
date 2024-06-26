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

## combine all splice types together
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
  count(Type, Uniprot_wrapped, Preference) %>%
  # add n = for first n
  mutate(n = as.character(n),
         n = case_when(Type == "A3SS" & Uniprot_wrapped == "Disulfide\nBond" & Preference == "Inclusion" ~ paste0("n = ",n),
                       TRUE ~ as.character(n)))

# Your ggplot code with adjustments
plot_dsp <- ggplot(psi_comb, aes(x = Uniprot_wrapped, y = dPSI * 100, fill = Preference)) +
  geom_boxplot(aes(color = Preference), position = position_dodge(width = 0.9), outlier.shape = NA, size = 0.5, alpha = 0.4) +
  ggforce::geom_sina(aes(color = Preference), pch = 16, size = 3, position = position_dodge(width = 0.9), alpha = 0.4) +
  facet_wrap(~Type, nrow = 1) +
  scale_color_manual(name = "Preference", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A")) +
  scale_fill_manual(name = "Preference", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A")) +
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

## subset by GOI
# gene list files
known_rbp_file <- file.path(input_dir,'RBP_known.txt')
known_epi_file <- file.path(input_dir,'epi_known.txt')

## subset with gene lists
genelist_ref_df <-read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuseData"))  %>%
  # remove cosmic census
  filter(type != "CosmicCensus") %>%
  mutate(type = gsub("CosmicCensus, |, CosmicCensus", "", type))

## other non-annoFuse gene lists, RBPs and epigenetic genes
known_rbp_inlist <- read.table(known_rbp_file,header=FALSE) %>% 
  dplyr::rename('Gene_Symbol' = V1) %>%
  filter(Gene_Symbol %in% genelist_ref_df$Gene_Symbol)

known_rbp_not_inlist <- read.table(known_rbp_file,header=FALSE) %>% 
  dplyr::rename('Gene_Symbol' = V1) %>%
  filter(!Gene_Symbol %in% genelist_ref_df$Gene_Symbol) %>%
  mutate(plot_type = "Other",
         plot_subtype = "RNA Binding Protein")

known_epi_inlist <- read.table(known_epi_file,header=FALSE) %>% 
  dplyr::rename('Gene_Symbol' = V1) %>%
  filter(Gene_Symbol %in% genelist_ref_df$Gene_Symbol)

known_epi_not_inlist <- read.table(known_epi_file,header=FALSE) %>% 
  dplyr::rename('Gene_Symbol' = V1) %>%
  filter(!Gene_Symbol %in% genelist_ref_df$Gene_Symbol) %>%
  mutate(plot_type = "Other",
         plot_subtype = "Epigenetic")


genelist_cat <- genelist_ref_df %>%
  mutate(type = ifelse(Gene_Symbol %in% known_rbp_inlist$Gene_Symbol, paste(type, "RNA Binding Protein", sep = ", "), type),
         type = ifelse(Gene_Symbol %in% known_epi_inlist$Gene_Symbol, paste(type, "Epigenetic", sep = ", "), type)) %>%
  # add epi/rbp for those which are already in list and recategorize
  mutate(plot_type = case_when(grepl("TumorSuppressorGene, Oncogene", type) ~ "Oncogene or Tumor Suppressor",
                               type %in% c("TumorSuppressorGene, TranscriptionFactor", "TumorSuppressorGene, Kinase", "TumorSuppressorGene", 
                                           "Kinase, TumorSuppressorGene", "TumorSuppressorGene, RNA Binding Protein",
                                           "TumorSuppressorGene, TranscriptionFactor, RNA Binding Protein, Epigenetic",
                                           "TumorSuppressorGene, RNA Binding Protein, Epigenetic", 
                                           "TumorSuppressorGene, TranscriptionFactor, Epigenetic", "TumorSuppressorGene, Epigenetic",
                                           "TumorSuppressorGene, Kinase, Epigenetic", "Kinase, TumorSuppressorGene, Epigenetic") ~ "Tumor Suppressor",
                               type %in% c("Kinase, Oncogene", "Oncogene", "Oncogene, Kinase", "Oncogene, TranscriptionFactor",
                                           "Oncogene, RNA Binding Protein", "Oncogene, TranscriptionFactor, RNA Binding Protein",
                                           "Oncogene, TranscriptionFactor, Epigenetic", "Oncogene, Kinase, Epigenetic",
                                           "Oncogene, Epigenetic", "Kinase, Oncogene, Epigenetic", "Oncogene, RNA Binding Protein, Epigenetic",
                                           "TranscriptionFactor, Oncogene") ~ "Oncogene",
                               TRUE ~ "Other"),
         plot_subtype = case_when(type %in% c("Kinase", "Kinase, TranscriptionFactor", "Kinase, Oncogene", "TumorSuppressorGene, Kinase", 
                                              "Kinase, TumorSuppressorGene, Oncogene", "Oncogene, Kinase", "Kinase, RNA Binding Protein",
                                              "TumorSuppressorGene, Oncogene, Kinase, Epigenetic", "TumorSuppressorGene, Kinase, Epigenetic",
                                              "Oncogene, Kinase, Epigenetic", "Kinase, TumorSuppressorGene, Oncogene, Epigenetic",
                                              "Kinase, TumorSuppressorGene, Epigenetic", "Kinase, Oncogene, Epigenetic", "Kinase, Epigenetic",
                                              "Kinase, TumorSuppressorGene", "TumorSuppressorGene, Oncogene, Kinase") ~ "Kinase",
                                  type %in% c("Oncogene, TranscriptionFactor", "TumorSuppressorGene, Oncogene, TranscriptionFactor",
                                              "TumorSuppressorGene, TranscriptionFactor", "TranscriptionFactor, Oncogene") ~ "Transcription Factor",
                                  type %in% c("Oncogene, RNA Binding Protein", "TranscriptionFactor, RNA Binding Protein", 
                                              "TumorSuppressorGene, Oncogene, RNA Binding Protein", "TranscriptionFactor, RNA Binding Protein, Epigenetic",
                                              "TumorSuppressorGene, RNA Binding Protein, Epigenetic", "TumorSuppressorGene, TranscriptionFactor, RNA Binding Protein",
                                              "TumorSuppressorGene, TranscriptionFactor, RNA Binding Protein, Epigenetic",
                                              "Oncogene, RNA Binding Protein, Epigenetic", "Oncogene, TranscriptionFactor, RNA Binding Protein",
                                              "TumorSuppressorGene, RNA Binding Protein") ~ "RNA Binding Protein",
                                  type %in% c("TumorSuppressorGene, TranscriptionFactor, Epigenetic", "TumorSuppressorGene, Oncogene, TranscriptionFactor, Epigenetic",
                                              "TumorSuppressorGene, Oncogene, Epigenetic", "TumorSuppressorGene, Epigenetic",
                                              "TranscriptionFactor, Epigenetic", "Oncogene, TranscriptionFactor, Epigenetic",
                                              "Oncogene, Epigenetic") ~ "Epigenetic",
                                  type == "TumorSuppressorGene" ~ "Other Tumor Suppressor",
                                  type == "Oncogene" ~ "Other Oncogene",
                                  type == "TumorSuppressorGene, Oncogene" ~ "Oncogene or Tumor Suppressor",
                                  type == "TranscriptionFactor" ~ "Transcription Factor",
                                  TRUE ~ type)) %>%
  dplyr::select(Gene_Symbol, plot_type, plot_subtype) %>%
  # add other RBP, Epi not already in the list
  bind_rows(known_rbp_not_inlist, known_epi_not_inlist) %>% 
  dplyr::rename('gene'=Gene_Symbol) %>%
  write_tsv(file.path(results_dir, "gene_categories.tsv"))

# check that all of the subgroups look right
table(genelist_cat$plot_subtype, genelist_cat$plot_type)

psi_comb_goi <- psi_comb %>% inner_join(genelist_cat, by="gene")   

# relevel the plot_type and subtype
psi_comb_goi$plot_type <- factor(psi_comb_goi$plot_type, levels = c("Oncogene", "Tumor Suppressor",
                                                                    "Oncogene or Tumor Suppressor", "Other"))
psi_comb_goi$plot_subtype <- factor(psi_comb_goi$plot_subtype, levels = c("Kinase", "Oncogene or Tumor Suppressor",
                                                                          "Transcription Factor", "RNA Binding Protein",
                                                                          "Epigenetic", "Other Oncogene", "Other Tumor Suppressor"))
unique(psi_comb_goi$plot_subtype)

psi_comb_goi_tab <- psi_comb_goi %>% dplyr::select(SpliceID,gene, dPSI,Type, Uniprot, Preference,plot_type, plot_subtype)

# write for supplemental 
write_tsv(psi_comb_goi_tab, file.path(results_dir, "differential_splice_by_goi_category.tsv"))

## plot num of hits per gene fam
plot_barplot_family <- ggplot(psi_comb_goi, aes(x = fct_rev(fct_infreq(plot_subtype)), fill= Preference)) +
  geom_bar(stat="count", position='dodge', color="black") + 
  facet_wrap(~plot_type, scales = "free_y", ncol = 1) +
  xlab("Gene Family")     + 
  ylab("Number of Genes Signficantly Mis-spliced") + 
  scale_fill_manual(name = "Preference",
                    values=c("#FFC20A","#0C7BDC")) + 
  geom_text(stat='count',aes(label=after_stat(count)), 
            position = position_dodge(width = 1),
            hjust = -0.5, size = 3.5) +
  theme_Publication() +
  coord_flip() +
  ylim(0,150)

# Save plot as PDF
pdf(file_dpsi_goi_plot, height = 8, width = 8, useDingbats = FALSE) 
plot_barplot_family
dev.off()


