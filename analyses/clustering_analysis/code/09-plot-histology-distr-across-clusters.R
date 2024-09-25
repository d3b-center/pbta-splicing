################################################################################
# 09-plot-histology-distr-across-clusters.R
# Plot the distribution of histologies across clusters
#
# Author: Ammar Naqvi, Jo Lynne Rokita
################################################################################

## libraries 
suppressPackageStartupMessages({
  library("tidyverse")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`


## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "clustering_analysis")

plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")
output_dir <- file.path(analysis_dir, "output")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## filepaths 
optimal_cluster_tsv <- file.path(output_dir, "optimal_clustering", "ccp_optimal_clusters.tsv")
cluster_membership_tsv <-  file.path(output_dir, "cluster_members_by_cancer_group_subtype.tsv")
subtype_hex <- file.path(input_dir, "subtype_hex.tsv")

# read in palette and cluster file
cluster_df <- read_tsv(optimal_cluster_tsv) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = sample)

histologies_df <- read_tsv(file.path(root_dir,"data", "histologies-plot-group.tsv"), guess_max = 100000) %>%
  dplyr::select(Kids_First_Biospecimen_ID, broad_histology, cancer_group, plot_group, plot_group_hex, molecular_subtype) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% cluster_df$Kids_First_Biospecimen_ID) %>%
  unique() %>%
  mutate(molecular_subtype_display = case_when(grepl("KIAA", molecular_subtype) ~ "KIAA1549--BRAF",
                                               grepl("LGG, BRAF V600|LGG, RTK, BRAF V600E", molecular_subtype) ~ "BRAF V600E",
                                               grepl("LGG, IDH|LGG, other MAPK, IDH", molecular_subtype) ~ "IDH",
                                               molecular_subtype %in% c("LGG, wildtype", "SEGA, wildtype") ~ "Wildtype",
                                               grepl("To be classified", molecular_subtype) ~ "To be classified",
                                               grepl("LGG|SEGA, RTK", molecular_subtype) ~ "Other MAPK",
                                               grepl("HGG, IDH", molecular_subtype) ~ "IDH",
                                               grepl("HGG, H3 wildtype", molecular_subtype) ~ "H3 wildtype",
                                               grepl("HGG, PXA", molecular_subtype) ~ "PXA",
                                               grepl("IHG", molecular_subtype) ~ "IHG",
                                               cancer_group == "Diffuse intrinsic pontine glioma" ~ "DIPG",
                                               grepl("H3 G35", molecular_subtype) ~ "H3 G35",
                                               grepl("H3 K28", molecular_subtype) ~ "H3 K28",
                                               is.na(molecular_subtype) & plot_group == "Other high-grade glioma" ~ "Oligodendroglioma",
                                               
         TRUE ~ molecular_subtype)) %>%
  dplyr::right_join(cluster_df)


color_df <- histologies_df %>%
  dplyr::select(plot_group_hex, plot_group) %>%
  dplyr::filter(!is.na(plot_group)) %>%
  unique()
cols <- as.character(color_df$plot_group_hex)
names(cols) <- as.character(color_df$plot_group)

# create plot
pdf(file.path(plots_dir, "cluster_membership.pdf"), height = 4, width = 7)
ggplot(histologies_df, aes(fill=plot_group, x= factor(cluster_assigned))) +
  geom_bar(stat="count", position="stack") + 
  xlab("Cluster") + ylab("Frequency") +
  scale_fill_manual("Histology", values = cols) + 
  theme_Publication()
dev.off()

## stratify by molecular subtype
# subset to only those plot groups: ATRT, MB, LGG, HGG, EPN
hist_subset <- histologies_df %>%
  filter(plot_group %in% c("Medulloblastoma", 
                           "Ependymoma",
                           "Low-grade glioma",
                           "Other high-grade glioma",
                           "DIPG or DMG")) %>%
  mutate(plot_group = case_when(plot_group %in% c("Other high-grade glioma",
                                                  "DIPG or DMG") ~ "High-grade glioma",
                                TRUE ~ plot_group))

subtype_hex_codes <- read_tsv(subtype_hex)
cols <- subtype_hex_codes$hex_code
names(cols) <- subtype_hex_codes$molecular_subtype_display

pdf(file.path(plots_dir, "cluster_membership-subtypes.pdf"), height = 6, width = 11)
hist_subset %>% 
  ggplot(aes(fill=molecular_subtype_display, x= factor(cluster_assigned))) +
  geom_bar(stat="count", position="stack", color = "black", size = 0.2) + 
  facet_wrap(~plot_group, nrow = 2, scales = "free_y") +
  xlab("Cluster") + 
  ylab("Frequency") +
  theme_Publication() +
  scale_fill_manual(name = "Molecular Subtype", values = cols) 
dev.off()


# write out bs id with cancer group, subtype, and cluster
histologies_df %>%
  dplyr::select(Kids_First_Biospecimen_ID, cancer_group, plot_group, molecular_subtype, molecular_subtype_display, cluster_assigned) %>%
  write_tsv(cluster_membership_tsv)
