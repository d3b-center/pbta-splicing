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


# 
fusion_info_df <- vroom(file.path(input_dir,"lgg-braf-fusion-breakpoint-annotation.tsv"))
histologies_w_fusion_df <- histologies_df %>%
  left_join(fusion_info_df, by = 'Kids_First_Biospecimen_ID') %>%
  # Replace NA values in 'breakpoint_group' and 'breakpoint_type' with 'non-fusion'
  mutate(breakpoint_group = ifelse(is.na(breakpoint_group), "non-fusion", breakpoint_group),
         breakpoint_type = ifelse(is.na(breakpoint_type), "non-fusion", breakpoint_type)) %>% 

  dplyr::select(Kids_First_Biospecimen_ID, plot_group, molecular_subtype_display, cluster_assigned,breakpoint_group, breakpoint_type) %>%
  filter(plot_group=="Low-grade glioma")


## plot LGGs with fusion
histologies_w_fusion_df <- histologies_w_fusion_df %>%  group_by(cluster_assigned, breakpoint_group) %>%
       summarise(count = n()) %>%
       mutate(total_in_cluster = sum(count), 
              percentage = (count / total_in_cluster) * 100) %>%
              arrange(cluster_assigned, desc(percentage))

pdf(file.path(plots_dir, "cluster_membership-subtypes-LGGs-fusion.pdf"), height = 6, width = 11)
ggplot(histologies_w_fusion_df, aes(x = factor(cluster_assigned), y = percentage, fill = breakpoint_group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Cluster Assigned", y = "Percentage", fill = "Breakpoint Group") +
  theme_Publication()
dev.off()

plot_enr <- function(df, var1, var2,
                     var1_names, var2_names,
                     padjust = FALSE) {
  
  enr <- matrix(0, nrow(unique(df[,var1])),
                nrow(unique(df[,var2])),
                dimnames = list(var1_names,
                                var2_names))
  pval <- enr
  ct <- enr
  
  # loop through cancer groups to calculate enrichment 
  for (i in 1:nrow(enr)){
    no_var1 <- sum(unlist(df[,var1]) == rownames(enr)[i])
    for (j in 1:ncol(enr)){
      no_var2 <- sum(unlist(df[,var2]) == colnames(enr)[j] & !is.na(unlist(df[,var2])))
      no_var1_var2 <- sum(unlist(df[,var1]) == rownames(enr)[i] & unlist(df[,var2]) == colnames(enr)[j])
      ct[i,j] <- no_var1_var2
      enr[i,j] <- (no_var1_var2/no_var2)/(no_var1/nrow(df))
      pval[i,j] <- phyper(no_var1_var2, no_var1, nrow(df) - no_var1, no_var2, lower.tail = F)
    }
  }
  
  if (padjust == TRUE) {
    
    fdr <- t(apply(pval, 1, function(x) p.adjust(x, "fdr")))
    sig_mat <- ifelse(fdr < 0.05 & enr > 1 & ct > 1, "*", "")
    
  } else {
    
    sig_mat <- ifelse(pval < 0.05 & enr > 1 & ct > 1, "*", "")
    
  }
  
  fill_mat <- matrix(glue::glue("{round(enr, 1)}{sig_mat}"), 
                     nrow(enr), ncol(enr))
  
  ct_enr_mat <- matrix(glue::glue("{ct}\n({fill_mat})"),
                       nrow(ct), ncol(ct))
  
  col_fun = colorRamp2(c(0, ceiling(max(enr))), c("white", "orangered"))
  
  region_ht <- Heatmap(enr,
                       name = "Odds ratio",
                       cluster_rows = F,
                       cluster_columns = F,
                       rect_gp = gpar(col = "black", lwd = 2),
                       col = col_fun,
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%s", ct_enr_mat[i, j]), x, y, gp = gpar(fontsize = 12))
                       },
                       column_names_side = "top",
                       column_names_rot = 25,
                       column_names_centered = FALSE)
  
  return(region_ht)
  
}

group_ht <- plot_enr(histologies_w_fusion_df, 
                     "cluster_assigned", "breakpoint_group",
                     var1_names = c("1", "2", "3", "6", "7", "8", "9", "10","11"),
                     var2_names = c("15:09", "16:09", "16:11", "18:10", "rare", "non-fusion"),
                     padjust = TRUE)

# plot enrichment results
pdf(file.path(plot_dir, "breakpoint_group_ancestry_ct_enr_heatmap.pdf"),
    height = 3, width = 4)

draw(group_ht)

invisible(dev.off())

draw(group_ht)
