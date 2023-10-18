################################################################################
# 01-generate-cohort-summary-circos-plot.R
# Generate circos summary plot of cohort used
# Author: Ammar Naqvi and Shehbeel Arif
# Usage: Rscript 01-generate-cohort-summary-circos-plot.R
################################################################################

# Load libraries
suppressPackageStartupMessages({
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
})

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "cohort_summary")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## output files for final plots
file_circos_plot <- file.path(analysis_dir, "plots", "cohort_circos.pdf")


# Load datasets and pre-process
hist_df <- vroom(file.path(data_dir,"histologies.tsv")) %>% 
  mutate(cancer_group = case_when(broad_histology %in% c("Non-tumor", "Benign Tumor") ~ NA_character_,
                                                                                     grepl("MPNST", pathology_diagnosis) ~ "Neurofibroma/Plexiform",
                                                                                     pathology_free_text_diagnosis == "meningioangiomatosis" ~ "Meningioma",
                                                                                     # move these from benign to other tumor
                                                                                     pathology_free_text_diagnosis %in% c("fibrous dysplasia",
                                                                                                                          "osteoblastoma",
                                                                                                                          "pituitary macroadenoma",
                                                                                                                          "pituitary adenoma",
                                                                                                                          "prolactinoma",
                                                                                                                          "perineuroma") ~ "Other tumor",
                                                                                     broad_histology == "Other tumor" ~ "Other tumor",
                                                                                     # pineoblastoma
                                                                                     pathology_free_text_diagnosis == "pnet - pineoblastoma with calcification" ~ "Pineoblastoma",
                                                                                     TRUE ~ as.character(cancer_group)))


# add cancer/plot group mapping file 
map_file <- read_tsv(file.path(input_dir, "plot-mapping.tsv")) %>%
  distinct()

# add plot mapping file and old plot groups
combined_hist_map <- hist_df %>%
  #left_join(hist_df) %>%
  left_join(map_file, by = c("broad_histology", "cancer_group")) 

## filter using independent specimens file
independent_specimens_df <- vroom(file.path(data_dir,"independent-specimens.rnaseqpanel.primary.tsv"))
independent_specimens_plus_df <- vroom(file.path(data_dir,"independent-specimens.rnaseqpanel.primary-plus.tsv"))
independent_specimens_rna_df <- rbind(independent_specimens_df,independent_specimens_plus_df) %>% 
  filter(cohort == "PBTA")

# Merge both meta datasets
hist_indep_df <- merge(combined_hist_map, independent_specimens_rna_df, by="Kids_First_Biospecimen_ID") %>% unique()

# Columns of interest: Kids_First_Biospecimen_ID, plot group, CNS_region, OS status, reported_gender
hist_indep_df <- hist_indep_df %>%
  dplyr::select(Kids_First_Biospecimen_ID, plot_group, CNS_region, OS_status, reported_gender) %>%
  remove_rownames() %>% 
  column_to_rownames("Kids_First_Biospecimen_ID") %>%
  arrange(plot_group, reported_gender) 
  

# Create histology factors for splitting the Circos plot
split <- factor(hist_indep_df$plot_group)

# color palette for cancer group
cols <- unique(map_file$plot_group_hex) 
names(cols) <- unique(as.character(map_file$plot_group))

gender_cols <- c("navy", "deeppink4")
names(gender_cols) <- c("Male", "Female")

os_cols <- c("red", "green","black")
names(os_cols) <- c("DECEASED", "LIVING", "NA")

cols_all <- c(cols, gender_cols, os_cols)

pdf(file_circos_plot, width = 6, height = 6)

# Reset Circos plot
circos.clear()

# Make base Circos plot
circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)
circos.heatmap(hist_indep_df,
               split=split, # Specify how to split plot into sectors
               col = unlist(cols_all), # Add in color list
               #show.sector.labels = T, # To show labels of the main sectors of plot
               #cell.border = "black", # To show individual sample cells
               track.height = 0.3, # Adjust the width of the plot
               
)


## adding legends
lgd_plot_group = Legend(title = "Plot group", 
                        at =  names(cols), 
                        legend_gp = gpar(fill = cols))

lgd_cns_group = Legend(title = "CNS region", 
                       at =  c("Posterior fossa","Other","Mixed","Hemispheric","Suprasellar","Spine","Midline","Ventricles","NA","Optic pathway"),
                       legend_gp = gpar(fill = 1:11))

lgd_gender = Legend(title = "Sex", 
                    at =  names(gender_cols), 
                    legend_gp = gpar(fill = gender_cols))

lgd_os = Legend(title = "OS Status", 
                    at =  names(os_cols), 
                    legend_gp = gpar(fill = os_cols))


## Placement of legend based on device size
h = dev.size()[2]

# Create legends
lgd_list1 = packLegend(lgd_plot_group, max_height = unit(0.5*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_cns_group, max_height = unit(0.5*h, "inch"), direction = "horizontal")
lgd_list3 = packLegend(lgd_gender, max_height = unit(0.5*h, "inch"), direction = "horizontal")
lgd_list4 = packLegend(lgd_os, max_height = unit(0.5*h, "inch"), direction = "horizontal")

# Add legends on plot
draw(lgd_list1, x = unit(110, "mm"), y = unit(118, "mm"), just = c("left", "top"))
draw(lgd_list2, x = unit(60, "mm"), y = unit(110, "mm"), just = c("right", "top"))
draw(lgd_list3, x = unit(225, "mm"), y = unit(110, "mm"), just = c("right", "top"))
draw(lgd_list4, x = unit(254, "mm"), y = unit(110, "mm"), just = c("right", "top"))
dev.off()



