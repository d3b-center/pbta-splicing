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
plots_dir <- file.path(analysis_dir, "plots")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## output files for final plots
file_circos_plot <- file.path(analysis_dir, "plots", "cohort_circos.pdf")


# Load datasets and pre-process
hist_df <- read_tsv(file.path(data_dir,"histologies.tsv"), guess_max = 100000) %>% 
  # filter
  filter(experimental_strategy == "RNA-Seq",
         cohort == "PBTA",
         !is.na(pathology_diagnosis),
         composition != "Derived Cell Line") %>%
  # collapse reported gender to 3 groups
  mutate(reported_gender = case_when(reported_gender == "Not Reported" ~ "Unknown",
                                     TRUE ~ reported_gender),
         cancer_group = case_when(broad_histology %in% c("Non-tumor", "Benign Tumor") ~ NA_character_,
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
independent_specimens_df <- read_tsv(file.path(data_dir,"independent-specimens.rnaseqpanel.primary-plus.tsv")) %>%
  filter(cohort == "PBTA",
         experimental_strategy == "RNA-Seq")

# Merge both meta datasets
hist_indep_df <- combined_hist_map %>%
  right_join(independent_specimens_df, by="Kids_First_Biospecimen_ID") %>% 
  unique() %>%
  # Columns of interest: Kids_First_Biospecimen_ID, plot group, CNS_region, OS status, reported_gender
  dplyr::select(Kids_First_Biospecimen_ID, plot_group, CNS_region, reported_gender) %>%
  remove_rownames() %>% 
  column_to_rownames("Kids_First_Biospecimen_ID") %>%
  arrange(plot_group, CNS_region, reported_gender) 
  
# Create histology factors for splitting the Circos plot
split <- factor(hist_indep_df$plot_group)

# color palette for cancer group
cols <- unique(map_file$plot_group_hex)
names(cols) <- unique(as.character(map_file$plot_group))
sorted_cols <- cols[order(names(cols))]

gender_cols <- c( "deeppink4", "navy", "lightgrey")
names(gender_cols) <- c("Female", "Male", "Unknown")

loc_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#882255", "#6699CC")
names(loc_cols) <- c(sort(unique(hist_indep_df$CNS_region)))

#os_cols <- c("red", "green","black")
#names(os_cols) <- c("DECEASED", "LIVING", "NA")

cols_all <- c(cols, loc_cols, gender_cols)

pdf(file_circos_plot, width = 6.5, height = 6.5)
# Reset Circos plot
circos.clear()

# Make base Circos plot
circos.par(points.overflow.warning = FALSE)
circos.heatmap(hist_indep_df,
               split=split, # Specify how to split plot into sectors
               col = unlist(cols_all))#, # Add in color list
        
               #show.sector.labels = T, # To show labels of the main sectors of plot
               #cell.border = "black", # To show individual sample cells
               #track.height = 0.3, # Adjust the width of the plot
      
# add border colors to sectors
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = 1)
  circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1],
              CELL_META$cell.xlim[2], CELL_META$cell.ylim[2],
              col = NA, border = unlist(cols_all)[sn])
}

## adding legends
lgd_plot_group = Legend(title = "Histology", 
                        at =  names(sorted_cols), 
                        legend_gp = gpar(fill = sorted_cols))

lgd_tum_loc = Legend(title = "Tumor location", 
                       at =  names(loc_cols),
                       legend_gp = gpar(fill = loc_cols))

lgd_gender = Legend(title = "Sex", 
                    at =  names(gender_cols), 
                    legend_gp = gpar(fill = gender_cols))

#lgd_os = Legend(title = "OS Status", 
 #                   at =  names(os_cols), 
  #                  legend_gp = gpar(fill = os_cols))


## Placement of legend based on device size
h = dev.size()[2]
circle_size = unit(1, "snpc")

# Create legends
lgd_list1 = packLegend(lgd_plot_group, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 = packLegend(lgd_tum_loc, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list3 = packLegend(lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal")
#lgd_list4 = packLegend(lgd_os, max_height = unit(0.5*h, "inch"), direction = "horizontal")

# Add legends on plot
draw(lgd_list1, x = unit(75, "mm"), y = unit(85, "mm"))
draw(lgd_list2, x = unit(115, "mm"), y = unit(90, "mm"))
draw(lgd_list3, x = unit(111, "mm"), y = unit(56, "mm"))
#draw(lgd_list4, x = unit(114, "mm"), y = unit(85, "mm"))
circos.clear()
dev.off()



