# R Corbett 2023
#
# Plot splicing survival models

# Load libraries
library(tidyverse)
library(survival)
library(patchwork)
library(colorblindr)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "survival")

# call survival models script
source(file.path(analysis_dir, "util", "survival_models.R"))

# Declare output directory
plots_dir <- file.path(analysis_dir, "plots")

# Input directory
input_dir <- file.path(analysis_dir, "results")

subtype_df <- read_tsv(file.path(analysis_dir, "results", "subtypes-for-survival.tsv"))
# Define histology groups and subtypes
groups <- c(unique(subtype_df$foldername),
            glue::glue("{subtype_df$foldername}, {subtype_df$subtype}"))

dir_names <- c(unique(subtype_df$foldername),
               subtype_df$foldername)
names(dir_names) <- groups

# Define identifiers for model files
file_names <- c(unique(subtype_df$names),
                glue::glue("{subtype_df$names}_{subtype_df$subtype_name}"))
names(file_names) <- groups

# Loop through histology groups and subtypes to generate Kaplan-Meier plots from models

for (group in groups){
  
  input_dir <- file.path(analysis_dir, "results", dir_names[group])
  plots_dir <- file.path(analysis_dir, "plots", dir_names[group])
  
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
    
  }
  
  km_os_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_OS_SIburden.RDS")
    )
  )
  
  km_efs_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_EFS_SIburden.RDS")
    )
  )
  
  km_output_pdf <- file.path(plots_dir, glue::glue("km_{file_names[group]}_OS_EFS_SIburden.pdf"))
  
  km_plot <- suppressWarnings(
                plotKM(model = list(km_os_result, km_efs_result),
                    variable = "SI_group",
                    combined = T, 
                    title = group)
  )
  
  ggsave(km_output_pdf, km_plot,
         width = 10, height = 6, units = "in",
         device = "pdf")
  
  
}

# Loop through histologies and subtypes to generate forest plots from coxph models with SI_group as predictor

# Forest plots for OS

for (group in groups){
  
  input_dir <- file.path(analysis_dir, "results", dir_names[group])
  plots_dir <- file.path(analysis_dir, "plots", dir_names[group])
  
  if (grepl("GNG|LGG", group)){
    os_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_resection_SIburden.RDS")
      ))
  }else{
    os_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_SIburden.RDS")
      ))
  }
  
  os_forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[group]}_OS_subtype_SIburden.pdf"))
  
  os_forest_plot <- suppressWarnings(
                        plotForest(os_survival_result)
  )
  
  ggsave(os_forest_pdf, os_forest_plot, width = 8, height = 3)
  
  # Forest plots for EFS
  
  if (grepl("GNG|LGG", group)){
    efs_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_EFS_additive_terms_subtype_resection_SIburden.RDS")
      ))
  }else{
    efs_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_EFS_additive_terms_subtype_SIburden.RDS")
      ))
  }
  
  efs_forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[group]}_EFS_subtype_SIburden.pdf"))
  
  efs_forest_plot <- suppressWarnings(
                          plotForest(efs_survival_result)
  )
  
  ggsave(efs_forest_pdf, efs_forest_plot, width = 8, height = 3)
  
}

# Loop through histologies and subtypes to generate forest plots from coxph models with SI as predictor 

for (group in groups){
  
  input_dir <- file.path(analysis_dir, "results", dir_names[group])
  plots_dir <- file.path(analysis_dir, "plots", dir_names[group])
  
  # Forest plots for OS
  
  if (grepl("GNG|LGG", group)){
    os_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_resection_SI.RDS")
      ))
  }else{
    os_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_SI.RDS")
      ))
  }
  
  os_forest_pdf <- file.path(plots_dir, 
                             glue::glue("forest_{file_names[group]}_OS_subtype_SI.pdf"))
  
  os_forest_plot <- suppressWarnings(
                        plotForest(os_survival_result)
  )
  
  ggsave(os_forest_pdf, os_forest_plot, width = 8, height = 3)
  
  # Forest plots for PFS
  
  if (grepl("GNG|LGG", group)){
    efs_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_EFS_additive_terms_subtype_resection_SI.RDS")
      ))
  }else{
    efs_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_EFS_additive_terms_subtype_SI.RDS")
      ))
  }
  
  efs_forest_pdf <- file.path(plots_dir, 
                              glue::glue("forest_{file_names[group]}_EFS_subtype_SI.pdf"))
  
  efs_forest_plot <- suppressWarnings(
                        plotForest(efs_survival_result)
  )
  
  ggsave(efs_forest_pdf, efs_forest_plot, width = 8, height = 3)
  
}

