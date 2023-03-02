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

# Define histology groups and subtypes
groups <- c("ATRT", "CPG", "EPN", "GNG",
            "HGG", "HGG WT", "HGG K28", 
            "LGG", "LGG WT", "LGG BRAF",
            "MB", "MB Group4", "MB SHH")

# Define identifiers for model files
file_names <- c("atrt", "cpg", "epn", 
                "gng", "hgg", "hgg_wildtype", 
                "hgg_DMG", 
                "lgg", "lgg_wildtype", "lgg_BRAF",
                "mb", "mb_Group4", "mb_SHH")
names(file_names) <- groups

# Loop through histology groups to generate Kaplan-Meier plots from models

for (group in groups){
  
  km_os_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_OS_SIburden.RDS")
    )
  )
  
  km_pfs_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_PFS_SIburden.RDS")
    )
  )
  
  km_output_pdf <- file.path(plots_dir, glue::glue("km_{file_names[group]}_OS_PFS_SIburden.pdf"))
  
  km_plot <- plotKM(model = list(km_os_result, km_pfs_result),
                    variable = "SI_group",
                    combined = T, 
                    title = group)
  
  ggsave(km_output_pdf, km_plot,
         width = 10, height = 6, units = "in",
         device = "pdf")
  
  
}

# Loop through broad histologies only (not subtypes) to generate forest plots from coxph models

groups <- groups[!grepl("HGG |LGG |MB ", groups)]

# Forest plots for OS

for (group in groups){
  
  # Forest plots for OS
  
  if (group == "LGG"){
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
  
  
  os_forest_plot <- plotForest(os_survival_result)
  
  ggsave(os_forest_pdf, os_forest_plot, width = 8, height = 3)
  
  # Forest plots for PFS
  
  if (group == "LGG"){
    pfs_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_PFS_additive_terms_subtype_resection_SIburden.RDS")
      ))
  }else{
    pfs_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_PFS_additive_terms_subtype_SIburden.RDS")
      ))
  }
  
  pfs_forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[group]}_PFS_subtype_SIburden.pdf"))
  
  
  pfs_forest_plot <- plotForest(pfs_survival_result)
  
  ggsave(pfs_forest_pdf, pfs_forest_plot, width = 8, height = 3)
  
  
}


# Define subtypes

subtypes <- c("HGG, H3 wildtype",
              "HGG, H3 K28",
              "LGG, wildtype",
              "LGG, BRAF mutant",
              "MB, Group4", 
              "MB, SHH")

file_names <- c("hgg_wildtype", "hgg_DMG", "lgg_wildtype", "lgg_BRAF",
                "mb_Group4", "mb_SHH")

names(file_names) <- subtypes


# Loop through subtypes to generate forest plots

for (subtype in subtypes){
  
  # OS
  
  if (grepl("LGG", subtype)){
    os_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[subtype]}_OS_additive_terms_subtype_resection_SIburden.RDS")
      ))
  }else{
    os_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[subtype]}_OS_additive_terms_subtype_SIburden.RDS")
      ))
  }
  
  os_forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[subtype]}_OS_SIburden.pdf"))
  
  
  forest_plot <- plotForest(os_survival_result)
  
  ggsave(os_forest_pdf, forest_plot, width = 8, height = 3)
  
  # PFS
  
  if (grepl("LGG", subtype)){
    pfs_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[subtype]}_PFS_additive_terms_subtype_resection_SIburden.RDS")
      ))
  }else{
    pfs_survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[subtype]}_PFS_additive_terms_subtype_SIburden.RDS")
      ))
  }
  
  pfs_forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[subtype]}_PFS_SIburden.pdf"))
  
  
  pfs_forest_plot <- plotForest(pfs_survival_result)
  
  ggsave(pfs_forest_pdf, pfs_forest_plot, width = 8, height = 3)
  
  
}



