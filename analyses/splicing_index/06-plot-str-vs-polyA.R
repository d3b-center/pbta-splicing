################################################################################
# 06-plot-str-vs-polyA.R
# 
# Plot scatter plot to assess possible differences in PSI calculations in polyA
# vs non-polyA 
# 
# written by Ammar Naqvi
#
# usage: Rscript 06-plot-str-vs-polyA.R
################################################################################

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_index")
plots_dir <- file.path(analysis_dir, "plots")

# Specify file paths
clin_file  <- file.path(data_dir,"histologies.tsv")
rmats_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")

## ouput
PT_RYMG3M91_scatter_path <- file.path(plots_dir,"PT_RYMG3M91_scatter.pdf")
PT_W5GP3F6B_scatter_path <- file.path(plots_dir,"PT_W5GP3F6B_scatter.pdf")


## Store data in dataframes
clin_df <- vroom(clin_file) %>% 
  filter(experimental_strategy=='RNA-Seq')

rmats_df <- fread(rmats_file)

## General plots assessing PSI in poly A vs Str for known samples PT_RYMG3M91 and PT_W5GP3F6B
clin_to_check_df <- clin_df %>% filter(Kids_First_Participant_ID=='PT_RYMG3M91' |
                                         Kids_First_Participant_ID=='PT_W5GP3F6B')

# Select two sample_ids for comparison
sample1_id <- "BS_68KX6A42"
sample2_id <- "BS_D7XRFE0R"

rmats_df <- fread(rmats_file)
psi_subset_df <- fread(rmats_file) %>% filter(sample_id==sample1_id| 
                                                  sample_id==sample2_id) 

psi_PT_RYMG3M91 <- psi_subset_df %>% filter(splicing_case=='SE') %>%
  dplyr::mutate(SpliceID = paste(geneSymbol, exonStart_0base, exonEnd, upstreamES,upstreamEE,downstreamES,downstreamEE, sep = ":") ) %>%
  select(sample_id,SpliceID,IncLevel1)

# Filter data for the selected sample_ids
sample1_data <- psi_PT_RYMG3M91 %>%
  filter(sample_id == sample1_id)
sample2_data <- psi_PT_RYMG3M91 %>%
  filter(sample_id == sample2_id)

# Merge data for the same SpliceID
merged_data <- merge(sample1_data, sample2_data, by = "SpliceID", suffixes = c("_sample1", "_sample2"))

# Create scatterplot using ggscatter
PT_RYMG3M91_scatter <- ggscatter(merged_data, x = "IncLevel1_sample1", y = "IncLevel1_sample2", 
               xlab = paste("Total PSI of", sample1_id,"(poly-A)"),
               ylab = paste("Total PSI of", sample2_id,"(stranded)"),
               title = "polyA vs stranded comparison of PSI ",
               add = "reg.line", 
               conf.int = TRUE, 
               cor.coef = TRUE, 
               cor.method = "pearson",
               add.params = list(color = "red",
                                 fill = "pink"),
               ticks = TRUE) + 
               theme_Publication()
              

# Save plot as pdf
pdf(PT_RYMG3M91_scatter_path, height = 4, width = 6, useDingbats = FALSE)
print(PT_RYMG3M91_scatter)
dev.off()

# Select two sample_ids for comparison from PT_W5GP3F6B (these are known same aliquots)
sample1_id <- "BS_7WM3MNZ0"
sample2_id <- "BS_KABQQA0T"

##  same with PT_W5GP3F6B
psi_PT_W5GP3F6B <- rmats_df %>% filter(sample_id==sample1_id | 
                                                  sample_id==sample2_id) %>% 
                                filter(splicing_case=='SE') %>%
                                dplyr::mutate(SpliceID = paste(geneSymbol, exonStart_0base, exonEnd, upstreamES,upstreamEE,downstreamES,downstreamEE, sep = ":") ) %>%
                                select(sample_id,SpliceID,IncLevel1)



# Filter data for the selected sample_ids
sample1_data <- psi_PT_W5GP3F6B %>%
  filter(sample_id == sample1_id)
sample2_data <- psi_PT_W5GP3F6B %>%
  filter(sample_id == sample2_id)

# Merge data for the same SpliceID
merged_data <- merge(sample1_data, sample2_data, by = "SpliceID", suffixes = c("_sample1", "_sample2"))

PT_W5GP3F6B_scatter <- ggscatter(merged_data, x = "IncLevel1_sample1", y = "IncLevel1_sample2", 
               xlab = paste("IncLevel1 for", sample1_id),
               ylab = paste("IncLevel1 for", sample2_id),
               title = "polyA vs stranded comparison of PSI ",
               add = "reg.line", 
               conf.int = TRUE, 
               cor.coef = TRUE, 
               cor.method = "pearson",
               add.params = list(color = "red",
                                 fill = "pink",
                                 alpha=0.05),
               ticks = TRUE) + 
              theme_Publication()

# Save plot as pdf
pdf(PT_W5GP3F6B_scatter_path, height = 4, width = 6, useDingbats = FALSE)
print(PT_W5GP3F6B_scatter)
dev.off()
