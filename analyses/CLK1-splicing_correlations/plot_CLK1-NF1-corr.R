################################################################################
# plot_CLK1-NF1-corr.R 
# Plot scatter and compute R for CLK1 expression and PSI values for HGGs
# written by Ammar S Naqvi, Jo Lynne Rokita
# Usage: Rscript plot_CLK1-NF1-corr.R
################################################################################

## load libraries
suppressPackageStartupMessages({
  library("ggplot2")
  library("tidyverse")
  library("ggpubr")
  library("vroom")
  library('data.table')
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing_correlations")
results_dir <- file.path(analysis_dir, "results")

plots_dir <- file.path(analysis_dir, "plots")
figures_dir <- file.path(root_dir, "figures")

if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

# theme for all plots
source(file.path(figures_dir, "theme_for_plots.R"))
source(file.path(analysis_dir, "util", "function-create-scatter-plot.R"))

## define input files
clin_file <- file.path(data_dir,"histologies-plot-group.tsv")
rsem_transc_counts <- file.path(data_dir,"rna-isoform-expression-rsem-tpm.rds")
rmats_clk1_file <- file.path(data_dir, "clk1-splice-events-rmats.tsv")
rmats_nf1_file <- file.path(results_dir, "nf1-splice-events-rmats.tsv")

## output file
plot_high_low_CLK1_file <- file.path(plots_dir,"high-low-CLK1_NF1.pdf")
plot_high_low_NF1_file <- file.path(plots_dir,"high-low-NF1_CLK1.pdf")
plot_scatter_file <- file.path(plots_dir,"scatter-CLK1-NF1-PSI.pdf")

## get CLK1 psi values in tumors and ctrls
indep_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
indep_df <- vroom(indep_file)

## read in histology file and count data
## filter histology file for all HGG, only stranded samples
histologies_df  <-  read_tsv(clin_file) %>%
  filter(cohort == "PBTA",
         experimental_strategy == "RNA-Seq",
         Kids_First_Biospecimen_ID %in% indep_df$Kids_First_Biospecimen_ID,
         RNA_library == 'stranded'
         )

hgg_bs_id <- histologies_df %>%
  filter(plot_group %in% c("DIPG or DMG", "Other high-grade glioma")) 

clk1_DMG_rmats <- histologies_df %>%
  filter(plot_group == "DIPG or DMG")

other_hgg_bs_id <- histologies_df %>%
  filter(plot_group == "Other high-grade glioma") 


NF1_rmats <- fread(rmats_nf1_file) %>%
  dplyr::filter(sample_id %in% hgg_bs_id$Kids_First_Biospecimen_ID,
         geneSymbol=="NF1",
         splicing_case =='SE',
         exonStart_0base=='31252937', 
         exonEnd=='31253000') %>% 
  dplyr::select(sample_id, IncLevel1)

## load rmats input for CLK1
clk1_DMG_rmats <- fread(rmats_clk1_file) %>%
  # filter for CLK1 and exon 4
  filter(sample_id %in% hgg_bs_id$Kids_First_Biospecimen_ID,
         geneSymbol=="CLK1",
         exonStart_0base=="200860124", 
         exonEnd=="200860215") %>% 
  # select minimal info
  select(sample_id, IncLevel1) %>%
  unique()

clk1_OtherHGG_rmats <- fread(rmats_clk1_file) %>%
  # filter for CLK1 and exon 4
  filter(sample_id %in% other_hgg_bs_id$Kids_First_Biospecimen_ID,
         geneSymbol=="CLK1",
         exonStart_0base=="200860124", 
         exonEnd=="200860215") %>% 
  # select minimal info
  select(sample_id, IncLevel1) %>%
  unique()

clk1_HGG_rmats <- fread(rmats_clk1_file) %>%
  # filter for CLK1 and exon 4
  filter(sample_id %in% hgg_bs_id$Kids_First_Biospecimen_ID,
         geneSymbol=="CLK1",
         exonStart_0base=="200860124", 
         exonEnd=="200860215") %>% 
  # select minimal info
  select(sample_id, IncLevel1) %>%
  unique()



## CLK1 PSI correlations
CLK1_NF1_rmats_HGG_df <- inner_join(clk1_HGG_rmats, NF1_rmats_df, by= "sample_id",suffix = c("_CLK1", "_NF1"))

## remove non-splice events (to remove ones that are 1 and 0)
CLK1_NF1_rmats_HGG_filter_df <- CLK1_NF1_rmats_HGG_df %>%
  filter(IncLevel1_NF1 < .95 & IncLevel1_NF1 > 0.05,
         IncLevel1_CLK1 < .95 & IncLevel1_CLK1 > 0.05)
  
## scatter plot of CLK1 vs NF1 PSI
scatter_CLK1_vs_NF1_incl <- ggscatter(CLK1_NF1_rmats_HGG_filter_df, 
                                      x="IncLevel1_CLK1", 
                                      y="IncLevel1_NF1", 
                                      add = "reg.line", 
                                      conf.int = TRUE, 
                                      cor.coef = TRUE, 
                                      cor.method = "spearman",
                                      add.params = list(color = "red",
                                                        fill = "pink"),
                                      ticks = TRUE) + 
  xlab(expression(bold(bolditalic("CLK1")~"Exon 4 Inclusion (PSI)"))) +
  ylab(expression(bold(bolditalic("NF1")~"Exon 23a Inclusion (PSI)"))) +
  theme_Publication()  

# save plot 
pdf(plot_scatter_file, width = 4.5, height = 4.5)
print(scatter_CLK1_vs_NF1_incl)
dev.off()

## separate by high vs low Ex4 PSI
## Compute quantiles to define high vs low Exon 4 PSI groups
quartiles_CLK1_ex4_psi <- quantile(CLK1_NF1_rmats_HGG_df$IncLevel1_CLK1, probs=c(.25, .75), na.rm = FALSE)

# Calculate IQR
IQR_psi <- IQR(CLK1_NF1_rmats_HGG_df$IncLevel1_CLK1)

# Get lower quantile (25%)
lower_psi <- quartiles_psi[1] 

# Get upper quantile (75%)
upper_psi <- quartiles_psi[2]

CLK1_NF1_rmats_HGG_level_df <- CLK1_NF1_rmats_HGG_df %>% 
  mutate(PSI_CLK1 = case_when(IncLevel1_CLK1 > upper_psi ~ "high",
                              IncLevel1_CLK1 < lower_psi ~ "low",
                               TRUE ~ NA_character_)) %>% 
  filter(!is.na(PSI_CLK1)) 


## Make box plot with stats
boxplot_CLK_vs_NF1_incl<- ggboxplot(CLK1_NF1_rmats_HGG_level_df,x = "PSI_CLK1", y = "IncLevel1_NF1") +
  xlab(expression(bold(bolditalic("CLK1")~"Exon 4 PSI Level"))) +
  ylab(expression(bold("NF1 Exon 23a Inclusion"))) +
  lims(y = c(0,1.1)) +
  ggforce::geom_sina(aes(color = PSI_CLK1), size = 2, shape = 21, fill = NA, stroke = 1) +
  scale_color_manual(name = "CLK1 PSI Level", values = c(high = "#FFC20A", low = "#0C7BDC")) +
  stat_compare_means(position = "identity", label.x = 1) +
  theme_Publication()


# save plot 
pdf(plot_high_low_CLK1_file, width = 4.5, height = 4.5)
print(boxplot_CLK_vs_NF1_incl)
dev.off()

## same idea but separate by high vs low NF1 PSI
## Compute quantiles to define high vs low Exon 4 PSI groups
quartiles_NF1_ex4_psi <- quantile(CLK1_NF1_rmats_HGG_df$IncLevel1_NF1, probs=c(.25, .75), na.rm = FALSE)

# Calculate IQR
IQR_psi <- IQR(CLK1_NF1_rmats_HGG_df$IncLevel1_NF1)

# Get lower quantile (25%)
lower_psi <- quartiles_psi[1] 

# Get upper quantile (75%)
upper_psi <- quartiles_psi[2]

CLK1_NF1_rmats_HGG_level_df <- CLK1_NF1_rmats_HGG_level_df %>% 
  mutate(PSI_NF1 = case_when(IncLevel1_NF1 > upper_psi ~ "high",
                              IncLevel1_NF1 < lower_psi ~ "low",
                              TRUE ~ NA_character_)) %>% 
  filter(!is.na(PSI_NF1)) 


## Make box plot with stats
boxplot_NF1_vs_CLK1_incl <- ggboxplot(CLK1_NF1_rmats_HGG_level_df,x = "PSI_NF1", y = "IncLevel1_CLK1") +
  xlab(expression(bold(bolditalic("NF1")~"Exon 23a PSI Level"))) +
  ylab(expression(bold("CLK1 Exon 4 Inclusion"))) +
  lims(y = c(0,1.1)) +
  ggforce::geom_sina(aes(color = PSI_NF1), size = 2, shape = 21, fill = NA, stroke = 1) +
  scale_color_manual(name = "NF1 PSI Level", values = c(high = "#FFC20A", low = "#0C7BDC")) +
  stat_compare_means(position = "identity", label.x = 1) +
  theme_Publication()

# save plot 
pdf(plot_high_low_NF1_file, width = 4.5, height = 4.5)
print(boxplot_NF1_vs_CLK1_incl)
dev.off()

## Transcript expression based
rsem_skipped_NF1_exon_df <- readRDS(rsem_transc_counts) %>%
  filter(grepl("ENST00000356175", transcript_id))

rsem_included_NF1_exon_df <- readRDS(rsem_transc_counts) %>%
  filter(grepl("ENST00000358273", transcript_id)) # ENST00000358273

## transcript we saw in morpho vs untreated comparison
rsem_ENST0000047157_df <- readRDS(rsem_transc_counts) %>%
  filter(grepl("ENST0000047157", transcript_id))


rmats_exp_CLK1_NF1_incl_HGG_df <- rsem_included_NF1_exon_df %>%
  select(all_of(clk1_HGG_rmats$sample_id)) %>% 
  pivot_longer(cols = tidyselect::everything(),
               names_to=c("sample_id"), 
               values_to="Expr") %>% 
  inner_join(clk1_HGG_rmats, by= "sample_id") %>%
  # need to add this so that we can use this in the plot axis
  #mutate(transcript = paste0(transcript_id)) %>%
  dplyr::filter(Expr > 1)


set.seed(2023)

##log2 expression
expr_psi_df_log <- rmats_exp_CLK1_NF1_incl_HGG_df %>%
  mutate(logExp = log(Expr, 2))

## create plot
scatterplot_included <- ggscatter(expr_psi_df_log, 
                         x="IncLevel1", 
                         y="Expr", 
                         add = "reg.line", 
                         conf.int = TRUE, 
                         cor.coef = TRUE, 
                         cor.method = "spearman",
                         add.params = list(color = "red",
                                           fill = "pink"),
                         ticks = TRUE) + 
  xlab(expression(bold(bolditalic("CLK1")~"Exon 4 Inclusion (PSI)"))) +
  ylab(substitute(bold("23a Included Transcript (RSEM)"))) +
  theme_Publication()  

# save plot 
pdf(file.path(plots_dir,"scatter-CLK1-PSI-ex23a-transcript.pdf"), width = 4.5, height = 4.5)
print(scatterplot_included)
dev.off()

## compare transcript expression of morpholiono included exon transcript (NMD-substrate) vs NF1 e23a included exon  
NF1_included_exon <- rsem_included_NF1_exon_df %>%
  select(all_of(clk1_HGG_rmats$sample_id)) %>% 
  pivot_longer(cols = tidyselect::everything(),
               names_to=c("sample_id"), 
               values_to="Expr") 

NF1_morpho_incl_exon_df <- rsem_ENST0000047157_df %>%
  select(all_of(clk1_HGG_rmats$sample_id)) %>% 
  pivot_longer(cols = tidyselect::everything(),
               names_to=c("sample_id"), 
               values_to="Expr") %>%
  inner_join(NF1_included_exon, by= "sample_id", suffix = c("_morpho","_ex23a" )) %>%
  # need to add this so that we can use this in the plot axis
  #mutate(transcript = paste0(gene)) %>%
  dplyr::filter(Expr_morpho > 1,
                Expr_ex23a > 1)

scatterplot_morpho_vs_incl <- ggscatter(NF1_morpho_incl_exon_df, 
                                  x="Expr_morpho", 
                                  y="Expr_ex23a", 
                                  add = "reg.line", 
                                  conf.int = TRUE, 
                                  cor.coef = TRUE, 
                                  cor.method = "spearman",
                                  add.params = list(color = "red",
                                                    fill = "pink"),
                                  ticks = TRUE) + 
  xlab(expression(bold("Morpholino Transcript (RSEM)"))) +
  ylab(substitute(bold("23a Included Transcript (RSEM)"))) +
  theme_Publication()  


# save plot 
pdf(file.path(plots_dir,"scatter-morpho-ex23a-transcript.pdf"), width = 4.5, height = 4.5)
print(scatterplot_morpho_vs_incl)
dev.off()


