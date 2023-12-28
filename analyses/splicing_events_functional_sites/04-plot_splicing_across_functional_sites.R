################################################################################
# 04-plot_splicing_across_functional_sites.R
# script that plots a table of splicing events overlapping uniprot sites
# written by Ammar Naqvi
#
# usage: Rscript 04-plot_splicing_across_functional_sites.R 
################################################################################

suppressPackageStartupMessages({
  library("ggstatsplot")
  library("vroom")
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("clusterProfiler")
  library("msigdbr")
  library("org.Hs.eg.db")
  library("ggrepel")
  library("cowplot")
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "splicing_events_functional_sites")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

##theme for all plots
# source function for theme for plots survival
figures_dir <- file.path(root_dir, "figures")
source(file.path(figures_dir, "theme_for_plots.R"))

## output files for final plots
file_dpsi_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites.HGG.pdf")
file_dpsi_kinase_plot <- file.path(analysis_dir, "plots", "dPSI_across_functional_sites_kinase.HGG.pdf")

## retrieve psi values from tables
file_psi_pos_func <-  file.path(results_dir,"splicing_events.total.HGG.pos.intersectUnip.ggplot.txt")
file_psi_neg_func <-  file.path(results_dir,"splicing_events.total.HGG.neg.intersectUnip.ggplot.txt")

## read table of recurrent functional splicing (skipping)
dpsi_unip_pos <- vroom(file_psi_pos_func) %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>%
  mutate(Preference='Inclusion')

## read table of recurrent functional splicing (inclusion) 
dpsi_unip_neg <- vroom(file_psi_neg_func) %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>% 
  mutate(Preference='Skipping')

psi_comb <- rbind(dpsi_unip_neg,dpsi_unip_pos)

## ggstatplot across functional sites
set.seed(123)
plot_dsp <-  ggplot(psi_comb,aes(Uniprot, dPSI*100) ) +  
  ylab(expression(bold("dPSI"))) +
  ggforce::geom_sina(aes(color = Preference, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  facet_wrap("Preference") +
  stat_compare_means() + 

  #ggforce::geom_sina(aes(color = Preference), size = 2,method="density") +
  scale_color_manual(name = "Preference", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A"))  + 
  theme_Publication() + 
  labs(y=expression('Percent Spliced In'), x= "Uniprot-defined Functional Site") + 
  theme(legend.position="none")

# Save plot as PDF
pdf(file_dpsi_plot, 
    width = 10, height = 4)
plot_dsp
dev.off()


# kinase gene list
known_kinase_file <- file.path(input_dir,'kinase_known.txt')
known_kinase_df <- vroom(known_kinase_file, delim = "\t", col_names = 'gene') 
psi_unip_kinase <- dplyr::inner_join(psi_comb, known_kinase_df, by='gene') 

## make sina plot
set.seed(45)
kinase_dpsi_plot <- ggplot(psi_unip_kinase,aes(Preference,dPSI*100) ) +  
  ylab(expression(bold("dPSI"))) +
  ggforce::geom_sina(aes(color = Preference, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  ggforce::geom_sina(aes(color = Preference), size = 2,method="density") +
  geom_label_repel(box.padding = 0.5, min.segment.length = 0.5,max.overlaps =Inf, aes(label = gene), data=psi_unip_kinase %>% subset(gene %in% c("CLK1")), size=2) +
  scale_color_manual(name = "Preference", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A")) + 
  theme_Publication() +
  labs(y=expression('Percent Spliced In')) + 
  
  theme(legend.position="none")

pdf(file_dpsi_kinase_plot, 
    width = 4, height = 4)
kinase_dpsi_plot
dev.off()

## over-repr analysis for kinases
# plot output path
ora_dotplot_path = file.path(plots_dir,"kinases-ora-plot.pdf")

# get gene sets relevant to H. sapiens
hs_msigdb_df <- msigdbr(species = "Homo sapiens")

## filter for kegg pathways that are included in the curated gene sets
hs_kegg_df <- hs_msigdb_df %>%
  dplyr::filter(gs_cat == "C2", gs_subcat == "CP:KEGG"
  )

hs_hm_df <- hs_msigdb_df %>%
  dplyr::filter(gs_cat == "H"  )

## skipping vs incl 
kinase_skip_pref <- psi_unip_kinase %>% dplyr::filter(Preference=='Skipping') 
kinase_incl_pref <- psi_unip_kinase %>% dplyr::filter(Preference=='Inclusion') 

## hallmark gene set
ora_incl_hm_results <- enricher(
  gene = kinase_incl_pref$gene, # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    hs_hm_df,
    gs_name,
    human_gene_symbol
  )
)

ora_skip_hm_results <- enricher(
  gene = kinase_skip_pref$gene, # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    hs_hm_df,
    gs_name,
    human_gene_symbol
  )
)

##kegg gene set
ora_incl_kegg_results <- enricher(
  gene = kinase_incl_pref$gene, # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    hs_kegg_df,
    gs_name,
    human_gene_symbol
  )
)

ora_skip_kegg_results <- enricher(
  gene = kinase_skip_pref$gene, # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    hs_kegg_df,
    gs_name,
    human_gene_symbol
  )
)

enrich_skip_hm_plot <- enrichplot::dotplot(ora_skip_hm_results) + 
  theme_Publication() 
enrich_incl_hm_plot <- enrichplot::dotplot(ora_incl_hm_results) + 
  theme_Publication() 

enrich_skip_kegg_plot <- enrichplot::dotplot(ora_skip_kegg_results) + 
  theme_Publication() 
enrich_incl_kegg_plot <- enrichplot::dotplot(ora_incl_kegg_results) + 
  theme_Publication() 

## since kegg pathways have and show more enrichment, we will save those
plot_pathways_kegg <- plot_grid(enrich_skip_kegg_plot,enrich_incl_kegg_plot,align="hv",labels = c('Skipping', 'Inclusion'))
plot_pathways_kegg

## save ORA dotplot as tiff
ggplot2::ggsave(ora_dotplot_path,
                width=20,
                height=6,
                device="pdf")

