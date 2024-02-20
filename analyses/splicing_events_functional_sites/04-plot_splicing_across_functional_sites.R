################################################################################
# 04-plot_splicing_across_functional_sites.R
# script that plots a table of splicing events overlapping uniprot sites
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 04-plot_splicing_across_functional_sites.R 
################################################################################

suppressPackageStartupMessages({
  library("ggpubr")
  library("vroom")
  library("ggplot2")
  library("dplyr")
  library("tidyverse")
  library("clusterProfiler")
  library("msigdbr")
  library("org.Hs.eg.db")
  library("ggrepel")
  library("cowplot")
  library("ggpubr")
  library("annoFuseData")
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
ora_dotplot_path <- file.path(plots_dir,"kinases-ora-plot.pdf")
kinases_functional_sites_skipped = file.path(results_dir,"kinases-functional_sites_skipped.txt")
kinases_functional_sites_incl = file.path(results_dir,"kinases-functional_sites_included.txt")

## retrieve psi values from tables
file_psi_pos_func <- file.path(results_dir,"splicing_events.total.HGG.pos.intersectUnip.ggplot.txt")
file_psi_neg_func <- file.path(results_dir,"splicing_events.total.HGG.neg.intersectUnip.ggplot.txt")

## read table of recurrent functional splicing (skipping)
dpsi_unip_pos <- vroom(file_psi_pos_func) %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>%
  mutate(Preference='Inclusion')

## read table of recurrent functional splicing (inclusion) 
dpsi_unip_neg <- vroom(file_psi_neg_func) %>% 
  mutate(gene=str_match(SpliceID, "(\\w+[\\.\\d]*)\\:")[, 2]) %>% 
  mutate(Preference='Skipping')

psi_comb <- rbind(dpsi_unip_neg,dpsi_unip_pos) %>% 
  mutate(Uniprot = case_when(Uniprot == 'DisulfBond' ~ "Disulfide Bond",
                             Uniprot == 'LocSignal' ~ "Localization Signal",
                             .default = Uniprot)
         )



## ggstatplot across functional sites
set.seed(123)
counts_psi_comb <- psi_comb %>% count(Preference, Uniprot )
plot_dsp <-  ggplot(psi_comb,aes(Uniprot, dPSI*100) ) +  
  ylab(expression(bold("dPSI"))) +
  ggforce::geom_sina(aes(color = Preference, alpha = 0.4), pch = 16, size = 5, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  facet_wrap("Preference") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Disulfide Bond", "Localization Signal"),
                                        c("Disulfide Bond", "Modifications"),
                                        c("Disulfide Bond", "Other"),
                                        c("Localization Signal", "Modifications"),
                                        c("Localization Signal", "Other"),
                                        c("Modifications", "Other"))) + 
  scale_color_manual(name = "Preference", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A"))  + 
  theme_Publication() + 
  
  labs(y="Percent Spliced In (PSI)", x= "Uniprot-defined Functional Site") + 
  geom_text(data = counts_psi_comb, aes(label = paste("n =",n), x = Uniprot, y = 0), vjust = 3, size = 4, hjust=.5) +
  theme(legend.position="none") +
  ylim(c(-1,170))

# Save plot as PDF
pdf(file_dpsi_plot, 
    width = 12, height = 4)
print (plot_dsp)
dev.off()

# get and filter for kinase genes
known_kinase_df <-read.delim(system.file("extdata", "genelistreference.txt", package = "annoFuseData")) %>%
  dplyr::rename(gene=Gene_Symbol) %>% 
  dplyr::filter(type=='Kinase')

psi_unip_kinase <- dplyr::inner_join(psi_comb, known_kinase_df, by='gene') 
counts_psi_unip_kinase <- psi_unip_kinase %>% count(Preference )

## make sina plot
set.seed(45)
kinase_dpsi_plot <- ggplot(psi_unip_kinase,aes(Preference,dPSI*100) ) +  
  ylab(expression(bold("dPSI"))) +
  ggforce::geom_sina(aes(color = Preference, alpha = 0.4), pch = 16, size = 4, method="density") +
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, coef = 0, aes(alpha = 0.4)) +
  geom_label_repel(box.padding = 0.5, min.segment.length = 0.5,max.overlaps =Inf, aes(label = gene), data=psi_unip_kinase %>% 
                     subset(gene %in% c("CLK1")), size=2) +
  scale_color_manual(name = "Preference", values = c(Skipping = "#0C7BDC", Inclusion = "#FFC20A")) + 
  theme_Publication() +
  labs(y="Percent Spliced In (PSI)") + 
  geom_text(data = counts_psi_unip_kinase, aes(label = paste("n =",n), x = Preference, y = 0), vjust = 3, size = 4, hjust=.5) +
  theme(legend.position="none")

pdf(file_dpsi_kinase_plot, 
    width = 4, height = 4)
kinase_dpsi_plot
dev.off()

## over-representation analysis for kinases
# get gene sets relevant to H. sapiens
hs_msigdb_df <- msigdbr(species = "Homo sapiens")

## filter for kegg pathways that are included in the curated gene sets
pathway_df <- hs_msigdb_df %>%
  dplyr::filter(gs_cat == "H" | gs_subcat %in% c("CP:KEGG", "CP:BIOCARTA", "TFT:GTRD"))

## skipping vs incl 
kinase_skip_pref <- psi_unip_kinase %>% 
  dplyr::filter(Preference=='Skipping') 
kinase_incl_pref <- psi_unip_kinase %>% 
  dplyr::filter(Preference=='Inclusion') 

##kegg gene set
ora_incl_results <- enricher(
  gene = unique(kinase_incl_pref$gene), # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)

ora_skip_results <- enricher(
  gene = unique(kinase_skip_pref$gene), # A vector of your genes of interest
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    pathway_df,
    gs_name,
    human_gene_symbol
  )
)


enrich_skip_plot <- enrichplot::dotplot(ora_skip_results, showCategory = 15) + 
  theme_Publication() +
  scale_color_gradient(name = "Adjusted p-value", 
                       low = "orange", high = "#0C7BDC") +  # Modify color range
  labs(color = "B-H adj p-value")  # Modify legend title 
enrich_incl_plot <- enrichplot::dotplot(ora_incl_results, showCategory = 15) +
  theme_Publication() +
  scale_color_gradient(name = "Adjusted p-value", 
                       low = "orange", high = "#0C7BDC") +  # Modify color range
  labs(color = "B-H adj p-value")  # Modify legend title

## since kegg pathways have and show more enrichment, we will save those
plot_pathways <- plot_grid(enrich_skip_plot,enrich_incl_plot,align="hv",
                                labels = c('Exon Skipping', 'Exon Inclusion'))
## save ORA dotplot
ggplot2::ggsave(filename = ora_dotplot_path,
                plot = plot_pathways, 
                width=17,
                height=7,
                device="pdf")

## write kinase results for table
write_lines(sort(unique(kinase_skip_pref$gene)), kinases_functional_sites_skipped)
write_lines(sort(unique(kinase_incl_pref$gene)), kinases_functional_sites_incl)

