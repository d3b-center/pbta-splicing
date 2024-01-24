suppressPackageStartupMessages({
  library("tidyverse")
  library("optparse")
  library("vroom")
  library("easyGgplot2")
  library(gridExtra)
  library(grid)
  library(ggvenn)
} )

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "grants-related")

input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plots_dir   <- file.path(analysis_dir, "plots")

histology_file <- file.path(input_dir,"v19_plus_20210311_pnoc_rna.primary.tsv")
prelim_file    <- file.path(results_dir,"splicing_neoepitope.pri.dpsi.se.tsv")
summ_dmg_file  <-  file.path(input_dir,"summary_table.tsv")
neoepi_file    <- file.path(results_dir,"splicing_neoepitope.primary.dpsi.se.unipLocExtra.summary.txt")


histology_df <- vroom(histology_file,col_names = FALSE) %>% 
  dplyr::rename('BSID'='X1', 'PatientID'='X4')

prelim_df  <-  vroom(prelim_file, comment = "#",delim="\t", col_names = TRUE, trim_ws = TRUE, show_col_types = FALSE) %>% 
  dplyr::rename(splice_id=Splice_ID)

summ_dmg_df  <-  vroom(summ_dmg_file, comment = "#",delim="\t", col_names = TRUE, trim_ws = TRUE, show_col_types = FALSE)

## code to get H3 subsets but did not finish completely
H3_wt_df <- summ_dmg_df %>% filter(H3_WT=="Y") %>% mutate(PatientID= paste0("PT_",str_split(match_id, "_", simplify = TRUE)[,2]))
H3_1_k28m_df <- summ_dmg_df %>% filter(H3.1_K28M_sum=="Y") %>% mutate(PatientID= paste0("PT_",str_split(match_id, "_", simplify = TRUE)[,2]))
H3_3_k28m_df <- summ_dmg_df %>% filter(H3.3_K28M_sum=="Y") %>% mutate(PatientID= paste0("PT_",str_split(match_id, "_", simplify = TRUE)[,2]))
histology_wt_df <- histology_df %>% inner_join(H3_wt_df, by="PatientID")
histology_H31_k28m <-  histology_df %>% inner_join(H3_1_k28m_df, by="PatientID")
histology_H33_k28m <- histology_df %>% inner_join(H3_3_k28m_df, by="PatientID")


neoepi_df  <-  vroom(neoepi_file, comment = "#",delim="\t", col_names = TRUE, trim_ws = TRUE, show_col_types = FALSE)
neoepi_skipping <- neoepi_known_df %>% filter(dPSI_mean>0)
neoepi_incl <- neoepi_known_df %>% filter(dPSI_mean<0)

neo_skip_scatter = ggplot(neoepi_skipping, aes(x= freq  , y = log2(TPM_mean_tumor), fill = dPSI_mean ))+
  geom_point(color = 'black',size = 3, pch=21)+
  scale_fill_continuous(low="lightblue", high="blue") + 
  xlab("Recurrence") + ylab("log2 (TPM)") + 
  ggtitle("HGG Neoepitope candidates (skipping)") + theme_Publication()

neo_incl_scatter = ggplot(neoepi_incl, aes(x= freq  , y = log2(TPM_mean_tumor), fill = dPSI_mean ))+
  geom_point(color = 'black',size = 3, pch=21)+
  scale_fill_continuous(low="red", high="darkred") + 
  xlab("Recurrence") + ylab("log2 (TPM)") + 
  ggtitle("HGG Neoepitope candidates (inclusion)") + theme_Publication()

## arrange all plots in one grid
tiff(file.path(plots_dir, "neoepitopes_scatter_total.tiff"), height = 1800, width = 2000, units = "px", res = 300)
grid.arrange(neo_skip_scatter,neo_incl_scatter, ncol = 1, nrow = 2) 
dev.off()


## venn diagram of mis-spliced genes
neoepi_incl     <- neoepi_df %>% filter(dPSI_mean < 0 ) 
neoepi_skipping <- neoepi_df %>% filter(dPSI_mean > 0 ) 

neoepi_skipping_genes_list <- str_match(neoepi_skipping$splice_id,"(\\w+)\\_")[,2]
neoepi_incl_genes_list <- str_match(neoepi_incl$splice_id,"(\\w+)\\_")[,2]

#Make the plot
mis_spliced_genes_total <- list(
  Skipping = neoepi_skipping_genes_list, 
  Inclusion = neoepi_incl_genes_list)

tiff(file.path(plots_dir, "genes_spliced_total.tiff"), height = 2000, width = 2000, units = "px", res = 300)
ggvenn(mis_spliced_genes_total, fill_color = c("lightblue", "red"), stroke_size = .5)
dev.off()


# Theme for all plots
theme_Publication <- function(base_size=12, base_family="Arial") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            # legend.margin = unit(0.5, "cm"),
            legend.margin = margin(5,5,5,5),
            legend.title = element_text(face="bold"),
            plot.margin=unit(c(10,5,5,10),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

