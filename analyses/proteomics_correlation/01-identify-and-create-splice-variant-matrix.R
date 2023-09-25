################################################################################
# 01-identify-and-create-splice-variant-matrix.R
#
# written by Ammar Naqvi
#
# usage: Rscript 01-identify-and-create-splice-variant-matrix.R
################################################################################

## libraries
library("vroom")
library("tidyverse")

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

##directory setup
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "proteomics_correlation")
results_dir <- file.path(analysis_dir, "results")

plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")
output_dir <- file.path(analysis_dir, "output")

plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

## input files
splice_events_skipping_file <-"splicing_events.total.pos.intersectUnip.wo.total.txt"
splice_events_inclusion_file <- "splicing_events.total.neg.intersectUnip.wo.total.txt"
histology_file <- "histologies.tsv"
brain_goi_file <- "brain-goi-list-new.txt"
splice_events_file = "splice-events-rmats.tsv.gz"

histologies_df <- vroom(file.path(data_dir,histology_file)) %>% 
  filter(short_histology=='HGAT') 

splice_events_skip_df <- vroom(file.path(input_dir,splice_events_skipping_file), col_names = c("chr","start","end","SpliceID","dPSI","strand")) %>%
  select("chr","start","end","SpliceID","dPSI","strand") %>% 
  separate(SpliceID, c("gene", "exon coordinates"), sep="_", remove = FALSE)

splice_events_incl_df <- vroom(file.path(input_dir,splice_events_inclusion_file), col_names = c("chr","start","end","SpliceID","dPSI","strand")) %>%
  select("chr","start","end","SpliceID","dPSI","strand") %>% 
  separate(SpliceID, c("gene", "exon coordinates"), sep="_", remove = FALSE)

brain_goi_df <-  read.table(file.path(input_dir, brain_goi_file),header=TRUE)

splice_events_mixed_df <- inner_join(splice_events_skip_df,splice_events_incl_df, by='SpliceID') %>% 
  mutate(chr = coalesce(chr.x, chr.y),
         start = coalesce(start.x, start.y),
         end = coalesce(end.x, end.y),
         gene = coalesce(gene.x, gene.y),
         strand = coalesce(strand.x, strand.y)) %>% 
         rename("dPSI_Skip"=dPSI.x) %>% rename("dPSI_Incl"=dPSI.y) %>% 
  select(gene,chr,start,end,strand,SpliceID,dPSI_Skip,dPSI_Incl) %>% 
  distinct() %>% 
  inner_join(brain_goi_df, by='gene')


write_tsv(splice_events_mixed_df,file = file.path(results_dir,"mixed_events.tsv"), quote = 'none')

  





