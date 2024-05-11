################################################################################
# 00-get-splice-transcripts
# written by Jo Lynne Rokita
#
# This script identifies the NF1 transcripts with splice events from the 
# morpholino experiment
#
# usage: Rscript --vanilla 00-get-splice-transcripts.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "CLK1-splicing-impact-morpholino")
results_dir <- file.path(analysis_dir, "results")

## create plots dir if it doesn't exist
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

# files
gtf_file <- file.path(root_dir, "data", "gencode.v39.primary_assembly.annotation.gtf.gz")
morph_results_file <- file.path(data_dir, "morpholno.merged.rmats.tsv")

# Import GTF file
gtf_df <- import(gtf_file) %>%
  as.data.frame()

# select gene of interest, here, we want to get the NF1 transcripts - can make a list
goi <- "NF1"
incl_level <- 0.1

# we are selecting only exons, since that is what we have from rMATs to match
exons <- gtf_df %>%
  filter(gene_name == goi) %>%
  mutate(transcript_name_meta = case_when(transcript_type == "protein_coding" ~ paste0(transcript_name, "_PC"),
                                          transcript_type == "processed_transcript" ~ paste0(transcript_name, "_PT"),
                                          transcript_type == "retained_intron" ~ paste0(transcript_name, "_RI"),
                                          transcript_type == "nonsense_mediated_decay" ~ paste0(transcript_name, "_NMD"),
                                          TRUE ~ NA_character_)) %>%
  filter(type == "exon") %>%
  select(start, end, transcript_name_meta, exon_number) %>%
  mutate(exon_number = as.numeric(exon_number)) %>%
  unique()
  
# do this for each gene of interest
for (each_gene in goi){
  
  # read in morph results for goi
  morph_results <- read_tsv(morph_results_file) %>%
    mutate(psi_diff = abs(IncLevelDifference)) %>%
    filter(geneSymbol == each_gene)
  
    # first, format exon of interest
  morp_exon_oi <- morph_results %>%
    filter(geneSymbol == each_gene) %>%
    mutate(start = ifelse(!is.na(exonStart_0base), exonStart_0base+1,riExonStart_0base+1), 
           end = ifelse(!is.na(exonEnd), exonEnd,riExonEnd),
           is_ri = ifelse(!is.na(riExonStart_0base), TRUE, FALSE)) %>%
    select(splicing_case, start, end, IncLevelDifference, psi_diff, is_ri) %>%
    left_join(exons, relationship = "many-to-many") %>%
    mutate(riExonStart_0base = ifelse(is_ri == TRUE, start-1, NA_integer_), 
           exonStart_0base = ifelse(is_ri == FALSE, start-1, NA_integer_), 
           riExonEnd = ifelse(is_ri == TRUE, end, NA_integer_), 
           exonEnd = ifelse(is_ri == FALSE, end, NA_integer_)) %>%
    select(-start, -end, -is_ri) %>%
    dplyr::rename(exon_number_eoi = exon_number) 
  
  # next, format upstream exon
  morp_exon_upstream <- morph_results %>%
    filter(geneSymbol == each_gene) %>%
    mutate(start = upstreamES+1, 
           end = upstreamEE) %>%
    select(splicing_case, start, end, IncLevelDifference, psi_diff,) %>%
    left_join(exons, relationship = "many-to-many") %>%
    mutate(upstreamES = start-1,
           upstreamEE = end) %>%
    select(-start, -end) %>%
    dplyr::rename(exon_number_upstream = exon_number)
  
  # next, downstream exon
  morp_exon_downstream <- morph_results %>%
    filter(geneSymbol == each_gene) %>%
    mutate(start = downstreamES+1, 
           end = downstreamEE) %>%
    select(splicing_case, start, end, IncLevelDifference, psi_diff,) %>%
    left_join(exons, relationship = "many-to-many") %>%
    mutate(downstreamES = start-1,
           downstreamEE = end) %>%
    select(-start, -end) %>%
    dplyr::rename(exon_number_downstream = exon_number)
  
  # finally, join them back with each other and the morph results for that gene
  morp_exons <- morp_exon_upstream %>%
    right_join(morp_exon_oi) %>%
    right_join(morph_results) %>%
    left_join(morp_exon_downstream) %>%
    mutate(is_sequential = (exon_number_eoi == exon_number_upstream + 1) & (exon_number_downstream == exon_number_eoi + 1)) %>%
    # filter for only those with sequential exons or NA (in the case of RI)
    filter(psi_diff >= incl_level) %>%
    write_tsv(file.path(paste0(results_dir, "/morpholino_splice_", each_gene, "_", incl_level, "_psi_diff_transcripts.tsv")))
}

  