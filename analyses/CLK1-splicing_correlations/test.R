# Load the package
library(rtracklayer)

# Import GTF file
# nf1 transcripts from GTF
gtf <- import(file.path(root_dir, "data", "gencode.v39.primary_assembly.annotation.gtf.gz")) %>%
  as.data.frame()

nf1 <- gtf %>%
  filter(gene_name == "NF1",
         type == "transcript") %>%
  mutate(transcript_name_meta = case_when(transcript_type == "protein_coding" ~ paste0(transcript_name, "_PC"),
                                          transcript_type == "processed_transcript" ~ paste0(transcript_name, "_PT"),
                                          transcript_type == "retained_intron" ~ paste0(transcript_name, "_RI"),
                                          transcript_type == "nonsense_mediated_decay" ~ paste0(transcript_name, "_NMD"),
                                          TRUE ~ NA_character_)) 

exon_23a_aka_31 <- gtf %>%
  filter(grepl("exon_number 31", V9) & grepl("ENST00000358273", V9))

head(nf1)
#


nf1_exons <- gtf %>%
  filter(gene_name == "NF1") %>%
  mutate(transcript_name_meta = case_when(transcript_type == "protein_coding" ~ paste0(transcript_name, "_PC"),
                                          transcript_type == "processed_transcript" ~ paste0(transcript_name, "_PT"),
                                          transcript_type == "retained_intron" ~ paste0(transcript_name, "_RI"),
                                          transcript_type == "nonsense_mediated_decay" ~ paste0(transcript_name, "_NMD"),
                                          TRUE ~ NA_character_)) %>%
  filter(type == "exon") %>%
  select(start, end, transcript_name_meta, exon_number) %>%
  mutate(exon_number = as.numeric(exon_number)) %>%
  unique()



morp_exon_oi <- read_tsv("../../data/morpholno.merged.rmats.tsv") %>%
  filter(geneSymbol == "NF1") %>%
  mutate(splice_type = ifelse(IncLevelDifference >0, "CLK1 Ex4 low", "CLK1 Ex4 high"),
         start = ifelse(!is.na(exonStart_0base), exonStart_0base+1,riExonStart_0base+1), 
         end = ifelse(!is.na(exonEnd), exonEnd,riExonEnd),
         is_ri = ifelse(!is.na(riExonStart_0base), TRUE, FALSE)) %>%
  select(splicing_case, start, end, IncLevelDifference, splice_type, is_ri) %>%
  left_join(nf1_exons) %>%
  mutate(riExonStart_0base = ifelse(is_ri == TRUE, start-1, NA_integer_), 
         exonStart_0base = ifelse(is_ri == FALSE, start-1, NA_integer_), 
         riExonEnd = ifelse(is_ri == TRUE, end, NA_integer_), 
         exonEnd = ifelse(is_ri == FALSE, end, NA_integer_)) %>%
  select(-start, -end, -is_ri) %>%
  dplyr::rename(exon_number_eoi = exon_number) 


morp_exon_upstream <- read_tsv("../../data/morpholno.merged.rmats.tsv") %>%
  filter(geneSymbol == "NF1") %>%
  mutate(splice_type = ifelse(IncLevelDifference >0, "CLK1 Ex4 low", "CLK1 Ex4 high"),
         start = upstreamES+1, 
         end = upstreamEE) %>%
  select(splicing_case, start, end, IncLevelDifference, splice_type) %>%
  left_join(nf1_exons) %>%
  mutate(upstreamES = start-1,
         upstreamEE = end) %>%
  select(-start, -end) %>%
  dplyr::rename(exon_number_upstream = exon_number)


morp_exon_downstream <- read_tsv("../../data/morpholno.merged.rmats.tsv") %>%
  filter(geneSymbol == "NF1") %>%
  mutate(splice_type = ifelse(IncLevelDifference >0, "CLK1 Ex4 low", "CLK1 Ex4 high"),
         start = downstreamES+1, 
         end = downstreamEE) %>%
  select(splicing_case, start, end, IncLevelDifference, splice_type) %>%
  left_join(nf1_exons) %>%
  mutate(downstreamES = start-1,
         downstreamEE = end) %>%
  select(-start, -end) %>%
  dplyr::rename(exon_number_downstream = exon_number)

morp_nf1_exons <- morp_exon_upstream %>%
  right_join(morp_exon_oi) %>%
  right_join(morp) %>%
  left_join(morp_exon_downstream) %>%
  mutate(is_sequential = (exon_number_eoi == exon_number_upstream + 1) & (exon_number_downstream == exon_number_eoi + 1)) %>%
  write_tsv("../CLK1-splicing-impact-morpholino/results/morph_nf1.tsv")


table(morp_nf1_exons$is_sequential)

  