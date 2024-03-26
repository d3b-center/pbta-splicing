################################################################################
# 01-oncoprint.R
# written by Ammar Naqvi, Jo Lynne Rokita
#
# usage: Rscript 01-oncoprint.R 
################################################################################

## libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(maftools)
  library(vroom)
  library('data.table')
  library(ComplexHeatmap)
  library(circlize)
  
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "oncoprint")
data_dir   <- file.path(root_dir, "data")

input_dir   <- file.path(analysis_dir, "input")
plots_dir   <- file.path(analysis_dir, "plots")
results_dir   <- file.path(analysis_dir, "results")

## check and create plots/results dir
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive=TRUE)
}

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive=TRUE)
}

## input files
cons_maf_file <- file.path(data_dir,"snv-consensus-plus-hotspots.maf.tsv.gz")
tumor_only_maf_file <- file.path(data_dir,"snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz")
clin_file <- file.path(data_dir, "histologies.tsv")
indep_rna_file <- file.path(data_dir, "independent-specimens.rnaseqpanel.primary.tsv")
diff_psi_file <- file.path(data_dir, "splice-events-rmats.tsv.gz")
goi_file <- file.path(input_dir,"oncoprint-goi-lists-OpenPedCan-gencode-v39.csv")
tmb_file <- file.path(input_dir, "snv-mutation-tmb-coding.tsv")
cnv_file <- file.path(data_dir, "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz")

## color for barplot
source(file.path(input_dir, "mutation-colors.R"))

## output file
plot_out <- file.path(plots_dir,"oncoprint.pdf")

# read in files
histologies_df <- read_tsv(clin_file, guess_max = 100000)

goi <- read_csv(goi_file) %>%
  pull(HGAT) %>%
  unique() %>%
  c("CLK1")

indep_rna_df <- vroom(indep_rna_file) %>% 
  dplyr::filter(cohort == 'PBTA') %>%
  # get match id for DNA samples
  left_join(histologies_df[,c("Kids_First_Biospecimen_ID", "match_id", "short_histology", "RNA_library")]) %>%
  filter(short_histology == "HGAT",
         RNA_library %in% c("stranded", "poly-A stranded"))

matched_dna_samples <- histologies_df %>%
  filter(experimental_strategy %in% c("WGS", "WXS", "Targeted Panel"),
         is.na(RNA_library),
         !is.na(pathology_diagnosis),
         match_id %in% indep_rna_df$match_id)

tmb_df <- read_tsv(tmb_file) %>%
  filter(Tumor_Sample_Barcode %in% matched_dna_samples$Kids_First_Biospecimen_ID) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  left_join(histologies_df[c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  mutate(tmb_status = case_when(tmb <10 ~ "Normal",
                                tmb >=10 & tmb < 100 ~ "Hypermutant",
                                tmb >= 100 ~ "Ultra-hypermutant")) %>%
  # remove WXS for a sample we have WGS for
  filter(Kids_First_Biospecimen_ID != "BS_QF7M4SHH")

splice_CLK1_df <-  fread(diff_psi_file) %>%
   dplyr::filter(geneSymbol=="CLK1",
                                 exonStart_0base=="200860124", 
                                 exonEnd=="200860215") %>%
  dplyr::rename('Kids_First_Biospecimen_ID'=sample_id,
                'gene'=geneSymbol,
                PSI = IncLevel1) %>%
    left_join(histologies_df[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  dplyr::select(match_id, PSI) %>%
  filter(match_id %in% indep_rna_df$match_id)

# read in cnv file and reformat to add to maf
cnv_df <- read_tsv(cnv_file) %>%
  # select only goi, DNA samples of interest
  filter(gene_symbol %in% goi,
         biospecimen_id %in% matched_dna_samples$Kids_First_Biospecimen_ID) %>%
  mutate(Variant_Classification = case_when(status == "amplification" ~ "Amp",
                                            status == "deep deletion" ~ "Del",
                                            TRUE ~ NA_character_)) %>%
  filter(!is.na(Variant_Classification)) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = biospecimen_id,
                Hugo_Symbol = gene_symbol) %>%
  select(Kids_First_Biospecimen_ID, Hugo_Symbol, Variant_Classification)


# maf cols to select
maf_cols <- c("Hugo_Symbol", 
                  "Chromosome", 
                  "Start_Position", 
                  "End_Position", 
                  "Reference_Allele", 
                  "Tumor_Seq_Allele2", 
                  "Variant_Classification", 
                  "Variant_Type",
                  "Tumor_Sample_Barcode",
                  "t_ref_count",
                  "t_alt_count")

# read in and combine MAFs
cons_maf <- data.table::fread(cons_maf_file, data.table = FALSE) %>%
  dplyr::select(all_of(maf_cols)) 

tumor_only_maf <- data.table::fread(tumor_only_maf_file, data.table = FALSE) %>%
  dplyr::select(all_of(maf_cols)) 

maf <- cons_maf %>%
  bind_rows(tumor_only_maf) %>% 
  dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count))

## filter maf for samples with RNA splicing + HGGs
maf_filtered <- maf %>%
  dplyr::filter(Hugo_Symbol %in% goi,
                Tumor_Sample_Barcode  %in% matched_dna_samples$Kids_First_Biospecimen_ID,
                Variant_Classification %in% names(colors))
                  
collapse_snv_dat <- maf_filtered %>%
select(Tumor_Sample_Barcode,Hugo_Symbol,Variant_Classification) %>%
  dplyr::group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  bind_rows(cnv_df) %>%
  dplyr::summarise(count = as.double(length(Variant_Classification[!is.na(Variant_Classification)])),
            Variant_Classification=str_c(unique(Variant_Classification),collapse = ",")) %>%
  left_join(matched_dna_samples[,c("Kids_First_Biospecimen_ID", "match_id")]) %>%
  select(-Kids_First_Biospecimen_ID)

# get genes in order of most to least mutations
gene_row_order <- collapse_snv_dat %>%
  count(Hugo_Symbol) %>%
  arrange(-n)

# complex heatmap
gene_matrix <- reshape2::acast(collapse_snv_dat,
                               Hugo_Symbol ~ match_id,
                               value.var = "Variant_Classification",
                               fun.aggregate = function(x) paste(unique(x), collapse = ", ")) %>%
  as.data.frame() %>%
  dplyr::mutate_if(is.character, ~replace_na(.,"")) %>%
  # add multi-hits
  mutate(across(everything(), ~if_else(str_detect(., ","), "Multi_Hit", .))) %>%
  rownames_to_column(var = "Hugo_Symbol") %>%
  mutate(Sort_Order = match(Hugo_Symbol, gene_row_order$Hugo_Symbol)) %>%
  arrange(Sort_Order)

rownames(gene_matrix) <- gene_matrix$Hugo_Symbol 

gene_matrix <- gene_matrix %>%
  select(-c(Sort_Order, Hugo_Symbol)) 

# mutate the hgg dataframe for plotting
histologies_df_sorted <- splice_CLK1_df %>%
  left_join(histologies_df, by = "match_id", relationship = "many-to-many") %>%
  select(match_id, cancer_group, reported_gender, molecular_subtype, CNS_region, PSI) %>%
  unique() %>%
  group_by(match_id, cancer_group, reported_gender, molecular_subtype) %>%
  summarise(CNS_region = str_c(unique(na.omit(CNS_region)), collapse = ","),
            PSI = mean(PSI)) %>%
  left_join(unique(tmb_df[,c("tmb_status", "match_id")])) %>%
  filter(match_id %in% names(gene_matrix)) %>%
  column_to_rownames("match_id") %>%
  arrange(PSI) %>%
  dplyr::mutate(molecular_subtype = gsub(", TP53", "", molecular_subtype),
                molecular_subtype = case_when(grepl("To be classified", molecular_subtype) ~ "To be classified",
                                              TRUE ~ molecular_subtype),
                CNS_region = case_when(CNS_region == "" ~ "Unknown",
                                    TRUE ~ CNS_region),
                tmb_status = case_when(is.na(tmb_status) ~ "Unknown",
                                       TRUE ~ tmb_status)) 

# get clk1 high/low
quantiles_clk1 <- quantile(histologies_df_sorted$PSI, probs=c(.25, .75), na.rm = TRUE)
lower_sbi <- quantiles_clk1[1]
upper_sbi <- quantiles_clk1[2]

histologies_df_sorted <- histologies_df_sorted %>%
  mutate(clk1_status = case_when(PSI > upper_sbi ~ "High",
                                 PSI < lower_sbi ~ "Low",
                                 TRUE ~ "Middle")) %>%
  select(reported_gender, cancer_group, molecular_subtype, CNS_region, tmb_status, clk1_status, PSI) %>%
  dplyr::rename("Gender"=reported_gender,
                "Cancer Group" = cancer_group,
                "Molecular Subtype"=molecular_subtype,
                "CNS Region"=CNS_region, 
                "Mutation Status"=tmb_status,
                "CLK1 status" = clk1_status,
                "CLK1 Ex4 PSI"=PSI) 

# write out metadata
histologies_df_sorted %>%
  rownames_to_column(var = "match_id") %>%
  left_join(unique(histologies_df[,c("match_id", "Kids_First_Participant_ID")])) %>%
  write_tsv(file.path(results_dir, "oncoprint_sample_metadata.tsv"))

loc_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
              "#44AA99", "#882255", "#6699CC")
names(loc_cols) <- c(sort(unique(histologies_df_sorted$`CNS Region`)))

ha = HeatmapAnnotation(name = "annotation", 
                       df = histologies_df_sorted,
                       col=list(
                         "Gender" = c("Male" = "#56B4E9",
                                      "Female" = "lavender",
                                      "Not Reported" = "whitesmoke"),
                         "Cancer Group" = c("Diffuse hemispheric glioma" = "springgreen4", 
                                            "Diffuse midline glioma" = "#ff40d9",
                                            "High-grade glioma" = "#ffccf5",
                                            "Infant-type hemispheric glioma" = "lightblue2",
                                            "Pleomorphic xanthoastrocytoma" = "navy"),
                         "Molecular Subtype" = c("DHG, H3 G35" = "springgreen4",
                                                 "DMG, H3 K28" = "#ff40d9",
                                                 "HGG, H3 wildtype" = "lightpink",
                                                 "HGG, IDH" = "indianred",
                                                 "IHG, NTRK-altered" = "cornflowerblue",
                                                 "IHG, ALK-altered" = "skyblue4",
                                                 "IHG, ROS1-altered" = "lightblue1",
                                                 "HGG, PXA" = "navy"),
                         "CNS Region" = loc_cols,
                         "Mutation Status" = c("Normal" = "grey80",
                                               "Hypermutant" = "orange",
                                               "Ultra-hypermutant" = "red",
                                               "Unknown" = "whitesmoke"),
                         "CLK1 status" = c("High" = "navy", "Low" = "#CAE1FF",  "Middle" = "grey80"),
                         "CLK1 Ex4 PSI" = colorRamp2(c(0, 0.25, 0.75, 1.0), c("whitesmoke", "#CAE1FF","cornflowerblue", "navy")),
                         annotation_name_side = "right", 
                         annotation_name_gp = gpar(fontsize = 9),
                         na_col = "whitesmoke"))
                       

col = colors
df = histologies_df_sorted

gene_matrix_sorted <- gene_matrix %>%
  select(all_of(rownames(df)))

# global option to increase space between heatmap and annotations
ht_opt$ROW_ANNO_PADDING = unit(1.25, "cm")

plot_oncoprint <- oncoPrint(gene_matrix_sorted[1:50,], get_type = function(x) strsplit(x, ",")[[1]],
          column_names_gp = gpar(fontsize = 9), show_column_names = F,
          alter_fun = list(
            background = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "whitesmoke",col="whitesmoke")),
            Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Missense_Mutation"]),col = NA)),
            Nonsense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Nonsense_Mutation"]),col = NA)),
            Frame_Shift_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Frame_Shift_Del"]), col = NA)),
            Frame_Shift_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Frame_Shift_Ins"]), col = NA)),
            Splice_Site = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Splice_Site"]), col = NA)),
            Translation_Start_Site = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Translation_Start_Site"]), col = NA)),
            Nonstop_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Nonstop_Mutation"]),col = NA)),
            In_Frame_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["In_Frame_Del"]),col = NA)),
            In_Frame_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["In_Frame_Ins"]), col = NA)),
            Stop_Codon_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Stop_Codon_Ins"]), col = NA)),
            Start_Codon_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Start_Codon_Del"]), col = NA)),
            Fusion = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Fusion"]),col = NA)),
            Multi_Hit_Fusion  = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Multi_Hit_Fusion"]),col = NA)),
            Multi_Hit = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Multi_Hit"]), col = NA)),
            Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Del"]), col = NA)),
            Amp = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Amp"]), col = NA)),
            Loss = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Loss"]), col = NA))),
          col = col,
          top_annotation = ha,
          alter_fun_is_vectorized = TRUE,
          #bottom_annotation = ha1,
          column_order =  colnames(gene_matrix_sorted))
                  
# Save plot as PDF
pdf(plot_out, width = 15, height = 8)
plot_oncoprint
dev.off()

# create df for enrichment
ids_clk1 <- histologies_df_sorted %>%
  rownames_to_column(var = "match_id") %>%
  select(match_id, `CLK1 status`)

total_high <- nrow(ids_clk1)/2
total_low <- nrow(ids_clk1)/2


alteration_counts <- collapse_snv_dat %>%
  full_join(ids_clk1) %>%
  filter(`CLK1 status` != "Middle") %>%
## group by junction and calculate means
select(Hugo_Symbol, `CLK1 status`) %>%
  group_by(Hugo_Symbol, `CLK1 status`) %>%
  count() %>%
  ungroup() %>%
  # Spread to wide format to get separate columns for "High" and "Low"
  pivot_wider(names_from = `CLK1 status`, values_from = n, values_fill = list(n = 0)) %>%
  rowwise() %>%
   mutate(
      Fisher_Test = list(
        fisher.test(
          matrix(
            c(High, total_high - High,  # Counts of High and the absence of High
              Low, total_low - Low),    # Counts of Low and the absence of Low
            nrow = 2
          )
        )
      ),
      P_Value = Fisher_Test$p.value
    ) %>%
    select(Hugo_Symbol, High, Low, P_Value) %>%
    ungroup() %>%
  arrange(P_Value) %>%
  write_tsv(file.path(results_dir, "clk1_high_low_mutation_counts.tsv"))


