#!/bin/sh

## create initial primary matrix
perl code/create-psi-matrix-primary.pl ../../data/splice-events-rmats.tsv.gz SE input/pan-cancer-SE-func.tsv.tmp

# Define the input and output directories
input_dir="../splicing_events_functional_sites/input/"
output_dir="input/"

# Define the bed files from uniprot
bed_files=("unipMod.hg38.col.bed" "unipOther.hg38.col.bed" "unipDisulfBond.hg38.col.bed" "unipLocSignal.hg38.col.bed" "unipDomain.hg38.col.bed")

# Define the splicing event files
#splicing_event_files=("psi-ri-primary.tsv.bed" "psi-a5ss-primary.tsv.bed" "psi-a3ss-primary.tsv.bed" "psi-se-primary.tsv.bed")
splicing_event_files=("pan-cancer-SE-func.tsv.tmp.bed")

# Loop through each splicing event file and bed file combination
for splicing_file in "${splicing_event_files[@]}"; do
  for bed_file in "${bed_files[@]}"; do
      output_file="${output_dir}${splicing_file%.*}.intersect${bed_file%.*}.wo.txt"
      bedtools intersect -wo -a "${output_dir}${splicing_file}" -b "${input_dir}${bed_file}" | sort -u > "$output_file"
  done
done

## format bedtools output
perl code/format-sites.pl

## filter primary matrices for functional sites from files outputted from above code
perl code/filter-sites.pl input/pan_cancer_SE.tmp.tsv input/psi-se.functonal-sites.tsv > input/pan-cancer-SE-func.tsv

rm *tmp*
