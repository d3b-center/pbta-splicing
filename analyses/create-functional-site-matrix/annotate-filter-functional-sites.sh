#!/bin/sh

# Define the directories
data_dir="../../data"
results_dir="results/"
bed_dir="../splicing_events_functional_sites/input/"

## create initial primary matrix
perl 01-create-psi-matrix-primary.pl $data_dir/splice-events-rmats.tsv.gz SE $results_dir/pan-cancer-SE-func.tsv.tmp

# Define the bed files from uniprot
bed_files=("unipMod.hg38.col.bed" "unipOther.hg38.col.bed" "unipDisulfBond.hg38.col.bed" "unipLocSignal.hg38.col.bed" "unipDomain.hg38.col.bed")

# Define the splicing event files
splicing_event_files=("pan-cancer-SE-func.tsv.tmp.bed")

# Loop through each splicing event file and bed file combination
for splicing_file in "${splicing_event_files[@]}"; do
  for bed_file in "${bed_files[@]}"; do
      output_file="${results_dir}${splicing_file%.*}.intersect${bed_file%.*}.wo.txt"
      bedtools intersect -wo -a "${results_dir}${splicing_file}" -b "${bed_dir}${bed_file}" | sort -u > "$output_file"
  done
done

## format bedtools output
perl 02-format-sites.pl "${results_dir}pan-cancer-SE-func.tsv.tmp.intersectunipLocSignal.hg38.col.wo.txt" "${results_dir}pan-cancer-SE-func.tsv.tmp.intersectunipDisulfBond.hg38.col.wo.txt" "${results_dir}pan-cancer-SE-func.tsv.tmp.intersectunipMod.hg38.col.wo.txt" "${results_dir}pan-cancer-SE-func.tsv.tmp.intersectunipOther.hg38.col.wo.txt" "${results_dir}pan-cancer-SE-func.list.tmp"

## filter primary matrices for functional sites from files outputted from above code
perl 03-filter-sites.pl "${results_dir}pan-cancer-SE-func.tsv.tmp" "${results_dir}pan-cancer-SE-func.list.tmp" > "${results_dir}pan-cancer-SE-func.tsv"
perl 03-filter-sites.pl "${results_dir}pan-cancer-SE-func.tsv.tmp" "${results_dir}pan-cancer-SE-func.list.tmp" > "${results_dir}pan-cancer-SE-func.tsv"
gzip ${results_dir}pan-cancer-SE-func.tsv

rm $results_dir/*tmp
rm $results_dir/*tmp*
