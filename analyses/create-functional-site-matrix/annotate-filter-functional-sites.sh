#!/bin/sh

# Define the directories
data_dir="../../data"
results_dir="../create-functional-site-matrix/results"
bed_dir="../splicing_events_functional_sites/input"

## create initial primary matrix
perl 01-create-psi-matrix-primary.pl $data_dir/splice-events-rmats.tsv.gz SE $results_dir/pan-cancer-SE-func.tsv.tmp

# Define the bed files from uniprot
bed_files=("$bed_dir/unipMod.hg38.col.bed" "$bed_dir/unipOther.hg38.col.bed" "$bed_dir/unipDisulfBond.hg38.col.bed" "$bed_dir/unipLocSignal.hg38.col.bed" "$bed_dir/unipDomain.hg38.col.bed")

# Define the splicing event files
splicing_event_files=("pan-cancer-SE-func.tsv.tmp.bed")

# Loop through each splicing event file and bed file combination
for splicing_file in "${splicing_event_files[@]}"; do
  for bed_file in "${bed_files[@]}"; do
      output_file="${output_dir}${splicing_file%.*}.intersect${bed_file%.*}.wo.txt"
      bedtools intersect -wo -a "${output_dir}${splicing_file}" -b "${input_dir}${bed_file}" | sort -u > "$output_file"
  done
done

## format bedtools output
perl 02-format-sites.pl

## filter primary matrices for functional sites from files outputted from above code
perl 03-filter-sites.pl $results_dir/pan-cancer-SE-func.tsv.tmp $results_dir/pan-cancer-SE-func.list.tmp > $results_dir/pan-cancer-SE-func.tsv

#rm $results_dir/*tmp
rm $results_dir/*tmp*
