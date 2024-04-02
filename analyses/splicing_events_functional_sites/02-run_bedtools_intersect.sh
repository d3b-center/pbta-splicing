  #!/bin/sh

## Cross reference splicing events with Unipro database (*hg38.col.bed)
#!/bin/sh

# Define the input and output directories
input_dir="input/"
output_dir="results/"

# Define the bed files
bed_files=("unipMod.hg38.col.bed" "unipOther.hg38.col.bed" "unipDisulfBond.hg38col.bed" "unipLocSignal.hg38.col.bed" "unipDomain.hg38.col.bed")

# Define the splicing event files
splicing_event_files=("splicing_events.SE.total.Other high-grade glioma.pos.bed" "splicing_events.SE.total.Other high-grade glioma.neg.bed" "splicing_events.SE.total.DMG.pos.bed" "splicing_events.SE.total.DMG.neg.bed" "splicing_events.SE.total.pos.bed" "splicing_events.SE.total.neg.bed")

# Loop through each splicing event file and bed file combination
for splicing_file in "${splicing_event_files[@]}"; do
    for bed_file in "${bed_files[@]}"; do
        output_file="${output_dir}${splicing_file%.*}.intersect${bed_file%.*}.wo.txt"
        bedtools intersect -wo -a "${output_dir}${splicing_file}" -b "${input_dir}${bed_file}" | sort -u > "$output_file"
    done
done
