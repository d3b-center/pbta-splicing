#!/bin/sh

# Define input and output directories
input_dir="input/"
output_dir="results/"

# Define the bed files
bed_files=("unipMod.hg38.col.bed" "unipOther.hg38.col.bed" "unipDisulfBond.hg38col.bed" "unipLocSignal.hg38.col.bed" "unipDomain.hg38.col.bed")
#!/bin/sh

# Define input and output directories
input_dir="input/"
output_dir="results/"

# Define the bed files
bed_files=("unipMod.hg38.col.bed" "unipOther.hg38.col.bed" "unipDisulfBond.hg38col.bed" "unipLocSignal.hg38.col.bed" "unipDomain.hg38.col.bed")

# Define corresponding functional site names
functional_sites=("Modifications" "Other" "DisulfBond" "LocSignal" "Domain")

# Define splicing event files
splicing_event_files=("splicing_events.SE.total.pos" "splicing_events.SE.total.neg" "splicing_events.SE.total.Other high-grade glioma.pos" "splicing_events.SE.total.Other high-grade glioma.neg" "splicing_events.SE.total.DMG.pos" "splicing_events.SE.total.DMG.neg")

# Loop through each splicing event file
for splicing_file in "${splicing_event_files[@]}"; do
    # Create a temporary file for each splicing event file
    tmp_file="${output_dir}${splicing_file}.intersectUnip.ggplot.tmp.txt"
    echo -e "SpliceID\tdPSI\tUniprot" > "$tmp_file"

    # Loop through each bed file
    for ((i=0; i<${#bed_files[@]}; i++)); do
        # Generate ggplot data and append to the temporary file
        cat "${output_dir}${splicing_file}.intersect${bed_files[i]%.*}.wo.txt" | awk -v var="${bed_files[i]}" -v site="${functional_sites[i]}" '{print $4"\t"$5"\t"site}' | sort -u >> "$tmp_file"
    done

    # Create the final ggplot file
    final_file="${output_dir}${splicing_file}.intersectUnip.ggplot.txt"
    cat "$tmp_file" > "$final_file"

    # Remove the temporary file
    rm "$tmp_file"
done
