#!/bin/sh

<<<<<<< HEAD
cat ../../data/morpholno.merged.rmats.tsv | grep ^SE |  awk -F "\t" '{if( ($36 < -.10) || ($36 >.10) ){ print $0}}' | awk -F "\t" '{if($33< 0.05){ print $0}}' |awk -F "\t" '{print $6"\t"$8"\t"$9"\t"$5":"$8"-"$9"\t"$36"\t"$7}' | perl -pe 's/\"//g' > input/morpho.diff.bed
bedtools intersect -wo -a input/morpho.diff.bed -b ../splicing_events_functional_sites/input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.morpho.intersectUnipMod.wo.txt
bedtools intersect -wo -a input/morpho.diff.bed -b ../splicing_events_functional_sites/input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.morpho.intersectUnipOther.wo.txt
bedtools intersect -wo -a input/morpho.diff.bed -b ../splicing_events_functional_sites/input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.morpho.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a input/morpho.diff.bed -b ../splicing_events_functional_sites/input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.morpho.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a input/morpho.diff.bed -b ../splicing_events_functional_sites/input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.morpho.intersectUnipDomain.wo.txt

## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.morpho.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.morpho.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.morpho.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.morpho.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.morpho.intersectUnip.ggplot.txt
=======
events=("SE" "RI" "A5SS" "A3SS")

for event in "${events[@]}"; do
    # Determine the column numbers based on the event type
    if [ "$event" = "SE" ]; then
        col1=6
        col2=8
        col3=9
    elif [ "$event" = "RI" ]; then
        col1=6
        col2=14
        col3=15
    else
        col1=6
        col2=20
        col3=21
    fi

    # Generate the bed file
    cat ../../data/morpholno.merged.rmats.tsv | grep "^$event" | awk -F "\t" '{if( ($36 < -.10) || ($36 >.10) ){ print $0}}' | awk -F "\t" '{if($33< 0.05){ print $0}}' |awk -F "\t" -v col1="$col1" -v col2="$col2" -v col3="$col3" '{print $col1"\t"$col2"\t"$col3"\t"$5":"$col2"-"$col3"\t"$36"\t"$7}' | perl -pe 's/\"//g' > input/morpho.diff.$event.bed

    # Intersect with functional site beds
    for functional_site in "Mod" "Other" "DisulfBond" "LocSignal" "Domain"; do
        bedtools intersect -wo -a input/morpho.diff.$event.bed -b ../splicing_events_functional_sites/input/unip$functional_site.hg38.col.bed | sort -u > results/splicing_events.morpho.$event.intersectUnip$functional_site.wo.txt
    done

    # Generate ggplot data
    echo -e "SpliceID\tdPSI\tUniprot\tType" > results/splicing_events.morpho.$event.intersectUnip.ggplot.txt
    for functional_site in "Mod" "Other" "DisulfBond" "LocSignal"; do
        cat results/splicing_events.morpho.$event.intersectUnip$functional_site.wo.txt | awk '{print $4"\t"$5"\t'$functional_site'\t'$event'"}' | sort -u >> results/splicing_events.morpho.$event.intersectUnip.ggplot.txt
    done
done
>>>>>>> main

rm results/*.wo.txt
rm input/*bed
