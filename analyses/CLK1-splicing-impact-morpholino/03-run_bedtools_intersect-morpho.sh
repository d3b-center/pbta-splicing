#!/bin/sh

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
    cat ../../data/morpholno.merged.rmats.tsv | grep "^$event" | awk -F "\t" '{if( ($36 < -.10) || ($36 >.10) ){ print $0}}' | awk -F "\t" '{if($32 < 0.05 && $33< 0.05){ print $0}}' |awk -F "\t" -v col1="$col1" -v col2="$col2" -v col3="$col3" '{print $col1"\t"$col2"\t"$col3"\t"$5":"$col2"-"$col3"\t"$36"\t"$7}' | perl -pe 's/\"//g' > input/morpho.diff.$event.bed

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

rm results/*.wo.txt
rm input/*bed
