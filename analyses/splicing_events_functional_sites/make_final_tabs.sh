#!/bin/sh

#perl scripts/extract_recurrent_dominant_lsvs.pl | grep "*" | awk '{print $6"\t"$1"\t"$2"\t"$7}'| perl -pe 's/(chr[\d+XY]+)\:/$1\t/g' | awk -F "-" '{print $1"\t"$2"-"$3"-"$4}' | awk '{if($6 == "+"){ print $0 } else {print $0"-" }}' | perl -pe 's/\*//' > output/dominant_changes_lsv.bed
#perl scripts/extract_recurrent_dominant_lsvs.pl | grep "*" | perl -pe 's/\*//'  > output/dominant_changes_lsv.txt
#perl scripts/extract_recurrent_dominant_lsvs.pl | grep "*" | perl scripts/reformat_for_bed.pl >  output/dominant_changes_lsv.bed

echo "bedtools intersect...\n";

## Unipro cross reference
bedtools intersect -wo -a results/splicing_events.dpsi10.jc10.rec2.bed -b unipMod.hg38.col.bed       |sort -u > results/dominant_events_lsvs.total.rec2.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.dpsi10.jc10.rec2.bed -b unipOther.hg38.col.bed     |sort -u > results/dominant_events_lsvs.total.rec2.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.dpsi10.jc10.rec2.bed -b unipDisulfBond.hg38col.bed |sort -u > results/dominant_events_lsvs.total.rec2.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.dpsi10.jc10.rec2.bed -b unipLocSignal.hg38.col.bed |sort -u > results/dominant_events_lsvs.total.rec2.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.dpsi10.jc10.rec2.bed -b unipDomain.hg38.col.bed    |sort -u > results/dominant_events_lsvs.total.rec2.intersectUnipDomain.wo.txt

echo "make tab for ggplot ...\n";
## generate for ggplot
cat results/dominant_events_lsvs.total.rec2.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u > results/dominant_events_lsvs.total.rec2.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/dominant_events_lsvs.total.rec2.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/dominant_events_lsvs.total.rec2.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/dominant_events_lsvs.total.rec2.intersectUnip.ggplot.txt

echo "dividing into positive and negative..."
cat results/dominant_events_lsvs.total.rec2.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u | awk '{if($2>0){ print $0}'} > results/dominant_events_lsvs.total.rec2.pos.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u | awk '{if($2>0){ print $0}'} >> results/dominant_events_lsvs.total.rec2.pos.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' | awk '{if($2>0){ print $0}'}| sort -u  >> results/dominant_events_lsvs.total.rec2.pos.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' | awk '{if($2>0){ print $0}'}| sort -u  >> results/dominant_events_lsvs.total.rec2.pos.intersectUnip.ggplot.txt

cat results/dominant_events_lsvs.total.rec2.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u | awk '{if($2<0){ print $0}'} > results/dominant_events_lsvs.total.rec2.neg.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u | awk '{if($2<0){ print $0}'} >> results/dominant_events_lsvs.total.rec2.neg.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' | awk '{if($2<0){ print $0}'}| sort -u  >> results/dominant_events_lsvs.total.rec2.neg.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' | awk '{if($2<0){ print $0}'}| sort -u  >> results/dominant_events_lsvs.total.rec2.neg.intersectUnip.ggplot.txt


# cat output/dominant_events_lsvs.total.intersectUnipDomain.wo.txt | awk '{print $4"\t"$5"\tDomain"}' >> output/dominant_events_lsvs.total.intersectUnip.ggplot.txt
