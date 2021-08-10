#!/bin/sh

echo "process rMATS with .10 dPSI and 10 junction read counts...\n";

## Process rMATS files given histologies file. Keep only DMG samples and storng splicing events
#perl extract_recurrent_splicing_events.pl ../../data/pbta-histologies.tsv

echo "bedtools intersect...\n";

## Cross reference splicing events with Unipro database (*hg38.col.bed)
bedtools intersect -wo -a results/splicing_events.total.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.total.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.total.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.total.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.total.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.intersectUnipDomain.wo.txt

echo "make tab for ggplot ...\n";

## Generate for ggplot for all events corresponding to functional sites
cat results/splicing_events.total.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u > results/splicing_events.total.intersectUnip.ggplot.txt
cat results/splicing_events.total.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.intersectUnip.ggplot.txt
cat results/splicing_events.total.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.intersectUnip.ggplot.txt
cat results/splicing_events.total.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.intersectUnip.ggplot.txt

## Divide into postivie and negative dPSIs or decreased/increased inclusions
echo "dividing into positive and negative..."
cat results/splicing_events.total.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u | awk '{if($2>0){ print $0}'} > results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u | awk '{if($2>0){ print $0}'} >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' | awk '{if($2>0){ print $0}'}| sort -u  >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' | awk '{if($2>0){ print $0}'}| sort -u  >> results/splicing_events.total.pos.intersectUnip.ggplot.txt

cat results/splicing_events.total.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u | awk '{if($2<0){ print $0}'} > results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u | awk '{if($2<0){ print $0}'} >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' | awk '{if($2<0){ print $0}'}| sort -u  >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' | awk '{if($2<0){ print $0}'}| sort -u  >> results/splicing_events.total.neg.intersectUnip.ggplot.txt

## Remove intermediatery files / off for now
rm results/splicing_events.total.*intersectUnipMod.wo.txt

# cat output/dominant_events_lsvs.total.intersectUnipDomain.wo.txt | awk '{print $4"\t"$5"\tDomain"}' >> output/dominant_events_lsvs.total.intersectUnip.ggplot.txt
