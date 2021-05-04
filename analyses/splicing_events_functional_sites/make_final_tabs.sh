#!/bin/sh

echo "process rMATS with .10 dPSI and 10 junction read counts...\n";

## Process rMATS files given histologies file. Keep only DMG samples and storng splicing events
perl extract_recurrent_splicing_events.pl ../../data/pbta-histologies.tsv

echo "bedtools intersect...\n";

## Cross reference splicing events with Unipro database (*hg38.col.bed)
bedtools intersect -wo -a results/splicing_events.dpsi10.jc10.rec2.bed -b unipMod.hg38.col.bed       |sort -u > results/dominant_events_lsvs.total.rec2.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.dpsi10.jc10.rec2.bed -b unipOther.hg38.col.bed     |sort -u > results/dominant_events_lsvs.total.rec2.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.dpsi10.jc10.rec2.bed -b unipDisulfBond.hg38col.bed |sort -u > results/dominant_events_lsvs.total.rec2.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.dpsi10.jc10.rec2.bed -b unipLocSignal.hg38.col.bed |sort -u > results/dominant_events_lsvs.total.rec2.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.dpsi10.jc10.rec2.bed -b unipDomain.hg38.col.bed    |sort -u > results/dominant_events_lsvs.total.rec2.intersectUnipDomain.wo.txt

echo "make tab for ggplot ...\n";

## Generate for ggplot for all events corresponding to functional sites
cat results/dominant_events_lsvs.total.rec2.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u > results/dominant_events_lsvs.total.rec2.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/dominant_events_lsvs.total.rec2.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/dominant_events_lsvs.total.rec2.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/dominant_events_lsvs.total.rec2.intersectUnip.ggplot.txt

## Divide into postivie and negative dPSIs or decreased/increased inclusions
echo "dividing into positive and negative..."
cat results/dominant_events_lsvs.total.rec2.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u | awk '{if($2>0){ print $0}'} > results/dominant_events_lsvs.total.rec2.pos.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u | awk '{if($2>0){ print $0}'} >> results/dominant_events_lsvs.total.rec2.pos.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' | awk '{if($2>0){ print $0}'}| sort -u  >> results/dominant_events_lsvs.total.rec2.pos.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' | awk '{if($2>0){ print $0}'}| sort -u  >> results/dominant_events_lsvs.total.rec2.pos.intersectUnip.ggplot.txt

cat results/dominant_events_lsvs.total.rec2.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u | awk '{if($2<0){ print $0}'} > results/dominant_events_lsvs.total.rec2.neg.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u | awk '{if($2<0){ print $0}'} >> results/dominant_events_lsvs.total.rec2.neg.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' | awk '{if($2<0){ print $0}'}| sort -u  >> results/dominant_events_lsvs.total.rec2.neg.intersectUnip.ggplot.txt
cat results/dominant_events_lsvs.total.rec2.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' | awk '{if($2<0){ print $0}'}| sort -u  >> results/dominant_events_lsvs.total.rec2.neg.intersectUnip.ggplot.txt

## Remove intermediatery files / off for now
#rm results/dominant_events_lsvs.total.rec2.*intersectUnipMod.wo.txt

# cat output/dominant_events_lsvs.total.intersectUnipDomain.wo.txt | awk '{print $4"\t"$5"\tDomain"}' >> output/dominant_events_lsvs.total.intersectUnip.ggplot.txt
