#!/bin/sh

echo "process rMATS with .20 dPSI and 10 junction read counts...";

## Process rMATS files given histologies file. Keep only HGG midlines samples and storng splicing events
perl extract_recurrent_splicing_events_hgg.pl ../../data/v19_plus_20210311_pnoc_rna.tsv

echo "bedtools intersect...";

## Cross reference splicing events with Unipro database (*hg38.col.bed)
bedtools intersect -wo -a results/splicing_events.total.pos.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.pos.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.total.pos.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.pos.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.total.pos.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.pos.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.total.pos.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.pos.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.total.pos.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.pos.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.total.neg.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.neg.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.total.neg.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.neg.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.total.neg.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.neg.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.total.neg.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.neg.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.total.neg.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.neg.intersectUnipDomain.wo.txt



echo "make tab for ggplot ...";

## Generate for ggplot for all events corresponding to functional sites
echo "Splice_ID\tdPSI\tUniprot" > results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.pos.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.pos.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.pos.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.pos.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.pos.intersectUnip.ggplot.txt

echo "Splice_ID\tdPSI\tUniprot" > results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.neg.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.neg.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.neg.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.neg.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.neg.intersectUnip.ggplot.txt



# ## Divide into postivie and negative dPSIs or decreased/increased inclusions
# echo "dividing into positive and negative..."
#
# echo "Splice_ID\tdPSI\tUniprot" > results/splicing_events.total.pos.intersectUnip.ggplot.txt
# cat results/splicing_events.total.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u | awk '{if($2>0){ print $0}'} >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
# cat results/splicing_events.total.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u | awk '{if($2>0){ print $0}'} >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
# cat results/splicing_events.total.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' | awk '{if($2>0){ print $0}'}| sort -u  >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
# cat results/splicing_events.total.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' | awk '{if($2>0){ print $0}'}| sort -u  >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
#
# echo "Splice_ID\tdPSI\tUniprot" > results/splicing_events.total.neg.intersectUnip.ggplot.txt
# cat results/splicing_events.total.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u | awk '{if($2<0){ print $0}'} >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
# cat results/splicing_events.total.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u | awk '{if($2<0){ print $0}'} >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
# cat results/splicing_events.total.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' | awk '{if($2<0){ print $0}'}| sort -u  >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
# cat results/splicing_events.total.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' | awk '{if($2<0){ print $0}'}| sort -u  >> results/splicing_events.total.neg.intersectUnip.ggplot.txt

## Remove intermediatery files / off for now
rm results/splicing_events.total.*intersectUnipMod.wo.txt

# cat output/dominant_events_lsvs.total.intersectUnipDomain.wo.txt | awk '{print $4"\t"$5"\tDomain"}' >> output/dominant_events_lsvs.total.intersectUnip.ggplot.txt
