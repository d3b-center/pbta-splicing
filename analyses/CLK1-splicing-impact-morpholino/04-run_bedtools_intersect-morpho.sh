#!/bin/sh

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

rm results/*.wo.txt
rm input/*bed
