#!/bin/sh

cat ../../data/morpholno.merged.rmats.tsv | grep ^SE |  awk -F "\t" '{if( ($36 < -.10) || ($36 >.10) ){ print $0}}' | awk -F "\t" '{if($33< 0.05){ print $0}}' |awk -F "\t" '{print $6"\t"$8"\t"$9"\t"$5":"$8"-"$9"\t"$36"\t"$7}' | perl -pe 's/\"//g' > input/morpho.diff.SE.bed
cat ../../data/morpholno.merged.rmats.tsv | grep ^RI |  awk -F "\t" '{if( ($36 < -.10) || ($36 >.10) ){ print $0}}' | awk -F "\t" '{if($33< 0.05){ print $0}}' |awk -F "\t" '{print $6"\t"$14"\t"$15"\t"$5":"$14"-"$15"\t"$36"\t"$7}' | perl -pe 's/\"//g' > input/morpho.diff.RI.bed
cat ../../data/morpholno.merged.rmats.tsv | grep ^A5SS |  awk -F "\t" '{if( ($36 < -.10) || ($36 >.10) ){ print $0}}' | awk -F "\t" '{if($33< 0.05){ print $0}}' |awk -F "\t" '{print $6"\t"$20"\t"$21"\t"$5":"$20"-"$21"\t"$36"\t"$7}' | perl -pe 's/\"//g' > input/morpho.diff.A5SS.bed
cat ../../data/morpholno.merged.rmats.tsv | grep ^A3SS |  awk -F "\t" '{if( ($36 < -.10) || ($36 >.10) ){ print $0}}' | awk -F "\t" '{if($33< 0.05){ print $0}}' |awk -F "\t" '{print $6"\t"$20"\t"$21"\t"$5":"$20"-"$21"\t"$36"\t"$7}' | perl -pe 's/\"//g' > input/morpho.diff.A3SS.bed

bedtools intersect -wo -a input/morpho.diff.SE.bed -b ../splicing_events_functional_sites/input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.morpho.SE.intersectUnipMod.wo.txt
bedtools intersect -wo -a input/morpho.diff.SE.bed -b ../splicing_events_functional_sites/input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.morpho.SE.intersectUnipOther.wo.txt
bedtools intersect -wo -a input/morpho.diff.SE.bed -b ../splicing_events_functional_sites/input/unipDisulfBond.hg38.col.bed |sort -u > results/splicing_events.morpho.SE.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a input/morpho.diff.SE.bed -b ../splicing_events_functional_sites/input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.morpho.SE.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a input/morpho.diff.SE.bed -b ../splicing_events_functional_sites/input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.morpho.SE.intersectUnipDomain.wo.txt

bedtools intersect -wo -a input/morpho.diff.RI.bed -b ../splicing_events_functional_sites/input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.morpho.RI.intersectUnipMod.wo.txt
bedtools intersect -wo -a input/morpho.diff.RI.bed -b ../splicing_events_functional_sites/input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.morpho.RI.intersectUnipOther.wo.txt
bedtools intersect -wo -a input/morpho.diff.RI.bed -b ../splicing_events_functional_sites/input/unipDisulfBond.hg38.col.bed |sort -u > results/splicing_events.morpho.RI.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a input/morpho.diff.RI.bed -b ../splicing_events_functional_sites/input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.morpho.RI.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a input/morpho.diff.RI.bed -b ../splicing_events_functional_sites/input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.morpho.RI.intersectUnipDomain.wo.txt

bedtools intersect -wo -a input/morpho.diff.A5SS.bed -b ../splicing_events_functional_sites/input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.morpho.A5SS.intersectUnipMod.wo.txt
bedtools intersect -wo -a input/morpho.diff.A5SS.bed -b ../splicing_events_functional_sites/input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.morpho.A5SS.intersectUnipOther.wo.txt
bedtools intersect -wo -a input/morpho.diff.A5SS.bed -b ../splicing_events_functional_sites/input/unipDisulfBond.hg38.col.bed |sort -u > results/splicing_events.morpho.A5SS.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a input/morpho.diff.A5SS.bed -b ../splicing_events_functional_sites/input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.morpho.A5SS.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a input/morpho.diff.A5SS.bed -b ../splicing_events_functional_sites/input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.morpho.A5SS.intersectUnipDomain.wo.txt

bedtools intersect -wo -a input/morpho.diff.A3SS.bed -b ../splicing_events_functional_sites/input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.morpho.A3SS.intersectUnipMod.wo.txt
bedtools intersect -wo -a input/morpho.diff.A3SS.bed -b ../splicing_events_functional_sites/input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.morpho.A3SS.intersectUnipOther.wo.txt
bedtools intersect -wo -a input/morpho.diff.A3SS.bed -b ../splicing_events_functional_sites/input/unipDisulfBond.hg38.col.bed |sort -u > results/splicing_events.morpho.A3SS.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a input/morpho.diff.A3SS.bed -b ../splicing_events_functional_sites/input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.morpho.A3SS.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a input/morpho.diff.A3SS.bed -b ../splicing_events_functional_sites/input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.morpho.A3SS.intersectUnipDomain.wo.txt



## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.morpho.SE.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.SE.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.morpho.SE.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.SE.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.morpho.SE.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.SE.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.morpho.SE.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.SE.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.morpho.SE.intersectUnip.ggplot.txt

## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.morpho.RI.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.RI.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.morpho.RI.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.RI.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.morpho.RI.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.RI.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.morpho.RI.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.RI.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.morpho.RI.intersectUnip.ggplot.txt

## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.morpho.A5SS.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.A5SS.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.morpho.A5SS.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.A5SS.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.morpho.A5SS.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.A5SS.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.morpho.A5SS.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.A5SS.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.morpho.A5SS.intersectUnip.ggplot.txt

## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.morpho.A3SS.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.A3SS.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.morpho.A3SS.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.A3SS.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.morpho.A3SS.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.A3SS.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.morpho.A3SS.intersectUnip.ggplot.txt
cat results/splicing_events.morpho.A3SS.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.morpho.A3SS.intersectUnip.ggplot.txt




rm results/*.wo.txt
rm input/*bed
