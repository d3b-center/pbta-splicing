#!/bin/sh

## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.HGG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.HGG.pos.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.HGG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.HGG.pos.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.HGG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.HGG.pos.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.HGG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.HGG.pos.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.HGG.pos.intersectUnip.ggplot.txt

echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.HGG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.HGG.neg.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.HGG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.HGG.neg.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.HGG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.HGG.neg.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.HGG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.HGG.neg.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.HGG.neg.intersectUnip.ggplot.txt

## H3.3
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.H33.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H33.pos.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.H33.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H33.pos.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.H33.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H33.pos.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.H33.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H33.pos.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.H33.pos.intersectUnip.ggplot.txt

echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.H33.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H33.neg.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.H33.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H33.neg.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.H33.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H33.neg.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.H33.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H33.neg.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.H33.neg.intersectUnip.ggplot.txt

### H3.1
## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.H31.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H31.pos.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.H31.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H31.pos.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.H31.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H31.pos.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.H31.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H31.pos.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.H31.pos.intersectUnip.ggplot.txt

echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.H31.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H31.neg.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.H31.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H31.neg.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.H31.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H31.neg.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.H31.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H31.neg.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.H31.neg.intersectUnip.ggplot.txt

##H3 WT
## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.H3WT.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H3WT.pos.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.H3WT.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H3WT.pos.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.H3WT.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H3WT.pos.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.H3WT.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.H3WT.pos.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.H3WT.pos.intersectUnip.ggplot.txt

echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.H3WT.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H3WT.neg.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.H3WT.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H3WT.neg.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.H3WT.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H3WT.neg.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.H3WT.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.H3WT.neg.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.H3WT.neg.intersectUnip.ggplot.txt
