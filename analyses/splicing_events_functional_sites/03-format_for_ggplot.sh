#!/bin/sh

## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.SE.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.SE.total.pos.intersectUnipMod.hg38.col.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.SE.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.SE.total.pos.intersectUnipOther.hg38.col.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.SE.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.SE.total.pos.intersectUnipDisulfBond.hg38.col.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.SE.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.SE.total.pos.intersectUnipLocSignal.hg38.col.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.SE.total.pos.intersectUnip.ggplot.txt

echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.SE.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.SE.total.neg.intersectUnipMod.hg38.col.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.SE.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.SE.total.neg.intersectUnipOther.hg38.col.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.SE.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.SE.total.neg.intersectUnipDisulfBond.hg38.col.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.SE.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.SE.total.neg.intersectUnipLocSignal.hg38.col.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.SE.total.neg.intersectUnip.ggplot.txt
