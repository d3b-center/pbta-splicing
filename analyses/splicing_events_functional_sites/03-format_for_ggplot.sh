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
