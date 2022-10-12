#!/bin/sh

## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.pos.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.pos.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.pos.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.pos.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.pos.intersectUnip.ggplot.txt

echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.neg.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.neg.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.neg.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.neg.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.neg.intersectUnip.ggplot.txt
