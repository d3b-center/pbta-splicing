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

## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.other-HGG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.other-HGG.pos.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.other-HGG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.other-HGG.pos.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.other-HGG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.other-HGG.pos.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.other-HGG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.other-HGG.pos.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.other-HGG.pos.intersectUnip.ggplot.txt

echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.other-HGG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.other-HGG.neg.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.other-HGG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.other-HGG.neg.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.other-HGG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.other-HGG.neg.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.other-HGG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.other-HGG.neg.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.other-HGG.neg.intersectUnip.ggplot.txt

## Generate for ggplot for all events corresponding to functional sites
echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.DMG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.DMG.pos.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.DMG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.DMG.pos.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.DMG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.DMG.pos.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.DMG.pos.intersectUnip.ggplot.txt
cat results/splicing_events.total.DMG.pos.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.DMG.pos.intersectUnip.ggplot.txt

echo -e "SpliceID\tdPSI\tUniprot" > results/splicing_events.total.DMG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.DMG.neg.intersectUnipMod.wo.txt |  awk '{print $4"\t"$5"\tModifications"}'|sort -u >> results/splicing_events.total.DMG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.DMG.neg.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' |sort -u >> results/splicing_events.total.DMG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.DMG.neg.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' |sort -u  >> results/splicing_events.total.DMG.neg.intersectUnip.ggplot.txt
cat results/splicing_events.total.DMG.neg.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' |sort -u  >> results/splicing_events.total.DMG.neg.intersectUnip.ggplot.txt
