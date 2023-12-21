  #!/bin/sh

## Cross reference splicing events with Unipro database (*hg38.col.bed)

bedtools intersect -wo -a results/splicing_events.SE.total.HGG.pos.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.HGG.pos.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.HGG.pos.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.HGG.pos.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.HGG.pos.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.HGG.pos.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.HGG.pos.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.HGG.pos.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.HGG.pos.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.HGG.pos.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.SE.total.HGG.neg.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.HGG.neg.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.HGG.neg.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.HGG.neg.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.HGG.neg.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.HGG.neg.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.HGG.neg.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.HGG.neg.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.HGG.neg.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.HGG.neg.intersectUnipDomain.wo.txt
