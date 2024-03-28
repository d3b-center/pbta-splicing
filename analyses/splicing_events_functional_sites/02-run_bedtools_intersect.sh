  #!/bin/sh

## Cross reference splicing events with Unipro database (*hg38.col.bed)

bedtools intersect -wo -a "results/splicing_events.SE.total.Other high-grade glioma.pos.bed" -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.other-HGG.pos.intersectUnipMod.wo.txt
bedtools intersect -wo -a "results/splicing_events.SE.total.Other high-grade glioma.pos.bed" -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.other-HGG.pos.intersectUnipOther.wo.txt
bedtools intersect -wo -a "results/splicing_events.SE.total.Other high-grade glioma.pos.bed" -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.other-HGG.pos.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a "results/splicing_events.SE.total.Other high-grade glioma.pos.bed" -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.other-HGG.pos.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a "results/splicing_events.SE.total.Other high-grade glioma.pos.bed" -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.other-HGG.pos.intersectUnipDomain.wo.txt

bedtools intersect -wo -a "results/splicing_events.SE.total.Other high-grade glioma.neg.bed"  -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.other-HGG.neg.intersectUnipMod.wo.txt
bedtools intersect -wo -a "results/splicing_events.SE.total.Other high-grade glioma.neg.bed"  -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.other-HGG.neg.intersectUnipOther.wo.txt
bedtools intersect -wo -a "results/splicing_events.SE.total.Other high-grade glioma.neg.bed" -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.other-HGG.neg.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a "results/splicing_events.SE.total.Other high-grade glioma.neg.bed" -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.other-HGG.neg.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a "results/splicing_events.SE.total.Other high-grade glioma.neg.bed" -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.other-HGG.neg.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.SE.total.DMG.pos.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.DMG.pos.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.DMG.pos.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.DMG.pos.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.DMG.pos.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.DMG.pos.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.DMG.pos.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.DMG.pos.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.DMG.pos.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.DMG.pos.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.SE.total.DMG.neg.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.DMG.neg.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.DMG.neg.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.DMG.neg.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.DMG.neg.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.DMG.neg.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.DMG.neg.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.DMG.neg.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.DMG.neg.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.DMG.neg.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.SE.total.pos.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.pos.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.pos.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.pos.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.pos.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.pos.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.pos.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.pos.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.pos.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.pos.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.SE.total.neg.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.neg.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.neg.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.neg.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.neg.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.neg.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.neg.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.neg.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.neg.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.neg.intersectUnipDomain.wo.txt
