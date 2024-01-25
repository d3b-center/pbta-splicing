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

bedtools intersect -wo -a results/splicing_events.SE.total.H33.pos.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.H33.pos.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H33.pos.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.H33.pos.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H33.pos.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.H33.pos.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H33.pos.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.H33.pos.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H33.pos.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.H33.pos.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.SE.total.H33.neg.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.H33.neg.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H33.neg.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.H33.neg.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H33.neg.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.H33.neg.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H33.neg.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.H33.neg.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H33.neg.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.H33.neg.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.SE.total.H31.pos.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.H31.pos.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H31.pos.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.H31.pos.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H31.pos.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.H31.pos.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H31.pos.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.H31.pos.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H31.pos.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.H31.pos.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.SE.total.H31.neg.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.H31.neg.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H31.neg.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.H31.neg.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H31.neg.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.H31.neg.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H31.neg.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.H31.neg.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H31.neg.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.H31.neg.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.SE.total.H3WT.pos.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.H3WT.pos.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H3WT.pos.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.H3WT.pos.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H3WT.pos.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.H3WT.pos.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H3WT.pos.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.H3WT.pos.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H3WT.pos.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.H3WT.pos.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.SE.total.H3WT.neg.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.H3WT.neg.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H3WT.neg.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.H3WT.neg.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H3WT.neg.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.H3WT.neg.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H3WT.neg.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.H3WT.neg.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.SE.total.H3WT.neg.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.H3WT.neg.intersectUnipDomain.wo.txt
