#!/bin/sh

## Cross reference splicing events with Unipro database (*hg38.col.bed)
bedtools intersect -wo -a results/splicing_events.total.pos.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.pos.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.total.pos.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.pos.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.total.pos.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.pos.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.total.pos.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.pos.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.total.pos.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.pos.intersectUnipDomain.wo.txt

bedtools intersect -wo -a results/splicing_events.total.neg.bed -b input/unipMod.hg38.col.bed       |sort -u > results/splicing_events.total.neg.intersectUnipMod.wo.txt
bedtools intersect -wo -a results/splicing_events.total.neg.bed -b input/unipOther.hg38.col.bed     |sort -u > results/splicing_events.total.neg.intersectUnipOther.wo.txt
bedtools intersect -wo -a results/splicing_events.total.neg.bed -b input/unipDisulfBond.hg38col.bed |sort -u > results/splicing_events.total.neg.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a results/splicing_events.total.neg.bed -b input/unipLocSignal.hg38.col.bed |sort -u > results/splicing_events.total.neg.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a results/splicing_events.total.neg.bed -b input/unipDomain.hg38.col.bed    |sort -u > results/splicing_events.total.neg.intersectUnipDomain.wo.txt
