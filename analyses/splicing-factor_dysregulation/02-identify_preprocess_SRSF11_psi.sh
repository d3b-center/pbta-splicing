#!/bin/sh

cat ../../data/rMATS_merged.comparison.tsv.gz | zcat | grep ^SE | grep SRSF11 | grep "70231858\t70231975" | grep "70228421\t70228555" | grep "70232267\t70232375" | awk '{print $2"\t"$35"\tinclusion"}' | grep BS | sort -u > input/se.srsf11.incl.txt
cat ../../data/rMATS_merged.comparison.tsv.gz | zcat | grep ^SE | grep SRSF11 | grep "70231858\t70231975" | grep "70228421\t70228555" | grep "70232267\t70232375" | awk '{print $2"\t"1-$35"\tskip"}' | grep BS | sort -u >> input/se.srsf11.incl.txt
cat ../../data/rMATS_merged.comparison.tsv.gz | zcat | grep ^SE | grep SRSF11 | grep "70231858\t70231975" | grep "70228421\t70228555" | grep "70232267\t70232375" | awk '{print $3"\t"$34"\tinclusion"}' | sort -u > input/se.srsf11.ctrl.incl.txt
cat ../../data/rMATS_merged.comparison.tsv.gz | zcat | grep ^SE | grep SRSF11 | grep "70231858\t70231975" | grep "70228421\t70228555" | grep "70232267\t70232375" | awk '{print $3"\t"1-$34"\tskip"}' | sort -u >> input/se.srsf11.ctrl.incl.txt

