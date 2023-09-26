#!/bin/sh

## get SRSF11  splicing event in tumor and ctrl samples
cat ../../data/rMATS_merged.comparison.tsv.gz | grep "^SE" | grep "SRSF11" | grep -P "70231858.0\\t70231975.0" | grep -P  "70228421.0\\t70228555.0" | grep -P  "70232267.0\\t70232375.0"| awk '{print $2"\t"$35"\tinclusion"}' | grep "BS" | sort -u > input/se.srsf11.incl.txt
cat ../../data/rMATS_merged.comparison.tsv.gz | grep "^SE" | grep "SRSF11" | grep  -P "70231858.0\\t70231975.0" | grep -P  "70228421.0\\t70228555.0" | grep -P "70232267.0\\t70232375.0" | awk '{print $2"\t"1-$35"\tskip"}' | grep BS | sort -u >> input/se.srsf11.incl.txt

cat ../../data/rMATS_merged.comparison.tsv.gz | grep "^SE" | grep "SRSF11" | grep -P "70231858.0\\t70231975.0" | grep -P "70228421.0\\t70228555.0" | grep  -P "70232267.0\\t70232375.0" | awk '{print "control-BS\t"$34"\tinclusion"}' | sort -nr | head -n 1 > input/se.srsf11.ctrl.incl.txt
cat ../../data/rMATS_merged.comparison.tsv.gz | grep "^SE" | grep "SRSF11" | grep -P "70231858.0\\t70231975.0" | grep  -P "70228421.0\\t70228555.0" | grep -P "70232267.0\\t70232375.0" | awk '{print "control-BS\t"1-$34"\tskip"}' | sort -n | head -n 1 >> input/se.srsf11.ctrl.incl.txt
