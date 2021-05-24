#!/bin/sh

perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep ATRT > results/perc_hist_as.30prev.ATRT.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep Craniopharyn > results/perc_hist_as.30prev.Cran.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep Epend > results/perc_hist_as.30prev.Epend.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep -i gang  > results/perc_hist_as.30prev.Gang.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep -i hgat  > results/perc_hist_as.30prev.HGAT.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep -i lgat  > results/perc_hist_as.30prev.LGAT.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep -i Medu  > results/perc_hist_as.30prev.medul.tsv

