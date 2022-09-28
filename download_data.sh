#!/usr/bin/bash

# download files and put in data folder

wget --no-check-certificate https://figshare.com/ndownloader/files/29171019 -O ./data/rMATS_merged.comparison.tsv
wget --no-check-certificate https://figshare.com/ndownloader/files/30452106 -O ./data/rMATS_merged.single.SE.tsv
wget --no-check-certificate https://figshare.com/ndownloader/files/35275717 -O ./data/pbta-histologies.RNA-Seq.initial.tsv
wget --no-check-certificate https://figshare.com/ndownloader/files/35275588 -O ./data/v19_plus_20210311_pnoc_rna.tsv
wget --no-check-certificate https://figshare.com/ndownloader/files/35361607 -O ./data/stranded_rsem_counts_tab.tsv
