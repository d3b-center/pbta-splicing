#!/bin/sh

cat input/7316_1763.CLK1.aln |  grep "chr2\t" | grep "\ttranscript\t" | awk -F "\t" '{print $NF}' | awk -F ";" '{print $2"\t"$5}' | awk '{print $2"\t"$4}' | perl -pe 's/\"//g' > input/7316_1763.CLK1.processed.txt
cat input/7316_1769.CLK1.aln | grep "chr2\t" | grep "\ttranscript\t" | awk -F "\t" '{print $NF}' | awk -F ";" '{print $2"\t"$5}' | awk '{print $2"\t"$4}' | perl -pe 's/\"//g' > input/7316_1769.CLK1.processed.txt
cat input/KNS42.CLK1.aln | grep "chr2\t" | grep "\ttranscript\t" | awk -F "\t" '{print $NF}' | awk -F ";" '{print $2"\t"$5}' | awk '{print $2"\t"$4}' | perl -pe 's/\"//g' > input/KNS42.CLK1.processed.txt
