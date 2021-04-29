#!/bin/sh

echo "parse out all relevant CLK1 inclusion events";
echo "Sample\tGene\tExon_st\tExon_end\tUps_st\tUps_end\tDown_st\tDown_end\tInclusion" > results/CLK1_IncLevel_brain_tumors.tab.v2.txt
grep "\"CLK1\"" ~/Desktop/pan_cancer_rmats/*freads10.txt |  awk -F "\t" '{ print $0} ' | grep "200860124\t200860215\t200859679\t200859746\t200861237\t200861466"| awk -F "\t" '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$21}' | perl -pe 'if($_=~/\.(BS\_\w+)\.non/){ print $1,"\t",$_}' | awk -F "\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$NF}' | perl -pe 's/\"//g' | grep -v "User" >> results/CLK1_IncLevel_brain_tumors.tab.v2.txt
