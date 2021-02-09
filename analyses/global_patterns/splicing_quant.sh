#!/bin/bash

out=$1
#echo $out
cat ~/Desktop/rmats_run/A5SS/*MATS.JC.txt | awk '{if(  ($19<=0.05) && ( ($NF>=.10) || ($NF <=-.10)) ){ print $0}}' | awk '{print $3"_"$6"-"$7"_"$8"-"$9"_"$10"-"$11}' | sort -u | wc -l > $out
cat ~/Desktop/rmats_run/A3SS/*MATS.JC.txt | awk '{if(  ($19<=0.05) && ( ($NF>=.10) || ($NF <=-.10)) ){ print $0}}' | awk '{print $3"_"$6"-"$7"_"$8"-"$9"_"$10"-"$11}' | sort -u | wc -l >> $out
cat ~/Desktop/rmats_run/RI/*MATS.JC.txt | awk '{if(  ($19<=0.05) && ( ($NF>=.10) || ($NF <=-.10)) ){ print $0}}' | awk '{print $3"_"$6"-"$7"_"$8"-"$9"_"$10"-"$11}' | sort -u | wc -l >> $out
cat ~/Desktop/rmats_run/SE/*MATS.JC.txt | awk '{if(  ($19<=0.05) && ( ($NF>=.10) || ($NF <=-.10)) ){ print $0}}' | awk '{print $3"_"$6"-"$7"_"$8"-"$9"_"$10"-"$11}' | sort -u | wc -l >> $out
cat ~/Desktop/rmats_run/MXE/*MATS.JC.txt | awk '{if(  ($21<=0.05) && ( ($NF>=.10) || ($NF <=-.10)) ){ print $0}}' | awk '{print $3"_"$6"-"$7"_"$8"-"$9"_"$10"-"$11}' | sort -u | wc -l >> $out
