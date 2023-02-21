#!/usr/bin/perl

#use strict;
use warnings;

## check to see if results directories are made
if (!-d "results")
{
  mkdir "results";
}

if (!-d "results/results_diff"){
 mkdir "results/results_diff";
}

while(<>)
{
  chomp;
  my @cols          = split/\_/;
  my $gene          = $cols[0];
  my $exon_coord    = $cols[1];
  my $upstream_ec   = $cols[2];
  my $downstream_ec = $cols[3];

  my ($exon_start,$exon_end) = split/\-/, $exon_coord;
  my $exon_start_minus = $exon_start-1;
  $upstream_ec     =~s/\-/\\t/;
  $downstream_ec   =~s/\-/\\t/;
  my $exon_coords_mod = $exon_start_minus."\\t".$exon_end;

  my $out_file = "results\/results_diff\/";
  #print $out_file,"\t",$gene,"\n";

	$out_file .= $gene."_".$exon_coords_mod."_".$upstream_ec."_".$downstream_ec."\.txt";
  $out_file =~s/\\t/\-/g;
  #system("cat ~/Desktop/pbta-splicing/data/merge_rMATS_splicing.SE.single.tsv | grep $gene | grep \"$exon_coords_mod\" | grep \"$upstream_ec\" | grep \"$downstream_ec\" | awk \'{print \$2\"\t\"\$34}\' > test_lsv.txt");
  system ("cat ../../data/rMATS_merged.single.SE.tsv.gz | zcat | grep $gene | grep \"$exon_coords_mod\" | grep \"$upstream_ec\" | grep \"$downstream_ec\" | awk \'{print \$2\"\\t\"\$34}\' > $out_file");
  print $out_file,"\n";
#	 print "\n";



}
