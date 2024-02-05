#!/usr/bin/perl

## input files
my $es_file = "results/splicing_events.total.HGG.neg.intersectUnip.ggplot.txt";
my $ei_file = "results/splicing_events.total.HGG.pos.intersectUnip.ggplot.txt";

## output files
open(OUT,">results/splice-functional-sites-table.tsv");

print OUT "Gene\tExon_coord\tUpstream_exon_coord\tDownstream_exon_coord\tdPSI\tSite\tPreference\n";

open(FIL, $es_file) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  next unless ($_=~/_/);
  my @cols = split "\t";
  my $splice_id = $cols[0];
  my $dpsi = $cols[1];
  my $site = $cols[2];

  my ($gene, $coords) = split/\:/, $splice_id;
  my ($exon, $upstream, $downstream) = split/\_/,$coords;

  print OUT $gene,"\t",$exon, "\t", $upstream,"\t",$downstream,"\t",$dpsi,"\tSkipping\n";
}


open(FIL, $ei_file) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  next unless ($_=~/_/);
  my @cols = split "\t";
  my $splice_id = $cols[0];
  my $dpsi = $cols[1];
  my $site = $cols[2];

  my ($gene, $coords) = split/\:/, $splice_id;
  my ($exon, $upstream, $downstream) = split/\_/,$coords;

  print OUT $gene,"\t",$exon, "\t", $upstream,"\t",$downstream,"\t",$dpsi,"\tInclusion\n";
}

__DATA__
PDLIM7	177491807-177491976_177491071-177491420_177492404-177492435	Schwannoma
