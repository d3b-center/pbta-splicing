#!/usr/bin/perl

## input files
my $es_file = "results/unique_events-es.tsv";
my $ei_file = "results/unique_events-ei.tsv";

## output files
open(OUT_ES,">results/hist-spec-es-table.tsv");
open(OUT_EI,">results/hist-spec-ei-table.tsv");

print OUT_ES "Gene\tExon_coord\tUpstream_exon_coord\tDownstream_exon_coord\thHistology\n";
print OUT_EI "Gene\tExon_coord\tUpstream_exon_coord\tDownstream_exon_coord\tHistology\n";

open(FIL, $es_file) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  next unless ($_=~/_/);
  my @cols = split "\t";
  my $splice_id = $cols[0];
  my $histology = $cols[1];

  my ($gene, $coords) = split/\:/, $splice_id;
  my ($exon, $upstream, $downstream) = split/\_/,$coords;

  print OUT_ES $gene,"\t",$exon, "\t", $upstream,"\t",$downstream,"\t",$histology,"\n";
}


open(FIL, $ei_file) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  next unless ($_=~/_/);
  my @cols = split "\t";
  my $splice_id = $cols[0];
  my $histology = $cols[1];

  my ($gene, $coords) = split/\:/, $splice_id;
  my ($exon, $upstream, $downstream) = split/\_/,$coords;

  print OUT_EI $gene,"\t",$exon, "\t", $upstream,"\t",$downstream,"\t",$histology,"\n";
}

__DATA__
PDLIM7	177491807-177491976_177491071-177491420_177492404-177492435	Schwannoma
