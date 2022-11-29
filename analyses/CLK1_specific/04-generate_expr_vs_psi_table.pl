#!/usr/bin/perl
#
###############################################################
# quant_expr_vs_psi_PCT.pl
# generate table of expression of pct-expr and PSI
#
# usage: perl quant_expr_vs_psi_PCT.pl <stranded_counts>
#        <polyA_counts>  <psi_file>
#
# author: Ammar Naqvi
###############################################################

#use warnings;
use strict;

## input files inluding rsem transcript counts and PSI tables previously made
my ($stranded_rsem, $polyA_rsem, $psi) = ($ARGV[0],$ARGV[1],$ARGV[2]);

##hard-coded list of splice ids to get expression for
my $input= "input/diffSplicing_cand.filterDiffExpr.txt";

my @samples;
my %sample_psi;
my @genes;
my @gene_list;
my %gene_list;
my %tpm_tumor_str;
my %sample_psi;
my %filter_samples;

## file of gene and their protein coding transcripts to be used for filtering below
my $pct_input = "input/gene_pct.names.txt";
my %ens_to_gene; # hash to convert transcript ID to gene names

## get samples with signficant CLK1 Exon 4 changes to filter on later
open(FIL, "input/CLK1_dpsi.midline_HGGs.sign.names.txt") || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  next unless $_=~/BS/;
  my $sample = split;
  #print $sample,"\n";
  $filter_samples{$_} = 1;
}

## extract gene name
if($psi=~/(\w+)\_\d+\-\d+\_\d+\-\d+\_\d+\-\d+\./)
{
  push @gene_list, $1;
  $gene_list{$1} = 1;
}
close(FIL);

## extract splice variant ID
my $splice_ID = "";
if($psi=~/(\w+\_\d+\-\d+\_\d+\-\d+\_\d+\-\d+)\./)
{
  $splice_ID = $1;
}

## make hash of gene names and protein coding transcripts
open(FIL, $pct_input) || ("Cannot Open File");
while(<FIL>)
{
  chomp;
  my @cols = split;
  my $parent = $cols[0];
  my $name = $cols[1];

  my ($title, $ens)  = split/\=/,$parent;
  $ens =~s/\.\d+//;

  my ($title, $gene) = split/\=/,$name;
  $gene =~s/\=//;

  next unless $gene_list{$gene};
  $ens_to_gene{$ens} = $gene;
}
close(FIL);

## save and store PSI values for each sample
open(PSI,$psi) || die("Cannot Open File");
while(<PSI>)
{
  chomp;
  next unless $_=~/BS/;
  my ($bs_id, $psi) = split "\t";

  next unless $filter_samples{$bs_id}; ## filter non-midline  HGGs
  $sample_psi{$bs_id} = $psi;
}
close(PSI);

## retrieve expression values and store (stranded)
my @samples;
open(TPM, $stranded_rsem) || die("Cannot Open File");
while(<TPM>)
{
  chomp;
  if($_=~/transcript_id/)
  {
    my @cols = split;
    shift @cols;

    @samples = @cols; ## store sample name in array
  }
  else
  {
    my @cols = split;
    shift @cols; ## remove blank column

    my $ens = shift @cols;
    if($ens=~/(ENST\d+)\./)
    {
      $ens = $1;
    }

    #skip unless gene  of interest provided by input
    next unless $ens_to_gene{$ens};
    my $gene = $ens_to_gene{$ens};
    my $i = 0;

    foreach my $tpm(@cols)
    {
      my $bs_id = $samples[$i];
      $i++;
      $bs_id=~s/\"//g; ##clean up sample ID

      ## skip sample if PSI event not reported
      next unless $sample_psi{$bs_id};

      ## store expression value by sample and gene
      $tpm_tumor_str{$bs_id}{$gene}+= $tpm;
    }
  }

}
close(TPM);

## check to see if results directories are made
if (!-d "scr")
{
  mkdir "scr";
}

my $outfile = "scr\/".$splice_ID."_vs_expr.tab.txt";
open(OUT, ">".$outfile);

##make table and report psi and expression info // class I or II
print OUT "sample\t";
foreach my $gene(@gene_list)
{
  print OUT $gene."\t";
}
print OUT "PSI\n";

foreach my $sample (sort keys %sample_psi)
{
  my $gene_name = $gene_list[0];
  next unless $tpm_tumor_str{$sample}{$gene_name} > 0;
  print OUT $sample,"\t";
  foreach my $gene(@gene_list)
  {
    if($tpm_tumor_str{$sample}{$gene} > 0)
    {
      print OUT $tpm_tumor_str{$sample}{$gene},"\t";
    }
    else{
      print OUT "0\t";
    }

  }
  print OUT $sample_psi{$sample},"\n";
}
