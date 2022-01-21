#!/usr/bin/perl
###############################################################
# quant_MHC_expr_vs_IR_PSI.pl
# generate table of expression of MHC-related genes and NLRC5
# intron retention
#
# usage: ./quant_MHC_expr_vs_IR_PSI.pl <stranded_counts>
#        <polyA_counts>  <psi_file>
#
# author: Ammar Naqvi
###############################################################

#use warnings;
use strict;

##input files
my ($stranded_rsem, $polyA_rsem, $psi) = ($ARGV[0],$ARGV[1],$ARGV[2]);

##hard-coded list of splice ids to get expression for
my $mhc_class1= "datasets/MHC_class1.other.txt";
my $input= "diffSplicing_cand.filterDiffExpr.txt";



my @samples;
my %sample_psi;
my @genes;
#my (@class1_genes,@class2_genes);
my @gene_list;

my %gene_list;
my %tpm_tumor_str;
my %sample_psi;

my %filter_samples;


my $pct_input = "input/gene_pct.names.txt";
my %ens_to_gene;


#
# #make  gene list
# open (FIL,$input) || die("Cannot Open File");
# while(<FIL>)
# {
#   chomp;
#   my (@cols) = split;
#   my $gene = $cols[0];
#   $gene=~s/\"//g;
#
#   ##save to array and hash to access later
#   push @gene_list, $gene;
#   #print $gene,"\n";
#   $gene_list{$gene} = 1;
# }
# close(FIL);



open(FIL, "input/CLK1_dpsi.midline_HGGs.sign.names.txt") || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  next unless $_=~/BS/;
  my $sample = split;
  #print $sample,"\n";
  $filter_samples{$_} = 1;
}

if($psi=~/(\w+)\_\d+\-\d+\_\d+\-\d+\_\d+\-\d+\./)
{
  #print $psi,"\t",$1."TEST\n";
  push @gene_list, $1;
  #print $gene,"\n";
  $gene_list{$1} = 1;

}
close(FIL);

my $splice_ID = "";
if($psi=~/(\w+\_\d+\-\d+\_\d+\-\d+\_\d+\-\d+)\./)
{
  #print $psi,"\t",$1."TEST\n";
  $splice_ID = $1;
  #print $gene,"\n";
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
  #my ($bs_id, $gene, $splcie_id, $psi) = split "\t";
  next unless $_=~/BS/;
  #next unless $_=~/Ependymoma/;

  my ($bs_id, $psi) = split "\t";
  next unless $filter_samples{$bs_id};

  #print $bs_id,"\t",$psi,"\n";
  $sample_psi{$bs_id} = $psi;
}
close(PSI);

## retrieve expression values and store (stranded) ~/Desktop/AS-DMG/data/stranded_rsem_counts_tab.tsv
my @samples;
open(TPM, $stranded_rsem) || die("Cannot Open File");
while(<TPM>)
{
  chomp;
  if($_=~/transcript_id/)
  {
    my @cols = split;
    shift @cols;
#    shift @cols;

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
    #print "ens\t",$ens,"\n";
    ## clean up gene name
    # $gene=~s/\.\d+//;
    # $gene=~s/ENSG\d+\_//;
    # $gene=~s/\"//g;

    #skip unless gene of interest provided by input
    next unless $ens_to_gene{$ens};

    my $gene = $ens_to_gene{$ens};

    my $i = 0;
    foreach my $tpm(@cols)
    {
      my $bs_id = $samples[$i];
      $i++;

      ##clean up sample ID
      $bs_id=~s/\"//g;

      ##skip sample if PSI event not reported
      next unless $sample_psi{$bs_id};

      ##store expression value by sample and gene
      #print $bs_id,"\t",$gene,"\t",$tpm,"\n";
      $tpm_tumor_str{$bs_id}{$gene}+= $tpm;
    }
  }

}
close(TPM);

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
  #next unless $sample_psi{$sample} > 0;
  #next if ($filter_samples{$sample});
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
    #print $tpm_tumor_str{$sample}{$gene},"\t";
    #print $tpm_tumor_str{$sample}{$gene},"\n";
  }
  print OUT $sample_psi{$sample},"\n";
}

__DATA__
foreach my $sample (keys %sample_psi)
{
  foreach my $gene(@genes)
  {
    next unless $tpm_tumor_str{$sample}{$gene};
    print $sample,"\t";
    print $tpm_tumor_str{$sample}{$gene},"\t";
    print $sample_psi{$sample},"\t";
    print $gene,"\n";
  }
}
