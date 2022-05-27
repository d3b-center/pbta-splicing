#!/usr/bin/perl

use Statistics::Lite qw(:all);

## input files
my $si_tab   = "results/splicing_index.total.txt"; # splicing index tab
my $rsem_tab = "../../data/stranded_rsem_counts_tab.tsv"; ## rsem counts

my %filter_samples;
my (@low_si_samples, @high_si_samples);
my @genes;
my %filter_samples_str;

open(FIL, $si_tab) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  next if ($_=~/HGAT/);

  my @cols = split;
  my $si = $cols[5];
  my $sample = $cols[0];
  $filter_samples{$sample} = 1;
  #print "**",$sample,"*".$si,"\n";

  #if ($si > 0.07229791)
  if ($si > 0.04789595) #non HGAT

  {
    push @high_si_samples, $sample;
  }

  #if ($si < 0.02699138)
  if ($si < 0.01242844) ## non HGAT
  {
    push @low_si_samples, $sample;
  }
}

my @samples;
open(TPM, $rsem_tab) || die("Cannot Open File");
while(<TPM>)
{
  chomp;
  if($_=~/gene_id/)
  {
    my @cols = split;
    shift @cols;
    shift @cols;

    @samples = @cols; ## store sample name in array
  }
  else
  {
    my @cols = split;
    shift @cols; ## remove blank column

    my $gene = shift @cols;

    ## clean up gene name
    #$gene=~s/\.\d+//;
    #$gene=~s/ENSG\d+\_//;
    $gene=~s/\"//g;
    push @genes, $gene;
    #skip unless gene of interest provided by input
    #next unless $gene_list{$gene};

    #print $gene;

    my $i = 0;
    foreach my $tpm(@cols)
    {
      my $bs_id = $samples[$i];
      $i++;

      ##clean up sample ID
      $bs_id=~s/\"//g;

      ##skip sample if IR event not reported
    #  print "CHECK1_".$bs_id."\n";
      #next unless $filter_samples{$bs_id};
    #  print "CHECK2_".$bs_id."\n";

      $filter_samples_str{$bs_id} = 1;
      ##store expression value by sample and gene
      #print "\t".$tpm;
      $tpm_tumor_str{$bs_id}{$gene}+= $tpm;
    }
    #print "\n";
  }

}
close(TPM);

print "gene";
foreach my $sample(@high_si_samples)
{
  print "\tHigh_", $sample if $filter_samples_str{$sample};
}

foreach my $sample(@low_si_samples)
{
  print "\tLow_", $sample if $filter_samples_str{$sample};
}
print "\n";

foreach my $gene(@genes)
{
  print $gene;
  foreach my $high_sample(@high_si_samples)
  {
    next unless $filter_samples_str{$high_sample};

    print "\t";
    print $tpm_tumor_str{$high_sample}{$gene};
  }
  foreach my $low_sample(@low_si_samples)
  {
    next unless $filter_samples_str{$low_sample};

    print "\t";
    print $tpm_tumor_str{$low_sample}{$gene};
  }
  print "\n";
}


__DATA__
results/splicing_index.total.txt
my $mean     = mean  (@{$splicing_psi{$splice_event}});
  my $count    = count (@{$splicing_psi{$splice_event}});
  my $std_dev  = stddev(@{$splicing_psi{$splice_event}});

