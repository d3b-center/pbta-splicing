#!/usr/bin/perl

use List::Util qw(max);

################################################################################
# create_matrix_of_PSI_SE_gene.pl
# create tsv file of PSI for each sample, to be used for consensus clustering
# written by Ammar Naqvi
#
# usage: ./create_matrix_of_PSI_removeDups.pl <histology file> <rMATs_file_paths.txt>
#                                             <primary tumor file> <primary plus file>
################################################################################
my ($histology,$rmats_tsv,$primary_tumor_dat,$primary_tumor_plus_dat) = ($ARGV[0], $ARGV[1], $ARGV[2],$ARGV[3],$ARGV[4]);
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);

my (%init_patientIDs, %progr_patientIDs, %recur_patientIDs, %unavail_patientIDs);
my %patient_id_count;

my %primary_initial_sample_list;

## store primary tumor samples
open(FIL,$primary_tumor_dat) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my @header = split "\t";
  my $bs_id = $header[1];
  my $cohort = $header[2];
  my $exp_strategy = $header[4];
  my $tumor_descr = $header[5];

  next unless ($cohort=~/PBTA/);

  $primary_initial_sample_list{$bs_id} = $bs_id;

}
close(FIL);

## store primary tumor samples
open(FIL,$primary_tumor_plus_dat) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my @header = split "\t";
  my $bs_id = $header[1];
  my $cohort = $header[2];
  my $exp_strategy = $header[4];
  my $tumor_descr = $header[5];

  next unless ($cohort=~/PBTA/);

  $primary_initial_sample_list{$bs_id} = $bs_id;

}
close(FIL);


# annotate histology file #
  # hash with BS and disease
  # make arrays of each histology of BS IDs

open(FIL,$histology) || die("Cannot Open File");
  while(<FIL>)
  {
    chomp;
    my @cols       = split "\t";
    my $hist       = $cols[53];
    my $bs_id      = $cols[0];
    my $patient_id = $cols[3];
    my $CNS_region = $cols[32];


    next unless ($primary_initial_sample_list{$bs_id});

    ## filter histologies of interests
    next unless ( ($hist=~/HGAT/)  ||
                  ($hist=~/LGAT/)  ||
                  #($hist=~/Oligodendroglioma/) ||
                  ($hist=~/Medulloblastoma/)   ||
                  ($hist=~/Ganglioglioma/)  ||
                  ($hist=~/Ependymoma/)||
                  ($hist=~/ATRT/)  ||
                  ($hist=~/Craniopharyngioma/) );

    ## convert histology names
    $hist =~s/Oligodendroglioma/OGG/;
    $hist =~s/Medulloblastoma/MB/;
    $hist =~s/Ganglioglioma/GNG/;
    $hist =~s/Ependymoma/EPN/;
    $hist =~s/Craniopharyngioma/CPG/;
    $hist =~s/HGAT/HGG/;
    $hist =~s/LGAT/LGG/;

  ## make an array and store histology information and BS IDs
  push @broad_hist, $hist;
  push @bs_ids, $bs_id;

  $bs_id_hist{$bs_id} = $hist;

  ## store total number of histologies
  $hist_count{$hist}++;
  push @{$histology_ids{$hist}}, $bs_id;

  $cns_regions{$bs_id} = $CNS_region;

  ## histology counter for downstream analysis
  $hist_count{$hist}++;

}

close(FIL);


my %inc_levels_gene;
my @genes;
my %hist_check_gene;

# store each library (using file name)
# for each line, make unique splice id
  # gene + chr + SE + UEx + DEx   (make array)
  # store inclusion level of SE   (hash of splice ids and library)
  ## process rMATS output (may take awhile) from merged input file
  print "\tprocessing rMATs results...\n";
  open(FIL, "gunzip -c $rmats_tsv |") || die ("can’t open $rmats_tsv");
  while(<FIL>)
  {
    chomp;

    ##only look at exon splicing events
    next unless($_=~/^SE/);

    my @cols = split "\t";
    my $bs_id = $cols[1];
    my $ctrl   = $cols[2];

    ## get histology for BS ID#
    my $hist_of_sample = $bs_id_hist{$bs_id};

    ## get gene name
    my $gene         = $cols[4];
    $gene=~s/\"//g; # remove quotation marks

    ## retrieve exon coordinates
    my $chr          = $cols[5];
    my $str          = $cols[6];
    my $exonStart    = $cols[7]+1; ## its 0-based so add 1 to overlap with uni-prot coord
    my $exonEnd      = $cols[8];
    my $upstreamES   = $cols[15];
    my $upstreamEE   = $cols[16];
    my $downstreamES = $cols[17];
    my $downstreamEE = $cols[18];

    ## retrieve inclusion level and junction count info
    my $inc_level_tumor = $cols[33];
    my $tumor_IJC        = $cols[25];
    my $tumor_SJC        = $cols[26];


    #get lengths of exons
    my $inc_len  = $cols[29];
    my $skip_len = $cols[30];


    ## only look at strong changes, remove any dPSI < .20 and tumor junction reads > 10 reads
  #  print "incl_level ".$inc_level_tumor,"\n";
    next unless ( ($inc_level_tumor >= .10) && ($inc_level_tumor <= .90) );
    next unless ( ($tumor_IJC + $tumor_SJC) >= 100 ) ;

    # annotate and re-name splice event IDs
    my $splice_id= $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;

    ## store inclusion levels of splicing event and BS ID#
    $inc_levels{$splice_id}{$bs_id} = $inc_level_tumor;

    ##store array of inclusion levels per gene to access later
    push @{$inc_levels_gene{$gene}{$bs_id}},$inc_level_tumor;
    push @genes, $gene;

    #print "gene:".$gene,"\t",$bs_id,"\n";
    ## get histology for BS ID#
    my $hist_of_sample = $bs_id_hist{$bs_id};

    ## store number of events per histology
    $hist_check{$splice_id}{$hist_of_sample}++;
    push @splicing_events, $splice_id;

    ##store gene / hsitology
    $hist_check_gene{$gene}{$hist_of_sample}++;

  }
close(FIL);


my @broad_hist_uniq      = do { my %seen; grep { !$seen{$_}++ } @broad_hist };
my @bs_ids_uniq          = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };
my @genes_uniq = do { my %seen; grep { !$seen{$_}++ } @genes };


# go through each splicing event
  # order by disease library
  # print splice id, inc level
my $out_file = "input/pan_cancer_splicing_SE.gene.dev.txt";
open(OUT,">",$out_file) || die("Cannot Open File");

## save and print header info for output file
print OUT "Splice_ID";
foreach my $hist (sort @broad_hist_uniq)
{
  my $i = 1;
  foreach my $sample (@{$histology_ids{$hist}})
  {
  #  last if($i>30);
    print OUT "\t";
    ##print histology name (renamed)
    #print $hist."_".$i;
    ##print sample name
    print OUT $sample;
    $i++;

  }
}
print OUT "\n";

## print each event and threshold value to output file
foreach my $event (@genes_uniq)
{
  ## @{$histology_ids{$broad_hist}}
  my $thresh_for_hist_prev = 1;
  #print "event:",$event,"\n";
  foreach my $hist (sort @broad_hist_uniq){

      ##threshold for prevalence per sample (%)
      #my $hist_prev_thr = .20*($hist_count{$hist});
      my $hist_prev_thr = 2;  ## must be recurrent event in histology (n>=2)
      my $prev_in_hist = 0;
      if($hist_check_gene{$event}{$hist}) { $prev_in_hist = $hist_check{$event}{$hist} };
      #print $prev_in_hist,"CHECK\t",$hist_prev_thr,"\n";
      if($prev_in_hist < $hist_prev_thr) {
        #$thresh_for_hist_prev = 0;
      }
  }
  next if ($thresh_for_hist_prev == 0);
  print OUT $event;

  foreach my $hist (sort @broad_hist_uniq)
  {

    my $i = 1;
    foreach my $sample (@{$histology_ids{$hist}})
    {

    #  last if($i>30); ## for testing purposes

      print OUT "\t";
      #print $sample;
      #my $greatest_psi=  ((sort {$b <=> $a} $inc_levels_gene{$event}{$sample})[0]);
      my $greatest_psi=  max(@{$inc_levels_gene{$event}{$sample}});

      if($greatest_psi>0)
      {
        #print $sample.":".$inc_levels{$event}{$sample};
        print OUT $greatest_psi;
      }
      else
      {
        # print $sample.":0";
        print OUT "0";
      }
      $i++;
    }
  }
  print OUT "\n";
}
