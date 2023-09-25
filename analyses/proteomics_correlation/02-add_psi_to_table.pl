#!/usr/bin/perl

my %splice_id_entry;
my %hgg_midline_samples;
my %genes_filter;
my %inc_levels;

my @samples;

my $histology = "../../data/histologies.tsv";
my $rmats_tsv = "../../data/splice-events-rmats.tsv.gz";
my $rmats_tsv_v2 = "../../data/test.tsv";


## open histology file and annotate DMG samples
open(FIL,$histology) || die("Cannot Open File");
while(<FIL>)
{
  chomp;

  ## only select RNA-seq and PBTA samples
  next unless ($_=~/RNA-Seq/);
  next unless ($_=~/PBTA/);

  my @cols       = split "\t";
  my $hist       = $cols[53];
  my $bs_id      = $cols[0];
  my $patient_id = $cols[3];
  my $CNS_region = $cols[25];

  ##keep only HGG midline samples
  next unless $hist=~/HGAT/;
  #next unless($CNS_region=~/Midline/);
  #print $hist,"\t",$CNS_region,"\t",$bs_id,"\n";

  $hgg_midline_samples{$bs_id} = $bs_id;
}
close(FIL);

open(FIL,"results/mixed_events.tsv");
while(<FIL>)
{
  chomp;
  next if($_=~/gene/);
  my @cols = split;
  my $gene = $cols[0];
  my $splice_id = $cols[5];
  #print $splice_id,"\t",$_,"\n";
  $splice_id_entry{$splice_id} = $_;
$genes_filter{$gene} = $gene;

}

## process rMATS output (may take awhile) from merged input file
print "processing rMATs results...\n";
open(FIL, "gunzip -c $rmats_tsv |") || die ("can’t open $rmats_tsv");
open(FIL, $rmats_tsv_v2 ) || die ("can’t open $rmats_tsv_v2");

while(<FIL>)
{
  chomp;

  ##only look at exon splicing events
  next unless($_=~/^SE/);

  my @cols = split "\t";
  my $bs_id = $cols[1];
  my $ctrl   = $cols[2];

  ## skip non-midline HGGs and non BS ctrl comparisons



  next unless $hgg_midline_samples{$bs_id};


  ## get gene name
  my $gene         = $cols[4];
  $gene=~s/\"//g; # remove quotation marks

  next unless $genes_filter{$gene};

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

   $exonStart    =~s/\.0//;
   $exonEnd      =~s/\.0//;
   $upstreamES   =~s/\.0//;
   $upstreamEE   =~s/\.0//;
   $downstreamES =~s/\.0//;
   $downstreamEE =~s/\.0//;


   ## only look at tumor junction reads > 10 reads
   #next unless ( ($tumor_IJC + $tumor_SJC) >= 10 ) ;

   # annotate and re-name splice event IDs
   my $splice_id= $gene."_".$exonStart."-".$exonEnd;
   print $splice_id,"\n";

   next unless($splice_id_entry{$splice_id});


   ## store inclusion levels of splicing event and BS ID#
   $inc_levels{$splice_id}{$bs_id} = $inc_level_tumor;

   ##store array of inclusion levels per gene to access later
   push @{$inc_levels_gene{$gene}{$bs_id}},$inc_level_tumor;
   push @genes, $gene;

   #print "gene:".$gene,"\t",$bs_id,"\t",$splice_id,"\n";
   my $hist_of_sample = $bs_id_hist{$bs_id};

   ## store number of events per histology
   push @samples, $bs_id;
   push @splicing_events, $splice_id;



 }
close(FIL);

## make arrays unique
my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

open(TAB, ">results/splicing_events-psi.tsv");


print TAB"gene\tchr\tstart\tend\tstrand\tSpliceID\tdPSI_Skip\tdPSI_Inc\t";
print TAB join "\t", @samples_uniq,"\n";

foreach my $splice_event(@splicing_events_uniq)
{
  print TAB$splice_id_entry{$splice_event};
  print TAB"\t";
  foreach my $sample(@samples_uniq)
  {
    my $psi = $inc_levels{$splice_event}{$sample};
    #print $sample.":";
    print TAB sprintf("%.2f",$psi),"\t";
  }
  print TAB "\n";
}
