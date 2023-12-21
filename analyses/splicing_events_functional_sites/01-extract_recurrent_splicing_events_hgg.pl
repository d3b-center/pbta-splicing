#!/usr/bin/perl

use Statistics::Lite qw(:all);
############################################################################################################
# 03-generate_splicing-index_and_diff-events_table.SE.pl
#
# Compute splicing index for each sample and generate splicing burden index tables
############################################################################################################
my ($histology,$rmats_tsv,$primary_tumor_file, $splice_case) = ($ARGV[0], $ARGV[1], $ARGV[2],$ARGV[3]);
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);
my @splicing_events;
my %splicing_psi;
my %cns_regions;

$splice_case=~s/\s+//;
print $splice_case,"\n";

## check to see if splicing case is correct
unless ($splice_case=~/SE$|A3SS$|A5SS$|RI$/)
{
  die("Splicing case does not exist, please use either 'SE', 'A5SS', 'A3SS', or 'RI'");
}
# annotate histology file #
  # hash with BS and disease
  # make arrays of each histology of BS IDs

  my %primary_initial_sample_list;


  ## store primary tumor samples
  open(FIL,$primary_tumor_file) || die("Cannot Open File");
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
      #my $hist       = $cols[53];
      my $hist = $cols[-2];
      my $bs_id      = $cols[0];
      my $CNS_region = $cols[25];


      next unless ($primary_initial_sample_list{$bs_id});



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


#Exon-specific
# For every alternatively spliced exon
#-------------------------------------
  #1. Compute mean of exon inclusion (PSI)
  #2. Compute STD

#-Look at each tumor --> is it 2 standard deviations from the mean?
  #--Yes --> aberrant
  #--No  --> non-aberrant

# store each library (using file name)
# for each line, make unique splice id
  # gene + chr + SE + UEx + DEx   (make array)
  # store total splicing events and those that pass threshold for each sample

my (%splice_totals_per_sample, %splice_totals,%chr, %str);
my %splice_event_hist;

## process rMATS output (may take awhile) from merged input file
print "processing rMATs results... ".(localtime)."\n";
open(FIL, "gunzip -c $rmats_tsv |") || die ("can’t open $rmats_tsv");
#open(FIL, "scr/events-rmats-head.tsv") || die ("can’t open $rmats_tsv");


while(<FIL>)
{
  chomp;

  ##only look at exon splicing events
  next unless($_=~/^$splice_case/);

  my @cols  = split "\t";
  my $bs_id = $cols[1];
  my $ctrl  = $cols[2];

  ## get gene name
  my $gene         = $cols[4];
  $gene=~s/\"//g; # remove quotation marks

  ## retrieve exon coordinates
  my $chr          = $cols[5];
  my $str          = $cols[6];

  my $Start    = "";
  my $End      = "";
  my $prevES   = "";
  my $prevEE   = "";
  my $nextES = "";
  my $nextEE = "";

  if($splice_case=~/SE/)
  {
    $Start    = $cols[7]+1; ## its 0-based so add 1
    $End      = $cols[8];
    $prevES   = $cols[15];
    $prevEE   = $cols[16];
    $nextES = $cols[17];
    $nextEE = $cols[18];
  }
  elsif($splice_case=~/RI/)
  {
    $Start    = $cols[13]+1; ## its 0-based so add 1
    $End      = $cols[14];
    $prevES   = $cols[15];
    $prevEE   = $cols[16];
    $nextES = $cols[17];
    $nextEE = $cols[18];
  }
  elsif($splice_case=~/A5SS/)
  {
    $Start    = $cols[19]+1; ## its 0-based so add 1
    $End      = $cols[20];
    $prevES   = $cols[21];
    $prevEE   = $cols[22];
    $nextES = $cols[23];
    $nextEE = $cols[24];

  }
  else{
    $Start    = $cols[19]+1; ## its 0-based so add 1
    $End      = $cols[20];
    $prevES   = $cols[21];
    $prevEE   = $cols[22];
    $nextES = $cols[23];
    $nextEE = $cols[24];
    #print $Start,"\t",$End,"\t",$prevES,"\t",$prevEE,"\t",$nextES,"\t",$nextEE,"\n";

  }


  ## retrieve inclusion level and junction count info
  my $inc_level = $cols[33];
  my $IJC        = $cols[25];
  my $SJC        = $cols[26];

  #get lengths of exons
  my $inc_len  = $cols[29];
  my $skip_len = $cols[30];

  ## only look at strong changes,  tumor junction reads > 10 reads
  next unless ($IJC >=10);
  next unless ($SJC >=10);

  ## create unique ID for splicing change
  #print $Start,"\t",$End,"\t",$prevES,"\t",$prevEE,"\t",$nextES,"\t",$nextEE,"\n";

  my $splice_id = $gene.":".$Start."-".$End."_".$prevES."-".$prevEE."_".$nextES."-".$nextEE;
  my $hist = $bs_id_hist{$bs_id};

  ##remove .0 in coord
  $splice_id=~s/\.0//g;

  $inc_levels{$splice_id}{$bs_id} = $inc_level;
  $chr{$splice_id} = $chr;
  $str{$splice_id} = $str;


  ## filter for HGG/histology of interest
  #next unless $hist=~/HGAT/;
  #next unless $cns_regions{$bs_id}=~/Midline/;

  push @splicing_events, $splice_id;
  push @{$splicing_psi{$splice_id}}, $inc_level;

  $splice_totals_per_sample{$bs_id}++;
  $splice_totals{$splice_id}++;

  $splice_event_hist{$splice_id}{$hist}=$hist;
  #print $splice_id,"\t",$hist,"\n";
}
close(FIL);

print "## make sure all arrays are unique...\n";
## make  all arrays are unique
my @broad_hist_uniq      = do { my %seen; grep { !$seen{$_}++ } @broad_hist };
my @bs_ids_uniq          = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };
my (%std_dev_psi,%count_psi,%mean_psi);

print "## calculate mean and sd for each splicing event...\n";
## calculate mean and sd for each splicing event

foreach my $splice_event(@splicing_events_uniq)
{


  my $mean     = mean  (@{$splicing_psi{$splice_event}});
  my $count    = count (@{$splicing_psi{$splice_event}});
  my $std_dev  = stddev(@{$splicing_psi{$splice_event}});

  $std_dev_psi{$splice_event} = $std_dev;
  $count_psi{$splice_event}   = $count;
  $mean_psi{$splice_event}    = $mean;
  #print $splice_event,"\t",$mean,"\t",$count,"\t",$std_dev,"\n";
}

## assess each tumor samples to identify if it is aberrant (2 standard deviations from the mean)
open(EVENTS,">results/splice_events.diff.".$splice_case.".txt");
print EVENTS "Splice ID\tCase\tSample\tHistology\tCNS\tType\n";
open(BEDPOS, ">results/splicing_events.SE.total.pos.bed");
open(BEDNEG, ">results/splicing_events.SE.total.neg.bed");

foreach my $sample(@bs_ids_uniq)
{
  #next unless $bs_id_hist{$sample}=~/HGAT/;
  foreach my $splice_event(@splicing_events_uniq)
  {
    ## look at splicing events that are present in sample and are recurrent
    next unless $inc_levels{$splice_event}{$sample};

    next if $count_psi{$splice_event} < 2;
    next unless $mean_psi{$splice_event};

    ## report only histology of interest if appropriate
    #next unless $bs_id_hist{$sample}=~/HGAT/;
    #print "second:",$bs_id_hist{$sample},"\n";
    #next unless $cns_regions{$sample}=~/Midline/;

    my $psi_tumor = $inc_levels{$splice_event}{$sample};
    my $std_psi   = $std_dev_psi{$splice_event};
    my $mean_psi  = $mean_psi{$splice_event};

    #> +2 z-scores
    #print $sample,"\t",$splice_event,"\tPSI:",$psi_tumor,"\tMeanPSI:",$mean_psi,"\tSTDPSI:",$std_psi,"\n";

    if($psi_tumor > ($mean_psi + ($std_psi + $std_psi)) )
    {

      print EVENTS $splice_event,"\t".$splice_case,"\t",$sample,"\t",$bs_id_hist{$sample},"\t",$cns_regions{$sample},"\tInclusion\n";
      print BEDPOS $chr{$splice_event},"\t";
      ## get exon coords for bed to intersect w Uniprot later
      if($splice_event=~/(\w+)\:(\d+\-\d+)/)
      {
        $exon_coord = $2;
        my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
        print BEDPOS $start_exon_coord,"\t",$end_exon_coord,"\t",$splice_event,"\t",$mean_psi,"\t",$str{$splice_event},"\n";

      }
    }
    # < -2 z-scores
    if($psi_tumor < ($mean_psi - ($std_psi + $std_psi)) )
    {
      #print $splice_event,"\tPSI:",$psi_tumor,"\tMeanPSI:",$mean_psi,"\tSTDPSI:",$std_psi,"\n";

      print EVENTS $splice_event,"\t".$splice_case,"\t",$sample,"\t",$bs_id_hist{$sample},"\t",$cns_regions{$sample},"\tSkipping\n";
      print BEDNEG $chr{$splice_event},"\t";
      ## get exon coords for bed to intersect w Uniprot later
      if($splice_event=~/(\w+)\:(\d+\-\d+)/)
      {
        $exon_coord = $2;
        my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
        print BEDNEG $start_exon_coord,"\t",$end_exon_coord,"\t",$splice_event,"\t",$mean_psi,"\t",$str{$splice_event},"\n";
      }
    }
  }
}
close(EVENTS);
close(BEDPOS);
close(BEDNEG);

## function to calculate stdev given an array of dPSI values
sub stdev{
        my($data) = @_;
        if(@$data == 1){

                return "1";
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

## function to calculate avg given an array of dPSI values
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
