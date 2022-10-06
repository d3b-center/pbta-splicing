#!/usr/bin/perl

use Statistics::Lite qw(:all);
############################################################################################################
# generate_splicing_index_tab_using_tumors.pl
#
# ./generate_splicing_index_tab_using_tumors.pl ../psi_clustering/input/pbta-histologies.RNA-Seq.initial.tsv
#                                   ~/Desktop/AS-DMG/analyses/merge_rMATS/merge_rMATS_splicing.SE.single.tsv
############################################################################################################
my ($histology,$rmats_tsv) = ($ARGV[0], $ARGV[1]);
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);
my @splicing_events;
my %splicing_psi;

# annotate histology file #
  # hash with BS and disease
  # make arrays of each histology of BS IDs
  open(FIL, $histology) || die("Cannot Open File $histology");
  while(<FIL>)
  {
    chomp;

    ## only select RNA-seq and PBTA samples
    next unless ($_=~/RNA-Seq/);
    next unless ($_=~/PBTA/);

    #filter out/skip those not solid tissue and second maligancy
    next unless ($_=~/BS/);
    next unless ($_=~/Tissue/);
    next if     ($_=~/Second\sMalignancy/);

    my @cols       = split "\t";
    my $hist       = $cols[44];
    my $bs_id      = $cols[0];
    my $patient_id = $cols[3];
    my $CNS_region = $cols[32];

    ## filter histologies of interests
    next unless ( ($hist=~/HGAT/)  ||
                  ($hist=~/LGAT/)  ||
                  ($hist=~/Oligodendroglioma/) ||
                  ($hist=~/Medulloblastoma/)   ||
                  ($hist=~/Ganglioglioma/)  ||
                  ($hist=~/Ependymoma/)||
                  ($hist=~/ATRT/)  ||
                  ($hist=~/Craniopharyngioma/) );

    ## store bs ids and histologies
    push @broad_hist, $hist;
    push @bs_ids, $bs_id;

    $bs_id_hist{$bs_id}  = $hist;
    $cns_regions{$bs_id} = $CNS_region;

    ## histology counter for downstream analysis
    $hist_count{$hist}++;

    ## hash to keep track of histology and bs id
    push @{$histology_ids{$hist}}, $bs_id;
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

my (%splice_totals_per_sample, %splice_totals);

## process rMATS output (may take awhile) from merged input file
print "processing rMATs results...\n";
open(FIL, "gunzip -c $rmats_tsv |") || die ("can’t open $rmats_tsv");
while(<FIL>)
{
  chomp;

  ##only look at exon splicing events
  next unless($_=~/^SE/);

  my @cols  = split "\t";
  my $bs_id = $cols[1];
  my $ctrl  = $cols[2];

  ## get gene name
  my $gene         = $cols[4];
  $gene=~s/\"//g; # remove quotation marks

  ## retrieve exon coordinates
  my $chr          = $cols[5];
  my $str          = $cols[6];
  my $exonStart    = $cols[7]+1; ## its 0-based so add 1
  my $exonEnd      = $cols[8];
  my $upstreamES   = $cols[15];
  my $upstreamEE   = $cols[16];
  my $downstreamES = $cols[17];
  my $downstreamEE = $cols[18];

  ## retrieve inclusion level and junction count info
  my $inc_level = $cols[33];
  my $IJC        = $cols[25];
  my $SJC        = $cols[26];

  #get lengths of exons
  my $inc_len  = $cols[29];
  my $skip_len = $cols[30];
  my $thr_diff = $cols[-1];

  ## only look at strong changes,  tumor junction reads > 10 reads
  next unless ($IJC >=10);
  next unless ($SJC >=10);

  ## create unique ID for splicing change
  my $splice_id = $gene.":".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;
  $inc_levels{$splice_id}{$bs_id} = $inc_level;

  #print $splice_id,"\t",$inc_level,"\n";
  push @splicing_events, $splice_id;
  push @{$splicing_psi{$splice_id}}, $inc_level;

  $splice_totals_per_sample{$bs_id}++;
  $splice_totals{$splice_id}++;
}
close(FIL);

print "## make sure all arrays are unique...\n";
## make  all arrays are unique
my @broad_hist_uniq      = do { my %seen; grep { !$seen{$_}++ } @broad_hist };
my @bs_ids_uniq          = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

print "## calculate mean and sd for each splicing event...\n";
## calculate mean and sd for each splicing event
my (%std_dev_psi, %count_psi, %mean_psi);
foreach my $splice_event(@splicing_events_uniq)
{
  my $mean     = mean  (@{$splicing_psi{$splice_event}});
  my $count    = count (@{$splicing_psi{$splice_event}});
  my $std_dev  = stddev(@{$splicing_psi{$splice_event}});

  $std_dev_psi{$splice_event} = $std_dev;
  $count_psi{$splice_event}   = $count;
  $mean_psi{$splice_event}    = $mean;
}

#-Look at each tumor --> is it 2 standard deviations from the mean?
#--Yes --> aberrant
#--No  --> non-aberrant
my %absplice_totals_per_sample;
my %absplice_totals_per_sample_pos;
my %absplice_totals_per_sample_neg;
foreach my $sample(@bs_ids_uniq)
{
  foreach my $splice_event(@splicing_events_uniq)
  {
    next unless $inc_levels{$splice_event}{$sample};
    my $psi_tumor = $inc_levels{$splice_event}{$sample};
    my $std_psi   = $std_dev_psi{$splice_event};
    my $mean_psi   = $mean_psi{$splice_event};
    #> +2 z-scores
    if($psi_tumor > ($mean_psi + ($std_psi + $std_psi)) )
    {
      $absplice_totals_per_sample_pos{$sample}++;
    }
    # < -2 z-scores
    if($psi_tumor < ($mean_psi - ($std_psi + $std_psi)) )
    {
      $absplice_totals_per_sample_neg{$sample}++;
    }
  }
}

#make table for plotting of splice_index
if (!-d "results")
{
  mkdir "results";
}

print "make table for plotting of splice_index...\n";
open(TAB,">results/splicing_index.total.txt");
print TAB "Sample\tTotal\tAS_neg\tAS_pos\tAS_total\tSI\tHistology\n";
foreach my $sample(@bs_ids_uniq)
{
  next unless $splice_totals_per_sample{$sample};

  print TAB $sample,"\t";
  print TAB $splice_totals_per_sample{$sample},"\t";
  print TAB $absplice_totals_per_sample_neg{$sample},"\t",$absplice_totals_per_sample_pos{$sample},"\t";

  my $total_absplice_totals_per_sample = $absplice_totals_per_sample_neg{$sample}+$absplice_totals_per_sample_pos{$sample};
  my $splice_index = $total_absplice_totals_per_sample/$splice_totals_per_sample{$sample};

  print TAB $total_absplice_totals_per_sample,"\t";
  print TAB $splice_index,"\t";
  print TAB $bs_id_hist{$sample},"\n";
}
close(TAB);
