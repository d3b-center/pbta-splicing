#!/usr/bin/perl

#
# generate_splicing_index_tab.pl
#
# ./generate_splicing_index_tab.pl pbta-histologies.tsv filtered_samples_files.v2.txt

my ($histology,$rmats) = ($ARGV[0], $ARGV[1]);
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);


# annotate histology file #
  # hash with BS and disease
  # make arrays of each histology of BS IDs
open(FIL, $histology) || die("Cannot Open File $histology");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $broad_hist = $cols[35];
  my $bs_id = $cols[0];

  ##filter for specific histologies
  #print $broad_hist,"\n";
  push @broad_hist, $broad_hist;
  push @bs_ids, $bs_id;

  $bs_id_hist{$bs_id} = $broad_hist;
  $hist_count{$broad_hist}++;

  #print "bs_id:".$bs_id,"\n";
  push @{$histology_ids{$broad_hist}}, $bs_id;
}
close(FIL);

# store each library (using file name)
# for each line, make unique splice id
  # gene + chr + SE + UEx + DEx   (make array)
  # store total splicing events and those that pass threshold for each sample

open(FIL,$rmats) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my $file = $_;
  my $sample = "";
  if($file=~/vs\_(BS\_\w+)\./)
  {
    #print "samp*".$1,"*\n";
    $sample = $1;
  }

  my $bs_id = "";
  if($file=~/\.(BS\_\w+)\./)
  {
    $bs_id = $1;

  }

  open(SE_FIL,$file) || die("Cannot Open File $file");
  while(<SE_FIL>)
  {
    #print $_,"HELL\n";
    my @cols = split "\t";
    my $gene         = $cols[2];
    my $chr          = $cols[3];
    my $exonStart    = $cols[5]+1;
    my $exonEnd      = $cols[6];
    my $upstreamES   = $cols[7];
    my $upstreamEE   = $cols[8];
    my $downstreamES = $cols[9];
    my $downstreamEE = $cols[10];
    my $inc_level    = $cols[20];

    ## create unique ID for splicing change
    $gene=~s/\"//g;
    my $splice_id= $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;
    $splice_totals_per_sample{$bs_id}++;

    next unless ($inc_level>=.10);

    $splice_filtered_totals_per_sample{$bs_id}++;
    $inc_levels{$splice_id}{$bs_id} = $inc_level;

    my $hist_of_sample = $bs_id_hist{$bs_id};
    $hist_check{$splice_id}{$hist_of_sample}++;
    push @splicing_events, $splice_id;


  }
  close(SE_FIL);
}
close(FIL);

## make sure all arrays are unique
my @broad_hist_uniq      = do { my %seen; grep { !$seen{$_}++ } @broad_hist };
my @bs_ids_uniq          = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

##make table for plotting of splice_index
print "sample\ttotal_splice_events\tfiltered_splice_events\tsplice_index\n";
foreach my $sample(@bs_ids_uniq)
{
  next unless $splice_totals_per_sample{$sample} > 0;
  print $sample,"\t";
  print $splice_totals_per_sample{$sample},"\t";
  print $splice_filtered_totals_per_sample{$sample},"\t";

  my $splice_index = $splice_filtered_totals_per_sample{$sample}/$splice_totals_per_sample{$sample};
  print $splice_index,"\t";
  print $bs_id_hist{$sample},"\n";

}

__DATA__
