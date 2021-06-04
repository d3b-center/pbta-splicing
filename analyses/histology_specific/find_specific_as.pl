#!/usr/bin/perl

my ($histology,$rmats) = ($ARGV[0], $ARGV[1]); # ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt
my (@broad_hist, @bs_id, @splicing_events);

my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);

my @ctrl_files =</Users/naqvia/Desktop/pan_cancer_rmats/*single_normals.SE.MATS.JC.txt>;

my %ctrl_inc_levels;
my %cns_regions;
my %cluster_group;

# annotate histology file #
  # hash with BS and disease
  # make arrays of each histology of BS IDs
open(FIL, $histology) || die("Cannot Open File $histology");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $broad_hist = $cols[35];
  my $bs_id      = $cols[0];
  my $CNS_region = $cols[1];
  my $cluster = $cols[-1];

  next if $_=~/Kids_First_Biospecimen_ID/;

  ##filter for specific histologies
  #print $broad_hist,"\n";
  next if $broad_hist =~/Embryonal/;
  next if $broad_hist =~/EWS/;
  next if $broad_hist =~/Oligodendroglioma/;

  push @broad_hist, $broad_hist;
  push @bs_ids, $bs_id;

  $bs_id_hist{$bs_id} = $broad_hist;
  $cns_regions{$bs_id} = $CNS_region;
  $cluster_group{$bs_id} = $cluster;

  #print "cns: ",$CNS_region,"\n";
  $hist_count{$broad_hist}++;

  #print "bs_id:".$bs_id,"\n";
  push @{$histology_ids{$broad_hist}}, $bs_id;
}
close(FIL);

## annotate normals/ctrl files
my %ctrl_path;

open(FIL,"input/ctrl_names.txt") || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my ($ctrl_name, $file_path) = split;
  #print "*".$file_path,"*\t*",$ctrl_name,"*\n";

  next if($ctrl_name=~/control-FB/);
  next if($ctrl_name=~/SRR4787052/);

  #print $file_path,"\t",$ctrl_name,"\n";
  $ctrl_path{$file_path} = $ctrl_name;
}
close(FIL);


## store "healthy" psi values/gene
my %ctrl_psi;
my %ctrl_event_filter ;
my $ctrl_splice_totals = 0;
my %ctrl_inc_levels;

foreach my $file(@ctrl_files)
{
  my @fields = split/\//,$file;
  my $file_name = $fields[-1];

  next unless $ctrl_path{$file_name};

  open(CTRL_SE,$file) || die("Cannot Open File $file");
  while(<CTRL_SE>)
  {
    chomp;
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

    my $ctrl_IJC     = $cols[12];
    my $ctrl_SJC     = $cols[13];

    my $inc_from_len  = $cols[16];
    my $skip_from_len = $cols[17];

    my $pval = $cols[18];
    my $thr_diff = $cols[-1];

    #next unless ($inc_level >=.10);

    ## only look thsoe with >=10 reads
    next unless ($ctrl_IJC >=10);
    next unless ($ctrl_SJC >=10);

    ## create unique ID for splicing change
    $gene=~s/\"//g;
    my $splice_id= $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;

    #print "ctrl_splice_id: ".$spliced_id,"\n";
    push @{$ctrl_inc_levels{$splice_id}}, $inc_level;

    ##store control
    my $ctrl_name = $ctrl_path{$file_name};
    #print "ctrl_name :",$ctrl_name,"\n";
    $ctrl_inc_levels{$ctrl_name}{$splice_id} = $inc_level;

    push @splicing_events, $splice_id;
  #  print "splicing event CTRL: ".$splice_id,"\n";
    $ctrl_event_filter{$splice_id}=1;
    $ctrl_splice_totals++;

  }
  close(CTRL_SE);
}

# store each library (using file name)
# for each line, make unique splice id
  # gene + chr + SE + UEx + DEx   (make array)
  # store total splicing events and those that pass threshold for each sample
my (%ctrl_splice_totals_per_sample);

## splicing event freq per histology hash
my %splicing_event_freq_hist;

open(FIL,$rmats) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my $file = $_;
  my $sample = "";

  my $bs_id = "";
  if($file=~/\.(BS\_\w+)\./)
  {
    $bs_id = $1;
  }

  open(SE_FIL,$file) || die("Cannot Open File $file");
  while(<SE_FIL>)
  {
    chomp;
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

    ##skip if not in healthy samples
    next unless $ctrl_event_filter{$splice_id};

    ## compute avg psi in healthy
    my $total_ctrl_psi_per_event = 0;
    my $ctrl_count = 0;
    #print "ctrl_psi_array: ",@{$ctrl_inc_levels{$splice_id}},"\n";

    foreach my $psi (@{$ctrl_inc_levels{$splice_id}})
    {
      #print "ctrl_psi:",$psi,"\n";
      $total_ctrl_psi_per_event = $total_ctrl_psi_per_event + $psi;
      $ctrl_count++;
    }

    my $avg_ctrl_inc = ($total_ctrl_psi_per_event / $ctrl_count);
    my $dpsi = 0;

    ## compare with approp controls, if not, compare to average of all ctrls
    if( ($bs_id_hist{$bs_id}=~/HGAT/) && ($cns_regions{$bs_id} =~/Midline/) )
    {
      $dpsi = $ctrl_inc_levels{"control-BS"}{$splice_id} - $inc_level;
      #print $_,"***",$dpsi,":",$ctrl_inc_levels{"control-BS"}{$splice_id},":",$inc_level,"\n";
    }
    elsif($bs_id_hist{$bs_id}=~/Medu/)
    {
      $dpsi = $ctrl_inc_levels{"control-BC"}{$splice_id} - $inc_level;
      print "Medullo_sample\t",$dpsi,"\t",$splice_id,"\n";
    }
    else{
      $dpsi = $avg_ctrl_inc-$inc_level;
    }

    #print $splice_id,"\tdpsi: ",$dpsi,"\n";
    $ctrl_splice_totals_per_sample{$bs_id}++;
    if(abs($dpsi)>=.20)
    {
      #print "ctrl avg: ",$avg_ctrl_inc,"\t",$inc_level,"\t";

      #print "TRUE\t";
      $splice_filtered_totals_per_sample{$bs_id}++;
      #print $splice_id,"\tdpsi: ",$dpsi,"\t",$bs_id,"\n";
      $inc_levels{$splice_id}{$bs_id} = $dpsi;

      ## store frequency in histology type
      my $hist = $bs_id_hist{$bs_id};
      #print $hist,"\n";
      $splicing_event_freq_hist{$hist}{$splice_id}++;

    }
    #my $hist_of_sample = $bs_id_hist{$bs_id};
    #$hist_check{$splice_id}{$hist_of_sample}++;
    #push @splicing_events, $splice_id;

  }
  close(SE_FIL);
}
close(FIL);

## make sure all arrays are unique
my @broad_hist_uniq      = do { my %seen; grep { !$seen{$_}++ } @broad_hist };
my @bs_ids_uniq          = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

my %splice_event_perc_hist;

## compute percentages of splicing event in each hist
foreach my $hist(@broad_hist_uniq)
{
  foreach my $event (@splicing_events_uniq)
  {
    next unless $splicing_event_freq_hist{$hist}{$event};
    my $perc_in_hist = ($splicing_event_freq_hist{$hist}{$event}/$hist_count{$hist}) * 100;
    if($perc_in_hist >= 30)
    {
      $splice_event_perc_hist{$hist}{$event} = $perc_in_hist;
    }
  }
  #$splicing_event_freq_hist{$hist}{$splice_id}++;
}

foreach my $hist(@broad_hist_uniq)
{
    foreach my $event (@splicing_events_uniq)
    {

      next unless $splice_event_perc_hist{$hist}{$event};

      print $hist,"\t";
      print $event,"\n";
    }
}

__DATA__
##make table for plotting of splice_index
print "sample\ttotal_splice_events\tfiltered_splice_events\tsplice_index\thist\n";
foreach my $sample(@bs_ids_uniq)
{
  #next unless $splice_filtered_totals_per_sample{$sample} > 0;
  #if($splice_filtered_totals_per_sample{$sample} < 1){ print "flag\t", $sample };
  next unless $splice_filtered_totals_per_sample{$sample} > 0;

  print $sample,"\t";
  print $ctrl_splice_totals_per_sample{$sample},"\t";
  print $splice_filtered_totals_per_sample{$sample},"\t";

  my $splice_index = $splice_filtered_totals_per_sample{$sample}/$ctrl_splice_totals_per_sample{$sample};
  print $splice_index,"\t";
  print $bs_id_hist{$sample},"\n";

}

perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep ATRT > results/perc_hist_as.30prev.ATRT.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep Craniopharyn > results/perc_hist_as.30prev.Cran.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep Epend > results/perc_hist_as.30prev.Epend.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep -i gang  > results/perc_hist_as.30prev.Gang.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep -i hgat  > results/perc_hist_as.30prev.HGAT.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep -i lgat  > results/perc_hist_as.30prev.LGAT.tsv
perl find_specific_as.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt | grep -i Medu  > results/perc_hist_as.30prev.medul.tsv
