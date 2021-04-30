#!/usr/bin/perl

######################################################################################
# extract_recurrent_splicing_events.pl
# descr: extract and append splicing ids from rMATS to be used to overlap with Uniprot
# usage: perl extract_recurrent_splicing_events.pl ../../data/pbta-histologies.tsv
######################################################################################

## input argument for histology file from OpenPBTA project
my $histology = $ARGV[0];

## directory of rMATS tsv files (SE type)
my @rmats_tsv = </Users/naqvia/Desktop/rmats_run/SE/*control-BS*JC.txt>;

## data structures to store values
my %dmg_samples;
my (@splicing_events, @splicing_event_deltapsis, @samples);
my %splicing_event_counts;
my %inc_levels;
my (%chr,%num_exons,%str);

print "processing histology...\n";

## open histology file and annotate DMG samples
open(FIL,$histology) || die("Cannot Open File");
while(<FIL>)
{
  chomp;

  ## skip unless RNA-seq
  next unless ($_=~/RNA\-Seq/);

  my @cols = split "\t";
  my $molecular_subtype = $cols[30];
  next unless ($molecular_subtype=~/DMG/);
  my $bs_id = $cols[0];
  $dmg_samples{$bs_id} = $bs_id;
}
close(FIL);

## process rMATS output
print "processing rMAT files...\n";
foreach my $file(@rmats_tsv)
{
  my $bs_id = "";
  if($file=~/vs\_(BS\_\w+)\./)
  {
    #print "samp*".$1,"*\n";
    $bs_id = $1;
  }
  next unless $dmg_samples{$bs_id};
  push @samples, $bs_id;

  open(FIL,$file) || die("Cannot Open File");
  while(<FIL>)
  {
    chomp;
    my @cols = split "\t";
    my $gene         = $cols[2];

    ## exon coordinates
    my $chr          = $cols[3];
    my $str          = $cols[4];
    my $exonStart    = $cols[5]+1;
    my $exonEnd      = $cols[6];
    my $upstreamES   = $cols[7];
    my $upstreamEE   = $cols[8];
    my $downstreamES = $cols[9];
    my $downstreamEE = $cols[10];

    ## inclusion level and junction count info
    my $inc_level    = $cols[20];
    my $ctrl_IJC      = $cols[12];
    my $ctrl_SJC      = $cols[13];
    my $tumor_IJC     = $cols[14];
    my $tumor_SJC     = $cols[15];

    #lengths of exons
    my $inc_len  = $cols[16];
    my $skip_len = $cols[17];

    ## p-value and dPSI values
    my $pval = $cols[18];
    my $thr_diff = $cols[-1];

    ## only look at strong changes, remove any dPSI < .10 and p-val >0.05, and IJC < 10 reads
    next unless (abs($thr_diff) >= .10);
    next unless ($pval<=0.05);
    next unless ( ($tumor_IJC >=10) || ($tumor_SJC >=10) );

    $gene=~s/\"//g;

    ## annotate and re-name splice event IDs
    #print $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE,"\t",$inc_level,"*\n";
    my $splice_id= $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;
    my $splice_id_short= $gene."_".$exonStart."-".$exonEnd;

    ##keep one splicing event per exon (choose larger dPSI), sometimes with different upstream/downstream exons
    my $related_event = 0;
    if($inc_levels{$splice_id_short}{$bs_id} > 0)
    {
      $related_event = $inc_levels{$splice_id_short}{$bs_id};
    }
    #print "splice_id_short:".$splice_id_short,"\tdPSI:",$thr_diff,"E\tRE:",$related_event,"\n";
    next unless(abs($thr_diff) > abs($related_event));
    #print "splice_id_short:".$splice_id_short,"\tdPSI:",$thr_diff,"E\tRE:",$related_event,"\n";

    $chr{$splice_id_short} = $chr;
    $str{$splice_id_short} = $str;
    $inc_levels{$splice_id_short}{$bs_id} = $thr_diff;

    ## store splicing events, inclusion levels, and splicing event counts/freqs
    push @splicing_events, $splice_id_short;
    push @{$splicing_event_deltapsis{$splice_id_short}}, $thr_diff;
    $splicing_event_counts{$splice_id_short}++;

  }
}
close(FIL);

## make arrays unique
my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

## write to files for summary tables (with a bed version for downstream analyses)
open(TAB, ">splicing_events.dpsi10.jc10.rec2.tsv");
open(BED, ">splicing_events.dpsi10.jc10.rec2.bed");

print "writing output...\n";
foreach my $event (@splicing_events_uniq)
{
  print "event: ",$event,"\t";
  print join "\t",@{$splicing_event_deltapsis{$event}},"\n";

  ## compute avg and stdevs of each lsv dPSI
  my $avg_dpsi = &average(\@{$splicing_event_deltapsis{$event}});
	my $std_dpsi = &stdev(\@{$splicing_event_deltapsis{$event}});

  my ($gene,$exon_coord) = split/\_/,$event;

  next unless $splicing_event_counts{$event} > 2;
  print TAB $gene,"\t",$event,"\t",$avg_dpsi,"\t",$std_dpsi,"\t";
  print TAB $splicing_event_counts{$event};

  print TAB "\t",$chr{$event},":",$exon_coord,"\t",$str{$event},"\n";

  print BED $chr{$event},"\t";
  if($event=~/(\w+)\_(\d+\-\d+)/)
  {
    #$lsv_coord = $1;
    $exon_coord = $2;
    #print BED $exon_coord." ";
    my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
    print BED $start_exon_coord,"\t",$end_exon_coord,"\t",$event,"\t",$avg_dpsi,"\t",$str{$event},"\n";
  }
}
close(TAB);
close(BED);

## calculate stdev given an array of PSI values
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

## calculate avg given an array of PSI values
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


__DATA__
