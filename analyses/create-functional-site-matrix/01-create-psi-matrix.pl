#!/usr/bin/perl

################################################################################
# create-psi-matrix-primary.pl
# create matrix file of PSI for each sample to be used for multimodal clustering
# written by Ammar Naqvi
#
# usage: ./create_matrix_of_PSI_removeDups.pl <rMATs_file_paths.txt>
#                                             <splice type> <outfile>
################################################################################

## get args-- input files and output file
my ($rmats_tsv, $splice_type ,$out_file) = @ARGV[0, 1, 2];

## check to see if splicing case is correct
unless ($splice_type=~/SE$|A3SS$|A5SS$|RI$/)
{
  die("Splicing case does not exist, please use either 'SE', 'A5SS', 'A3SS', or 'RI'");
}

## data structures
my (@broad_hist, @bs_id, @splicing_events);
my (%str, %chr);

## process rMATS output (may take awhile) from merged input file
open(FIL, "gunzip -c $rmats_tsv |") || die ("can’t open $rmats_tsv");
while (<FIL>) {
    chomp;

    next unless($_=~/^$splice_type/);

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

    if($splice_type=~/SE/)
    {
      $Start    = $cols[7]; ## its 0-based so add 1
      $End      = $cols[8];
      $prevES   = $cols[15];
      $prevEE   = $cols[16];
      $nextES = $cols[17];
      $nextEE = $cols[18];
    }
    elsif($splice_type=~/RI/)
    {
      $Start    = $cols[13]+1; ## its 0-based so add 1
      $End      = $cols[14];
      $prevES   = $cols[15];
      $prevEE   = $cols[16];
      $nextES = $cols[17];
      $nextEE = $cols[18];
    }
    elsif($splice_type=~/A5SS/)
    {
      $Start    = $cols[19]; ## its 0-based so add 1
      $End      = $cols[20];
      $prevES   = $cols[21];
      $prevEE   = $cols[22];
      $nextES = $cols[23];
      $nextEE = $cols[24];

    }
    else{
      $Start    = $cols[19]; ## its 0-based so add 1
      $End      = $cols[20];
      $prevES   = $cols[21];
      $prevEE   = $cols[22];
      $nextES = $cols[23];
      $nextEE = $cols[24];
    }

    ## retrieve inclusion level and junction count info
    my $inc_level  = $cols[33];
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

    my $splice_id = $gene.":".$Start."-".$End."_".$prevES."-".$prevEE."_".$nextES."-".$nextEE;
    $splice_id=~s/\.0//g;

    $str{$splice_id} = $str;
    $chr{$splice_id} = $chr;

    $inc_levels{$splice_id}{$bs_id} = $inc_level;

    ## get gene name
    my $gene = $cols[4];
    $gene =~ s/\"//g; # remove quotation marks

    next unless (($IJC >=10 ) &&
                 ($SJC >=10 ));


    ## store inclusion levels of splicing event and BS ID#
    $inc_levels{$splice_id}{$bs_id} = $inc_level;

    ## store array of inclusion levels per gene to access later
    push @genes, $gene;
    push @splicing_events, $splice_id;
    push @bs_ids, $bs_id;

}
close(FIL);

my @bs_ids_uniq = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };
my @genes_uniq = do { my %seen; grep { !$seen{$_}++ } @genes };


open(OUT, ">", $out_file) || die("Cannot Open File");
open(OUTBED, ">", $out_file.".bed") || die("Cannot Open File");

# go through each splicing event
print OUT "Splice_ID";
foreach my $sample (sort @bs_ids_uniq) {
  print OUT "\t";
  print OUT $sample;
}
print OUT "\n";

## print each event and threshold value to output file and in bed file
foreach my $event (@splicing_events_uniq) {
  print OUT $event;
  foreach my $sample (@bs_ids_uniq) {
      print OUT "\t";
      my $psi = $inc_levels{$event}{$sample};
      if ($psi > 0) {
          print OUT $psi;
      } else {
          print OUT "0";
      }
  }

  if($event=~/(\w+)\:(\d+\-\d+)\_\d+\-\d+\_\d+\-\d+/)
  {
    $exon_coord = $2;
    my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
    print OUTBED $chr{$event},"\t",$start_exon_coord,"\t",$end_exon_coord,"\t",$event,"\t1\t",$str{$event},"\n";
  }

  print OUT "\n";
}

close(OUT);
close(OUTBED);
