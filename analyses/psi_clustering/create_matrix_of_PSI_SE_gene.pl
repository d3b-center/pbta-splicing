#!/usr/bin/perl

use List::Util qw(max);

################################################################################
# create_matrix_of_PSI.pl
# create tsv file of PSI for each sample, to be used for consensus clustering
# written by Ammar Naqvi
#
# usage: ./create_matrix_of_PSI_removeDups.pl <histology file> <rMATs_file_paths.txt>
################################################################################
my ($histology,$rmats_tsv, $abridged) = ($ARGV[0], $ARGV[1], $ARGV[2]); ## args: histology, rmats path, and abridged histology
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);

my (%init_patientIDs, %progr_patientIDs, %recur_patientIDs, %unavail_patientIDs);
my %patient_id_count;

## store tumor descriptor info for filtering later
open(FIL,$histology) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  next unless ($_=~/RNA-Seq/);

  my $patient_id = $cols[3];
  my $bs_id      = $cols[0];

  #filter out/skip those not solid tissue and second maligancy
  next unless ($_=~/BS/);
  next unless ($_=~/Tissue/);
  next if     ($_=~/Second/);

  $patient_id_count{$patient_id}++;

  ##store tumor descriptor
  if($_=~/Initial/)
  {
    $init_patientIDs{$patient_id} = $bs_id;
  }
  elsif($_=~/Progressive/)
  {
    $progr_patientIDs{$patient_id} = $bs_id;
  }
  elsif($_=~/Recurrence/)
  {
    $recur_patientIDs{$patient_id} = $bs_id;
  }
  elsif($_=~/Unavail/)
  {
    $unavail_patientIDs{$patient_id} = $bs_id;
  }
  else{}

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
    my $hist       = $cols[38]; #short histology in ../psi_clustering/input/v19_plus_20210311_pnoc.tsv
    my $harm_dx    = $cols[36]; ##harmonized diagnosis
    my $bs_id      = $cols[0];
    my $patient_id = $cols[3];
    my $CNS_region = $cols[32];

    #my $cluster    = $cols[-1]; ## if

    ## skip non RNA-seq samples
    next if $_=~/Kids_First_Biospecimen_ID/;
    next unless $_=~/RNA-Seq/;

    ##filter for specific histologies
    #next if($hist=~/EWS/);
    #next if($hist=~/Embryonal/);

    next unless ( ($hist=~/HGAT/) || ($hist=~/LGAT/) || ($hist=~/Oligo/) || ($hist=~/Med/) || ($hist=~/Gang/) || ($$hist=~/Epend/) || ($hist=~/ATRT/) || ($hist=~/Crani/) );


  ## if patient has more than one RNA-seq dataset
  if($patient_id_count{$patient_id} > 1)
  {

    ## prioritize and keep only initial RNA-seq set
    if($init_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $init_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      #print $bs_id,"\t",$patient_id,"\t";
      #print $broad_hist,"\n";
    }

    ## if no initial, is it progressive? if so, prioritize/keep that only
    elsif($progr_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $progr_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      #print $bs_id,"\t",$patient_id,"\t";
      #print $broad_hist,"\n";
    }

    ## if no initial or progressive? prioritize/keep recurrence tumors  only
    elsif($recur_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $recur_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      #print $bs_id,"\t",$patient_id,"\t";
      #print $broad_hist,"\n";
    }

    ## if unknown, just keep it
    elsif($unavail_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $unavail_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      #print $bs_id,"\t",$patient_id,"\t";
      #print $broad_hist,"\n";
    }
    else{
      #print $_,"***\n";
    }
  }

  ## make an array and store histology information and BS IDs
  push @broad_hist, $hist;
  push @bs_ids, $bs_id;

  $bs_id_hist{$bs_id} = $hist;

  ## store total number of histologies
  $hist_count{$hist}++;
  push @{$histology_ids{$hist}}, $bs_id;
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
  open(FIL,$rmats_tsv) || die("Cannot Open File");
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
my $out_file = "input/pan_cancer_splicing_SE.gene.txt";
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

__DATA__
28612	"ENSG00000079308.19"	"TNS1"	chr2	-	217813214	217813304	217812367	217812445	217813684	217813816	28612	365	17			239	149	NA	NA	0.93		NA

		upstreamES	upstreamEE	downstreamES	downstreamEE
Kids_First_Biospecimen_ID	CNS_region	sample_id	aliquot_id	Kids_First_Participant_ID	experimental_strategy	sample_type	composition	tumor_descriptor	primary_site	reported_gender	race	ethnicity	age_at_diagnosis_days	pathology_diagnosis	RNA_library	OS_days	OS_status	cohort	age_last_update_days	seq_center	parent_aliquot_id	cancer_predispositions	pathology_free_text_diagnosis	cohort_participant_id	extent_of_tumor_resection	germline_sex_estimate	normal_fraction	tumor_fraction	tumor_ploidy	molecular_subtype	integrated_diagnosis	Notes	harmonized_diagnosis	broad_histology	short_histology

cat pbta-histologies.RNA-Seq.tsv | grep "Low-grade astrocytic" > pbta-histologies.RNA-Seq.initial.tsv
cat pbta-histologies.RNA-Seq.tsv | grep "Diffuse astrocytic and oligodendroglial tumor" >> pbta-histologies.RNA-Seq.initial.tsv
cat pbta-histologies.RNA-Seq.tsv | grep "Medulloblastoma" >> pbta-histologies.RNA-Seq.initial.tsv
cat pbta-histologies.RNA-Seq.tsv | grep "Ependymoma" >> pbta-histologies.RNA-Seq.initial.tsv
cat pbta-histologies.RNA-Seq.tsv | grep "Craniopharyngioma" >> pbta-histologies.RNA-Seq.initial.tsv
