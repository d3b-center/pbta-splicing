#!/usr/bin/perl
################################################################################
# create_matrix_of_PSI.pl
# create tsv file of PSI for each sample, to be used for consensus clustering
# written by Ammar Naqvi
#
# usage: ./create_matrix_of_PSI.pl <histology file> <rMATs_file_paths.txt>
################################################################################
my ($histology,$rmats, $abridged) = ($ARGV[0], $ARGV[1], $ARGV[2]);
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);

my (%init_patientIDs, %progr_patientIDs, %recur_patientIDs, %unavail_patientIDs);
my %patient_id_count;

open(FIL,$abridged) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $patient_id = $cols[1];
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
  elsif($_=~/Recu/)
  {
    $recur_patientIDs{$patient_id} = $bs_id;
  }
  elsif($_=~/Unavail/)
  {
    $unavail_patientIDs{$patient_id} = $bs_id;
  }
  else{}

}
# annotate histology file #
  # hash with BS and disease
  # make arrays of each histology of BS IDs

open(FIL, $histology) || die("Cannot Open File $histology");
while(<FIL>)
{
  chomp;
  next if ($_=~/Kids/);
  my @cols = split "\t";
  my $broad_hist = $cols[35];
  my $bs_id = $cols[0];
  my $patient_id = $cols[4];

  ##filter for specific histologies
  next if($broad_hist=~/EWS/);
  next if($broad_hist=~/Embryonal/);

  if($patient_id_count{$patient_id} > 1)
  {
    #print $bs_id,"\tdouble",$patient_id,"\t";
    #print $broad_hist,"\n";

    if($init_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $init_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      print $bs_id,"\t",$patient_id,"\t";
      print $broad_hist,"\n";
    }
    elsif($progr_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $progr_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      print $bs_id,"\t",$patient_id,"\t";
      print $broad_hist,"\n";
    }
    elsif($recur_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $recur_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      print $bs_id,"\t",$patient_id,"\t";
      print $broad_hist,"\n";
    }
    elsif($unavail_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $unavail_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      print $bs_id,"\t",$patient_id,"\t";
      print $broad_hist,"\n";
    }
    else{ print $_,"***\n"; }
  }



  push @broad_hist, $broad_hist;
  push @bs_ids, $bs_id;

  $bs_id_hist{$bs_id} = $broad_hist;
  $hist_count{$broad_hist}++;

  #print "bs_id:".$bs_id,"\n";
  push @{$histology_ids{$broad_hist}}, $bs_id;
}
close(FIL);
#die;

# store each library (using file name)
# for each line, make unique splice id
  # gene + chr + SE + UEx + DEx   (make array)
  # store inclusion level of SE   (hash of splice ids and library)

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

    next unless ($inc_level >=.10);
    #$inc_level    = $cols[12];

    $gene=~s/\"//g;

    #print $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE,"\t",$inc_level,"*\n";
    my $splice_id= $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;
    $inc_levels{$splice_id}{$bs_id} = $inc_level;
    my $hist_of_sample = $bs_id_hist{$bs_id};

    ##store number of events per histology
    #print $hist_of_sample,":hist\n";
    $hist_check{$splice_id}{$hist_of_sample}++;
    push @splicing_events, $splice_id;

  }
  close(SE_FIL);
}
close(FIL);

my @broad_hist_uniq      = do { my %seen; grep { !$seen{$_}++ } @broad_hist };
my @bs_ids_uniq          = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };


# go through each splicing event
  # order by disease library
  # print splice id, inc level

#   foreach my $hist (sort @broad_hist_uniq)
#   {
#     my $i = 1;
#     foreach my $sample (@{$histology_ids{$hist}})
#     {
#       #last if($i>10);
#       print $hist,"\n";
#       $i++;
#
#     }
#   }
#   print "\n";
# die;


my $out_file = "results/pan_cancer_splicing.thr10.report_select.remDup.txt";
open(OUT,">",$out_file) || die("Cannot Open File");
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





foreach my $event (@splicing_events_uniq)
{
  ## @{$histology_ids{$broad_hist}}
  my $thresh_for_hist_prev = 1;
  #print "event:",$event,"\n";
  foreach my $hist (sort @broad_hist_uniq){

      ##threshold for prevalence per sample (%)
      #my $hist_prev_thr = .20*($hist_count{$hist});
      my $hist_prev_thr = 2;
      my $prev_in_hist = 0;
      if($hist_check{$event}{$hist}) { $prev_in_hist = $hist_check{$event}{$hist} };
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
      if($inc_levels{$event}{$sample}>0)
      {
        #print $sample.":".$inc_levels{$event}{$sample};
        print OUT $inc_levels{$event}{$sample};
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
