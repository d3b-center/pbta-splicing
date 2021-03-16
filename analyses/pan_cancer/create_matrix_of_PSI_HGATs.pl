#!/usr/bin/perl

#
# create_matrix_of_PSI.pl
#
# ./create_matrix_of_PSI.pl pbta-histologies.tsv filtered_samples_files.txt

my ($histology,$rmats) = ($ARGV[0], $ARGV[1]);
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);

my (%status, %file_paths_controlBR, %file_paths_controlBR);

##add only dmg/hgat $samples for downstream analyses (normal vs tumor)
my $hist = "../../data/pbta-histologies.addedv16.dat";
open(FIL,$hist) || die("Cannot Open File $hist");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $id = $cols[0];

  my @cols = split "\t";
  my $id = $cols[0];
  my $molec_subtype = $cols[-1];

  if($_=~/Diffuse\sintrinsic\spontine\sglioma/)
  {
    if($molec_subtype=~/K28/)
    {
      $status{$id} = "H3K28";
    }
    else{
      $status{$id} = "DMG_WT";
    }
  }
  else{

    if($molec_subtype=~/K28/){
      $status{$id} = "H3K28";
    }
    elsif($molec_subtype=~/G35/){
      $status{$id} = "G35";
    }
    else
    {
      #print "*id:".$id,"\n";
      $status{$id} = "H3WT_HGAT";
    }
  }
}

##make sure to retreive correct comparison rmats file (dmg/BS and hgg/BR)
my (%file_paths_controlBS,%file_paths_controlBR);
open(FIL,"../../data/controlBS_vs_BS.SE.rmat_files.txt") || die("Cannot Open File");
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
  $file_paths_controlBS{$sample} = $file;

}
close(FIL);

open(FIL,"../../data/controlBR_vs_BS.SE.rmat_files.txt") || die("Cannot Open File");
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
  $file_paths_controlBR{$sample} = $file;

}
close(FIL);


# annotate histology file #
  # hash with BS and disease
  # make arrays of each histology of BS IDs

open(FIL,"../../data/pbta-histologies.RNA-Seq.tsv") || die("Cannot Open File $histology");
while(<FIL>)
{
  chomp;
  next if ($_=~/Kids/);
  my @cols = split "\t";
  my $broad_hist = $cols[35];
  my $bs_id = $cols[0];

  ##filter for specific histologies
  next unless $broad_hist=~/HGAT/;

  #print $broad_hist,"\n";

  push @broad_hist, $broad_hist;
  push @bs_ids, $bs_id;

  $bs_id_hist{$bs_id} = $broad_hist;
  $hist_count{$broad_hist}++;

  print "bs_id:".$bs_id,"\t", "broad_hist:",$broad_hist,"\n";
  push @{$histology_ids{$broad_hist}}, $bs_id;
}
close(FIL);

# store each library (using file name)
# for each line, make unique splice id
  # gene + chr + SE + UEx + DEx   (make array)
  # store inclusion level of SE   (hash of splice ids and library)
open(INP,"../../data/controlBS_vs_BS.SE.rmat_files.txt") || die("Cannot Open File");
while(<INP>)
{
  chomp;
  my $file = $_;
  my $sample = "";
  if($file=~/vs\_(BS\_\w+)\./)
  {
    #print "samp*".$1,"*\n";
    $sample = $1;
  }

  ## find rmat file for processing
  if( ($status{$sample}=~/H3K28/) || ($status{$sample}=~/DMG_WT/) )
  {
    $file = $file_paths_controlBS{$sample};
  }
  elsif( $status{$sample}=~/H3WT_HGAT/ )
  {
    $file = $file_paths_controlBR{$sample};
  }

  #print "file:".$file,"\n";
  my $bs_id = $sample;

  open(SE_FIL,$file) || die("Cannot Open File $file");
  while(<SE_FIL>)
  {
    #print $_,"HELL\n";
    my @cols = split "\t";
    my $gene         = $cols[2];
    my $chr          = $cols[3];
    my $exonStart    = $cols[5]+1; ## bc it 0-based
    my $exonEnd      = $cols[6];
    my $upstreamES   = $cols[7];
    my $upstreamEE   = $cols[8];
    my $downstreamES = $cols[9];
    my $downstreamEE = $cols[10];
    my $inc_level    = $cols[20];

    my $tumor_IJC     = $cols[14];
    my $tumor_SJC     = $cols[15];
    my $inc_from_len  = $cols[16];
    my $skip_from_len = $cols[17];
    my $pval = $cols[18];
    my $thr_diff = $cols[-1];

    #next unless ($inc_level >=.10);



    ## only look at strong changes
    next unless (abs($thr_diff) >= .10);
    next unless ($pval<=0.05);
    next unless ($tumor_IJC >=10);
    next unless ($tumor_SJC >=10);

    #$inc_level    = $cols[12];

    $gene=~s/\"//g;

    #print $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE,"\t",$inc_level,"*\n";
    my $splice_id= $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;

    ##store inclusion level differences between healthy vs tumor
    $thr_diff=~s/\s+//; $thr_diff=~s/\n+//;

    $inc_levels{$splice_id}{$bs_id} = $thr_diff;

    my $hist_of_sample = $bs_id_hist{$bs_id};

    ##store number of events per histology
    #print $hist_of_sample,":hist\n";
    $hist_check{$splice_id}{$hist_of_sample}++;
    push @splicing_events, $splice_id;

  }
  close(SE_FIL);
}

close(INP);

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

print "Splice_ID";
foreach my $hist (sort @broad_hist_uniq)
{
  foreach my $sample (@{$histology_ids{$hist}})
  {
    print "\t";

    # print histology name (renamed)
    # print $hist."_".$i;

    # print sample name
    print $sample;
    $i++;

  }
}
#print "\n";

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
        $thresh_for_hist_prev = 0;
      }
  }
  next if ($thresh_for_hist_prev == 0);
  print "\n",$event;
  foreach my $hist (sort @broad_hist_uniq)
  {
    foreach my $sample (@{$histology_ids{$hist}})
    {
      print "\t";
      #print $sample,":",$hist,":";
      my $dPSI = 0;
      if($inc_levels{$event}{$sample}) { $dPSI = $inc_levels{$event}{$sample} };

      print $dPSI;
    #  if(abs($inc_levels{$event}{$sample})>0)
    #  {
        #print $sample.":".$inc_levels{$event}{$sample};
        #print $inc_levels{$event}{$sample};
    #  }
    #  else
    #  {
        # print $sample.":0";
    #    print "0";
    #  }

    }
  }
  #print "\n";
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
