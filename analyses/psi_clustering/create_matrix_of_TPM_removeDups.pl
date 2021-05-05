#!/usr/bin/perl
################################################################################
# create_matrix_of_TPM_removeDups.pl
# create tsv file of TPM for each sample, to be used for consensus clustering
# written by Ammar Naqvi
#
# usage: ./create_matrix_of_TPM_removeDups.pl <histology file> <polyA_gene_counts> <stranded_gene_counts> <abridged hist>
################################################################################
my ($histology,$polyA, $stranded, $abridged) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4]);
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);

my (%init_patientIDs, %progr_patientIDs, %recur_patientIDs, %unavail_patientIDs);
my %patient_id_count;
my @gene_id;

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

      #print $bs_id,"\t",$patient_id,"\t";
      #print $broad_hist,"\n";
    }
    elsif($progr_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $progr_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      #print $bs_id,"\t",$patient_id,"\t";
      #print $broad_hist,"\n";
    }
    elsif($recur_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $recur_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      #print $bs_id,"\t",$patient_id,"\t";
      #print $broad_hist,"\n";
    }
    elsif($unavail_patientIDs{$patient_id})
    {
      my $bs_id_to_check = $unavail_patientIDs{$patient_id};
      next unless ($bs_id=~/$bs_id_to_check/);

      #print $bs_id,"\t",$patient_id,"\t";
      #print $broad_hist,"\n";
    }
    else{ #print $_,"***\n";
    }
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
  my @samples = "";

open(FIL,$polyA) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  #my $file = $_;
  if($_=~/gene_id/)
  {
    #print "line:",$_,"\n";
    @samples = split "\t";
    shift @samples;
    shift @samples;
    #print $bs_id[0];
  }
  else{
      my @tpm = split "\t";
      #print "line:",$_,"\n";

      shift @tpm;
      my $gene = shift @tpm;
      $gene=~s/\"//g;

      my $k = 0;
      foreach my $tpm(@tpm)
      {
        my $sample = $samples[$k];
        $sample=~s/\"//g;

        #print $gene,"\t",$sample,"\t",$tpm,"\n";
        $tpm_counts{$gene}{$sample} = $tpm;
        my $hist_of_sample = $bs_id_hist{$sample};

        ##store number of events per histology
        $hist_check{$gene}{$hist_of_sample}++;
        $k++;
      }
      push @gene_id, $gene;
  }
}
close(FIL);

my @samples = "";
open(FIL,$stranded) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  #my $file = $_;
  my @bs_id = "";
  if($_=~/gene_id/)
  {
    @samples = split "\t";
    shift @samples;
    shift @samples;
  }
  else{
      my @tpm = split "\t";
      shift @tpm;
      my $gene = shift @tpm;
      $gene=~s/\"//g;

      my $k = 0;
      foreach my $tpm(@tpm)
      {
        my $sample = $samples[$k];
        $sample=~s/\"//g;

        $tpm_counts{$gene}{$sample} = $tpm;
        my $hist_of_sample = $bs_id_hist{$sample};

        ##store number of events per histology
        $hist_check{$gene}{$hist_of_sample}++;
        $k++;
      }
      push @gene_id, $gene;
  }
}
close(FIL);


my @broad_hist_uniq      = do { my %seen; grep { !$seen{$_}++ } @broad_hist };
my @bs_ids_uniq          = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @gene_id_uniq = do { my %seen; grep { !$seen{$_}++ } @gene_id };

#my $out_file = "results/pan_cancer_splicing.thr10.report_select.remDup.txt";
#open(OUT,">",$out_file) || die("Cannot Open File");
print  "Gene_ID";

  foreach my $hist (sort @broad_hist_uniq)
  {
    my $i = 1;
    foreach my $sample (@{$histology_ids{$hist}})
    {
    #  last if($i>30);
      print "\t";

      ##print histology name (renamed)
      #print $hist."_".$i;

      ##print sample name
      print $sample;
      $i++;

    }
  }
print "\n";


foreach my $event (@gene_id_uniq)
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
  print $event;

  foreach my $hist (sort @broad_hist_uniq)
  {

    my $i = 1;
    foreach my $sample (@{$histology_ids{$hist}})
    {

    #  last if($i>30); ## for testing purposes

      print  "\t";
      #print $sample;
      if($tpm_counts{$event}{$sample}>0)
      {
        #print $sample.":".$inc_levels{$event}{$sample};
        print  $tpm_counts{$event}{$sample};
      }
      else
      {
        # print $sample.":0";
        print  "0";
      }
      $i++;
    }
  }
  print "\n";
}

__DATA__
perl create_matrix_of_TPM_removeDups.pl ../../data/pbta-histologies.RNA-Seq.initial.tsv ../../data/polyA_gene_counts_tab.tsv ../../data/stranded_gene_counts_tab.tsv ../../data/pbta-histologies.RNA-Seq.abridged.tsv  > results/pan_cancer_expr_tpm_report_select.remDup.txt
