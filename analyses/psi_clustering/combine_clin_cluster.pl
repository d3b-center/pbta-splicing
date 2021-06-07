#!/usr/bin/perl
################################################################################
# combine_clin_cluster.pl.pl
# create tsv file of TPM for each sample, to be used for consensus clustering
# written by Ammar Naqvi
#
# usage: ./create_matrix_of_TPM_removeDups.pl <clsuter file> <histology file>
################################################################################
my ($cluster, $hist) = ($ARGV[0], $ARGV[1]);
my %sample_cluster;

open(CLU, $cluster) || die("Cannot Open File");
while(<CLU>)
{
  chomp;
  my ($bs_id, $cluster) = split "\t";
  $bs_id=~s/\"//g;

  $sample_cluster{$bs_id} = $cluster;
}
close(CLU);

open(OUT,">input/pbta-histologies_w_clusters.v2.tsv");

print  "Kids_First_Biospecimen_ID\tCNS_region\tsample_id\taliquot_id\tKids_First_Participant_ID\texperimental_strategy\tsample_type\tcomposition\ttumor_descriptor\tprimary_site\treported_gender\trace\tethnicity\tage_at_diagnosis_days\tpathology_diagnosis\tRNA_library\tOS_days\tOS_status\tcohort\tage_last_update_days\tseq_center\tparent_aliquot_id\tcancer_predispositions\tpathology_free_text_diagnosis\tcohort_participant_id\textent_of_tumor_resection\tgermline_sex_estimate\tnormal_fraction\ttumor_fraction\ttumor_ploidy\tmolecular_subtype\tintegrated_diagnosis\tNotes\tharmonized_diagnosis\tbroad_histology\tshort_histology\tCluster\tColor\n";

open(HIS, $hist) || die("Cannot Open File");
while(<HIS>)
{
  chomp;
  my (@cols) = split "\t";
  my $sample = $cols[0];
  next unless $sample_cluster{$sample};
  print  $_;
  print  "\tCluster", $sample_cluster{$sample},"\t";
  if($sample_cluster{$sample}=~/1/) {  print "FF0000\n";  }
  if($sample_cluster{$sample}=~/2/) {  print "00CC00\n";  }
  if($sample_cluster{$sample}=~/3/) {  print "0000FF\n";  }



  #$sample_cluster{$bs_id} = $cluster;
}
close(HIS);
close(OUT);

__DATA__
Kids_First_Biospecimen_ID	Cluster
BS_15ETQ0E4	1
BS_1QNTEZPS	1
BS_31TF0P1N	1
BS_4E1G7E3G	2
