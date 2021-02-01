#!/usr/bin/perl

#
# reformat_and_filter_psi_table_for_oncoprint.pl
#

my ($matrix_fil, $id_mappings_rnaseq, $id_mappings_wgs) = ($ARGV[0], $ARGV[1], $ARGV[2]);
my %sample_link_rnaseq;
my %rnaseq_to_wgs_sample;

## link RNA to WGS Samples
open(RNA,$id_mappings_rnaseq) || die("Cannot Open File $id_mappings_rnaseq");
while(<RNA>)
{
  chomp;
  if($_=~/BS/)
  {
    my ($biospecimen, $sample) = split;
    #print $_,"\n";
    $sample_link_rnaseq{$sample} = $biospecimen;
  }
}
close(RNA);

open(DNA,$id_mappings_wgs) || die("Cannot Open File $id_mappings_wgs");
while(<DNA>)
{
  chomp;
  if($_=~/BS/)
  {
    my ($biospecimen, $sample) = split;
    my $rnaseq_sample = $sample_link_rnaseq{$sample};
    if ($rnaseq_sample)
    {
      $rnaseq_to_wgs_sample{$rnaseq_sample} = $biospecimen;
    }
  }
}
close(DNA);


my %sample_gene;

my @samples = "";
open(FIL,$matrix_fil) || die("Cannot Open File");
while(<FIL>)
{
  chomp;

  if($_=~/(BS_\w+)/)
  {
    @samples = split;
  }
  else{
    @psi = split;
    my $i = 0;
    my $gene = "";
    foreach my $psi(@psi)
    {
      if ($i == 0)
      {
        my (@splice_id) = split/\_/,$psi;
        $gene = $splice_id[0];
        $i++;
        next;
      }
      else{
        if($psi >.10)
        {
          my $sample = $samples[$i];

          ##check to see if there corresponding WGS, skip if not
          unless ($rnaseq_to_wgs_sample{$sample})
          {
            $i++;
            next;
          }

          ## already found a splicing event in gene?
          next if($sample_gene{$sample}{$gene});

          my $wgs_id = $rnaseq_to_wgs_sample{$sample};

          #print $samples[$i],"\t";

          #pring wgs id to match maf
          print $wgs_id,"\t";
          print $gene, "\t";
          print "Splicing\tOTHER\n";
          $sample_gene{$sample}{$gene} = 1;
        }
      }
      $i++;
    }
  }
}
