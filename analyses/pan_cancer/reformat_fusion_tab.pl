#
# reformat_fusion_tab.pl
#

my ($matrix_fil, $id_mappings_rnaseq, $id_mappings_wgs) = ($ARGV[0], $ARGV[1], $ARGV[2]);


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

open(FIL, $matrix_fil) || die("Cannot Open File $matrix_fil");
while(<FIL>)
{
  chomp;
  my ($sample,$fusion_name) = split;
  my ($geneA, $geneB) = split/\-\-/,$fusion_name;
  next unless ($rnaseq_to_wgs_sample{$sample});
  print $rnaseq_to_wgs_sample{$sample},"\t",$geneA,"\tFusion\tOTHER\n";
  print $rnaseq_to_wgs_sample{$sample},"\t",$geneB,"\tFusion\tOTHER\n";

}
close(FIL);


__DATA__
Tumor_Sample_Barcode   Hugo_Symbol Variant_Classification Variant_Type
