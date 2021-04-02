#!/usr/bin/perl

my ($samples_hist, $clk1_tab) = ($ARGV[0], $ARGV[1]);
my %sample_hist;

open(HIST,$samples_hist) || die("Cannot Open File");
while(<HIST>)
{
  chomp;
  next if($_=~/sample/);
  my ($sample, $hist) = split "\t";
  $sample_hist{$sample} = $hist;
}
close(HIST);

print "Sample\tInclusion\tHist\n";

open(CLK1, $clk1_tab) || die("Cannot Open File");
while(<CLK1>)
{
  chomp;
  next if($_=~/Sample/);
  my @cols = split "\t";
  my $sample = $cols[0];
  my $inc_level  = $cols[-1];
  next unless $sample_hist{$sample};
  print $sample,"\t",$inc_level,"\t", $sample_hist{$sample},"\n";
}


__DATA__
grep "\"CLK1\"" ~/Desktop/rmats_run/SE/*control-BS*JC.txt | awk '{if($19 <= 0.05) { print $0 }}' | awk -F "\t" '{ if( ($NF<=-.10) || ($NF>=.10))  { print $0}}' | grep "200860124\t200860215\t200859679\t200859746\t200861237\t200861466"| awk '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' | perl -pe 'if($_=~/\_(BS\_\w+)\.SE/){ print $1,"\t",$_}' | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | perl -pe 's/\"//g' | grep -v "User"
grep "\"CLK1\"" ~/Desktop/pan_cancer_rmats/*freads10.txt |  awk -F "\t" '{ print $0} ' | grep "200860124\t200860215\t200859679\t200859746\t200861237\t200861466"| awk -F "\t" '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$21}' | perl -pe 'if($_=~/\.(BS\_\w+)\.non/){ print $1,"\t",$_}' | awk -F "\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$NF}' | perl -pe 's/\"//g' | grep -v "User" > CLK1_IncLevel_brain_tumors.tab.txt
