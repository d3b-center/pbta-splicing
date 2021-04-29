#!/usr/bin/perl

# match inclusion levels from rMATS and MAJIQ to assess correlations
# usage: ./corr_rmats_majiq.pl CLK1_IncLevel_brain_tumors.tab.txt CLK1_IncLevel_brain_tumors.hgg.majiq.tab.txt
#

my ($rmats, $majiq) = ($ARGV[0], $ARGV[1]);
my (%rmats_inc,%majiq_inc,%sample_rmats, %sample_majiq);
my @samples;

open(FIL, $majiq) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  if($_=~/vs_(BS_\w+)\./)
  {
    # print $1;
    my @cols = split;
    my $sample = $1;
    my $inc_levels    = $cols[-1];
    my ($first, $second) = split/\;/,$inc_levels;

    $majiq_inc{$sample} = $second;
    $sample_majiq{$sample} = $sample;
    push @samples, $sample;
    #print "majiq: ",$sample,"*\n";

  }

}
close(FIL);

open(FIL, $rmats) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my @cols = split;
  my $sample = $cols[0];
  my $inc    = $cols[-1];
  #print $sample;

  $rmats_inc{$sample} = $inc;
  $sample_rmats{$sample} = $sample;
  push @samples, $sample;
  #print "rmats: ",$sample,"*\n";

}
close(FIL);

my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
foreach my $sample(@samples_uniq)
{
  #print $sample,"\n";
  if($sample_rmats{$sample} && $sample_majiq{$sample})
  {
    print $sample,"\t",$rmats_inc{$sample},"\t",$majiq_inc{$sample},"\n";
  }
}
