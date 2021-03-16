#!/usr/bin/perl

my ($cluster, $hist) = ($ARGV[0], $ARGV[1]);

open(CLU, $cluster) || die("Cannot Open File");
while(<CLU>)
{
  chomp;
  my ($bs_id, $cluster) = split "\t";
  $sample_cluster{$bs_id} = $cluster;
}
close(CLU);

open(HIS, $hist) || die("Cannot Open File");
while(<HIS>)
{
  chomp;
  my (@cols) = split "\t";
  my $sample = $cols[2];
  print $_;
  print "\tCluster", $sample_cluster{$sample},"\n";

  #$sample_cluster{$bs_id} = $cluster;
}
close(HIS);

__DATA__
Kids_First_Biospecimen_ID	Cluster
BS_15ETQ0E4	1
BS_1QNTEZPS	1
BS_31TF0P1N	1
BS_4E1G7E3G	2
