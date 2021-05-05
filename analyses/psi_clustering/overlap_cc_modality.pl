#!/usr/bin/perl

#
# overlap_cc_modality.pl
#

my ($splicing, $gene_expr) = ($ARGV[0], $ARGV[1]);
my %splice_cc_ids;

open(SPL,$splicing) || die("Cannot Open File");
while(<SPL>)
{
  chomp;
  my ($bs_id, $cluster) = split;
  $cc_membership_psi{$bs_id} = $cluster;
  $splice_cc_ids{$bs_id} = $bs_id;
}
close(SPL);

open(EXP,$gene_expr) || die("Cannot Open File");
while(<EXP>)
{
  chomp;
  my (@cols) = split;
  my $cluster  = $cols[1];
  my $bs_id    = $cols[4];

  #print $bs_id,"\n";
  next unless $splice_cc_ids{$bs_id};
  $cc_membership_expr{$bs_id} = $cluster;

  print $bs_id,"\t",$cluster,"\t",$cc_membership_psi{$bs_id},"\n";
}
close(EXP);

__DATA__
results/CC_groups.remDups.txt
BS_850BAHH9	1
