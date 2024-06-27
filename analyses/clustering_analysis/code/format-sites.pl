#!/usr/bin/perl

use strict;
use warnings;


# Process files
open(OUT, '>', 'input/input/pan-cancer-SE-func.tsv.tmp') or die "Cannot open file: $!";
print OUT "SpliceID\t\tUniprot\n";
open(IN, "input/pan-cancer-SE-func.tsv.tmp.intersectunipLocSignal.hg38.col.wo.txt") or die "Cannot open file: $!";
my %seen;
while (my $line = <IN>)
{
    chomp $line;
    my ($splice_id) = (split /\s+/, $line)[3];
    next if $seen{$splice_id}++; # Skip duplicates
    print OUT $splice_id."\tLocSignal\n";
}
close(IN);

open(IN, "input/pan-cancer-SE-func.tsv.tmp.intersectunipDisulfBond.hg38.col.wo.txt") or die "Cannot open file: $!";
my %seen;
while (my $line = <IN>)
{
    chomp $line;
    my ($splice_id) = (split /\s+/, $line)[3];
    next if $seen{$splice_id}++; # Skip duplicates
    print OUT $splice_id."\tDisulfBond\n";
}
close(IN);

open(IN, "input/pan-cancer-SE-func.tsv.tmp.intersectunipMod.hg38.col.wo.txt") or die "Cannot open file: $!";
my %seen;
while (my $line = <IN>)
{
    chomp $line;
    my ($splice_id) = (split /\s+/, $line)[3];
    next if $seen{$splice_id}++; # Skip duplicates
    print OUT $splice_id."\tMod\n";
}
close(IN);

open(IN, "input/pan-cancer-SE-func.tsv.tmp.intersectunipOther.hg38.col.wo.txt") or die "Cannot open file: $!";
my %seen;
while (my $line = <IN>)
{
    chomp $line;
    my ($splice_id) = (split /\s+/, $line)[3];
    next if $seen{$splice_id}++; # Skip duplicates
    print OUT $splice_id."\tOther\n";
}
close(IN);


close(OUT);
