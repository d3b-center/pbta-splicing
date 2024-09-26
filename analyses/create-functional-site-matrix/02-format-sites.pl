#!/usr/bin/perl

use strict;
use warnings;

# Check if the correct number of arguments is provided
my $num_args = 5;
if (@ARGV < $num_args) {
    die "Usage: $0 signal_file disulf_file modif_file other_file output_file\n";
}

# Assign command-line arguments to variables
my ($file_signal, $file_disulf, $file_mod, $file_other, $output_file) = @ARGV[0..$num_args-1];

# Process files
open(OUT, '>', $output_file) or die "Cannot open file: $!";
print OUT "SpliceID\t\tUniprot\n";
open(IN, $file_signal) or die "Cannot open file: $!";
my %seen;
while (my $line = <IN>)
{
    chomp $line;
    my ($splice_id) = (split /\s+/, $line)[3];
    next if $seen{$splice_id}++; # Skip duplicates
    print OUT $splice_id."\tLocSignal\n";
}
close(IN);

open(IN, $file_disulf) or die "Cannot open file: $!";
%seen=();
while (my $line = <IN>)
{
    chomp $line;
    my ($splice_id) = (split /\s+/, $line)[3];
    next if $seen{$splice_id}++; # Skip duplicates
    print OUT $splice_id."\tDisulfBond\n";
}
close(IN);

open(IN, $file_mod) or die "Cannot open file: $!";
%seen=();
while (my $line = <IN>)
{
    chomp $line;
    my ($splice_id) = (split /\s+/, $line)[3];
    next if $seen{$splice_id}++; # Skip duplicates
    print OUT $splice_id."\tMod\n";
}
close(IN);

open(IN, $file_other) or die "Cannot open file: $!";
%seen=();
while (my $line = <IN>)
{
    chomp $line;
    my ($splice_id) = (split /\s+/, $line)[3];
    next if $seen{$splice_id}++; # Skip duplicates
    print OUT $splice_id."\tOther\n";
}
close(IN);


close(OUT);
