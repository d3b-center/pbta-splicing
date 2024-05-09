#!/usr/bin/perl

use strict;
use warnings;

# Generate for ggplot for all events corresponding to functional sites
open(my $out_pos, '>', 'results/splicing_events.SE.total.pos.intersectunip.ggplot.txt') or die "Cannot open file: $!";
print $out_pos "SpliceID\tdPSI\tUniprot\n";

open(my $out_neg, '>', 'results/splicing_events.SE.total.neg.intersectunip.ggplot.txt') or die "Cannot open file: $!";
print $out_neg "SpliceID\tdPSI\tUniprot\n";

sub process_file {
    my ($prefix, $suffix, $category, $output) = @_;

    open(my $in, '<', "results/splicing_events.SE.total.$prefix.intersectunip$suffix.hg38.col.wo.txt") or die "Cannot open file: $!";
    my %seen;
    while (my $line = <$in>) {
        chomp $line;
        my ($splice_id, $dpsi) = (split /\s+/, $line)[3, 4];
        next if $seen{$splice_id}++; # Skip duplicates
        print $output "$splice_id\t$dpsi\t$category\n";
    }
    close($in);
}

# Process positive files
process_file('pos', 'Mod', 'Modifications', $out_pos);
process_file('pos', 'Other', 'Other', $out_pos);
process_file('pos', 'DisulfBond', 'DisulfBond', $out_pos);
process_file('pos', 'LocSignal', 'LocSignal', $out_pos);

# Process negative files
process_file('neg', 'Mod', 'Modifications', $out_neg);
process_file('neg', 'Other', 'Other', $out_neg);
process_file('neg', 'DisulfBond', 'DisulfBond', $out_neg);
process_file('neg', 'LocSignal', 'LocSignal', $out_neg);

close($out_pos);
close($out_neg);
