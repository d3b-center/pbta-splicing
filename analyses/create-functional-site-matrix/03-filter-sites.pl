#!/usr/bin/perl

my ($pri_file, $func_file) = @ARGV[0, 1];
my %filter_event;

open(FIL, $func_file) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  next if($_=~/SpliceID/);
  my @cols = split "\t";
  my $splice_id = $cols[0];
  $filter_event{$splice_id} = $cols[1];
}
close(FIL);

open(IN,$pri_file) || die("Cannot Open File");
while(<IN>)
{
  chomp;
  if($_=~/Splice_ID/)
  {
    print $_,"\n";
  }
  my @cols = split "\t";
  my $splice_id = $cols[0];

  next unless $filter_event{$splice_id};
  print $_,"\n";
}

close(IN);
