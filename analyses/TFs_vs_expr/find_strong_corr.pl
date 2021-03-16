#!/usr/bin/perl

#
#
# usage: perl scripts/find_strong_corr.pl pval_SFs_v2.txt rcorr_SFs_v2.txt
#

my ($file_pval, $file_r2) = ($ARGV[0], $ARGV[1]);
my @genes_expr;
my (%rcorrs, %pvals);

open(FIL,$file_pval) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  if ($_=~/^\s/)
  {
    chomp;
    my @cols = split "\t";
    shift @cols;
    foreach $item(@cols)
    {
      #next unless ($item=~/expr/);
      #print $item,"\t";
      push @genes_expr,$item;
    }
  }
  else
  {
    my @cols = split "\t";
    my $name = shift @cols;
    #next if ($name=~/expr/);
    #print "name:".$name,"\t";
    my $i=0;
    foreach my $pval(@cols)
    {
    #  print $genes_expr[$i],":";
      my $against_item = $genes_expr[$i];
    #  print $pval,"\t";
      push @pvals,$pval;
      $pvals{$name}{$against_item} = $pval;
      #$rcorrs{$name}{$against_item} = $pval;
      $i++;
    }
    print "\n";
  }
}

open(FIL,$file_r2) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  if ($_=~/^\s/)
  {
    chomp;
    my @cols = split "\t";
    shift @cols;

    foreach $item(@cols)
    {
      #next unless ($item=~/expr/);
      #print $item,"\t";
      push @genes_expr,$item;
    }
  }
  else
  {
    my @cols = split "\t";
    my $name = shift @cols;
    #next if ($name=~/expr/);
    #print "name:".$name,"\t";
    my $i=0;
    foreach my $pval(@cols)
    {
      #print $genes_expr[$i],":";
      my $against_item = $genes_expr[$i];
      #print $pval,"\t";
      push @pvals,$pval;
      #$pvals{$name}{$against_item} = $pval;
      $rcorrs{$name}{$against_item} = $pval;
      $i++;
    }
    #print "\n";
  }
}

foreach my $item(keys %pvals)
{
  #next unless $item=~/expr/;
  foreach my $second_item(keys %{$pvals{$item}})
  {
    #next unless($pvals{$item}{$second_item} <= 0.05);
    #next if $second_item=~/expr/;

    print $item,"\t",$second_item,"\t",$rcorrs{$item}{$second_item},"\t",$pvals{$item}{$second_item},"\n";
  }
}
