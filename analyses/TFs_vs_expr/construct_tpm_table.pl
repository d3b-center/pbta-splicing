#
# construct_tpm_table.pl
#

## input data
# normal expression table, tumor expression table, histology file
#~/Desktop/AS-DMG/analyses/TFs_vs_expr/gene_tpm.brain_normals.txt
my ($norm_expr, $tumor_expr, $hist) = ($ARGV[0], $ARGV[1], $ARGV[2]);

## hashes to store tpm values
my (%norm_expr, %tumor_polya, %tumor_stranded);

## store genes
my @genes;
my %gene_exists;

## go through normal tpm value table and store in hash
open(FIL,$norm_expr) || die("Cannot Open File $gene_expr");
while(<FIL>)
{
  chomp;
  if($_=~/EN/)
  {
    my @cols = split "\t";
    my $gene = shift @cols;
    $gene=~s/ENSG\d+\_//;

    next if $gene_exists{$gene};
    $gene_exists{$gene} =1;

    push @genes, $gene;
    foreach my $tpm (@cols)
    {
      my $rounded = sprintf "%.0f", $tpm;

      #print $gene,"\t",$rounded,"\n";
      push @{$norm_expr{$gene}}, $rounded;
    }
  }
}
# ENSG00000000003_TSPAN6	5.10	17.90	6.23	6.80	5.51	0.00	0.00	0.00	8.10

## go through tumor tpm value table and store in hash
## divide into stranded vs polyA // tpm_ordered_by_polyA_dmg_only.csv
open(FIL,$tumor_expr) || die("Cannot Open File $gene_expr");
while(<FIL>)
{
  chomp;
  if($_=~/gene/)
  {
    my @cols = split ",";
    print join(",", @cols);
    print ",healthy_1,healthy_2,healthy_3,healthy_4,healthy_5,healthy_6,healthy_7,healthy_8,healthy_9";
  }
  else
  {
    my @cols = split ",";
    print join(",", @cols);
    print ",";
    my $gene = $cols[0];
    if(@{$norm_expr{$gene}})
    {
      print join(",", @{$norm_expr{$gene}});
    }
    else{
      print "0,0,0,0,0,0,0,0,0";
    }
  }
  print "\n";
}

## make gene list unique
## go through histology file
## store histology info, polyA status, and subtype

## report tpm values for each gene and make table with normals and DMG (def) or all HGGs
