#
# construct_CLK1_pos_neg_tpm_table.pl
#

## input data
# normal expression table, tumor expression table, histology file
#~/Desktop/AS-DMG/analyses/TFs_vs_expr/gene_tpm.brain_normals.txt
my ($gene_exp) = ($ARGV[0]);
my (%CLK1_neg, %CLK1_pos);

##store top 5 CLK1-AS associated Samples
open(POS,"CLK1_samples_pos_v2.txt") || die("Cannot Open File");
while(<POS>)
{
  chomp;
  my $sample = $_;
  $sample =~s/\s+//;
  $CLK1_pos{$sample} = $sample;
}
close(POS);

open(NEG,"CLK1_samples_neg_v2.txt") || die("Cannot Open File");
while(<NEG>)
{
  chomp;
  my $sample = $_;
  $sample =~s/\s+//;
  $CLK1_neg{$sample} = $sample;
}
close(NEG);

## hashes to store tpm values
my (%norm_expr, %tumor_polya, %tumor_stranded);

## store genes
my @genes;
my %gene_exists;

## go through tpm value table and store in hash
open(FIL,$gene_exp) || die("Cannot Open File");
while(<FIL>)
{
  chomp;

  if($_=~/gene/)
  {
    my @cols = split ",";
    shift @cols;
    @header = @cols;
    print "gene,";

    foreach my $sample(@header)
    {

      $sample=~s/nonpolyA_dmg_wt\://;
      $sample=~s/nonpolyA_dmg_h3k28\://;

      $sample=~s/polyA_dmg_wt\://;
      $sample=~s/polyA_dmg_h3k28\://;

      #print "sample: ",$sample,"\n";
      if($CLK1_pos{$sample})
      {
        my $sample_rename = $sample."_pos";
        print $sample_rename;
        print ",";

      }

      if($CLK1_neg{$sample})
      {
        $sample_rename = $sample."_neg";
        print $sample_rename;
        print ",";

      }
      push @samples, $sample;
    }
  }
  else
  {
    my @cols = split ",";
    my $gene = shift @cols;
    #print "gene: ",$gene,"\n";
    #next unless $target_list{$gene};

    push @genes, $gene;
    print "\n";
    print $gene,",";

    my $i = 0;
    foreach my $tpm (@cols)
    {
      my $bs_id = $samples[$i];
      #print $bs_id,"\n";
      $bs_id=~s/polyA_dmg_wt\://;
      $bs_id=~s/nonpolyA_dmg_wt\://;
      $bs_id=~s/polyA_dmg_h3k28\://;
      $bs_id=~s/nonpolyA_dmg_h3k28\://;

      #next unless ($CLK1_neg{$bs_id} || $CLK1_pos{$bs_id});
      if($CLK1_pos{$bs_id})
      {
        $bs_id = $bs_id."_pos";
      #  $gene_tpm{$bs_id}{$gene} = $tpm;
        print $tpm;
        print ",";

      }

      if($CLK1_neg{$bs_id})
      {
        $bs_id = $bs_id."_neg";
        #$gene_tpm{$bs_id}{$gene} = $tpm;
        print $tpm;
        print ",";

      }

      #print "bs_id: ",$bs_id,"\n";
      #$gene_tpm{$bs_id}{$gene} = $tpm;
      $i++;
    }
  }
}
close(FIL);

__DATA__

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
