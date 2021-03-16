#!/usr/bin/perl

my ($gene_exp, $psi_vals, $gene_name, $target_list) = ($ARGV[0], $ARGV[1],$ARGV[2],$ARGV[3]);

my @header;
my @samples;
my @genes;
my (%gene_tpm, %gene_psi);
my %target_list;

open(TARG,$target_list) ||  die("Cannot Open File");
while(<TARG>)
{
  chomp;
  my $target = $_;
  $target=~s/\s+//;
  #print "target: ",$target,"\n";
  $target_list{$target} = $target;
}
close(TARG);

open(FIL,$gene_exp) || die("Cannot Open File");
while(<FIL>)
{
  chomp;

  if($_=~/gene/)
  {
    my @cols = split ",";
    shift @cols;
    @header = @cols;
    foreach my $sample(@header)
    {
      $sample=~s/nonpolyA_dmg_wt\://;
      $sample=~s/nonpolyA_dmg_h3k28\://;

      $sample=~s/polyA_dmg_wt\://;
      $sample=~s/polyA_dmg_h3k28\://;

      #print "sample: ",$sample,"\n";

      push @samples, $sample;
    }
  }
  else
  {
    my @cols = split ",";
    my $gene = shift @cols;
    #print "gene: ",$gene,"\n";
    next unless $target_list{$gene};

    push @genes, $gene;

    my $i = 0;
    foreach my $tpm (@cols)
    {
      my $bs_id = $samples[$i];
      $bs_id=~s/polyA_dmg_wt\://;
      $bs_id=~s/nonpolyA_dmg_wt\://;
      $bs_id=~s/polyA_dmg_h3k28\://;
      $bs_id=~s/nonpolyA_dmg_h3k28\://;

      #print "bs_id: ",$bs_id,"\n";
      $gene_tpm{$bs_id}{$gene} = $tpm;
      $i++;
    }
  }
}
close(FIL);

open(FIL,$psi_vals) || die("Cannot Open File");
while(<FIL>)
{
  my @cols = split;
  my $path = $cols[0];
  my $gene = $cols[1];
  my $psi  = $cols[-1];
  my $sample = "";

  if($path=~/vs_(BS_\w+)\./)
  {
    $sample = $1;
  }
  else{
    $sample = $path;
  }
  my @psi_vals = split/\;/,$psi;

  #print $sample,"*\t",$gene,"\t",$psi_vals[0],"\t", $psi,"\n";
  #next unless $target_list{$gene};

  $gene_psi{$sample}{$gene} = $psi_vals[0];

  ##rmats version
  $gene_psi{$sample}{$gene} = $psi;
}

print "PSI\n";
foreach my $sample(@samples)
{
 # print $sample,"\t";
 if($gene_psi{$sample}{"NFIC"})
{
  print $gene_psi{$sample}{"NFIC"},"\n";
}
else{print "0\n"; }
}

die;

foreach my $gene (@genes)
{
  print $gene,"_expr,";
}
#print $gene,"TAF1_PSI,";
print "\n";

foreach my $sample(@samples)
{
  foreach my $gene (@genes)
  {
    if($gene_tpm{$sample}{$gene}){
    print $gene_tpm{$sample}{$gene};
  }
  else{
    print "0";
  }
  print ",";
  }
print "\n";
}








__DATA__
TAF1_as_event.tsv
tpm_ordered_by_polyA_dmg_only.csv
gene,polyA_dmg_wt:BS_F0JB4EAK,polyA_dmg_wt:BS_G3NN392N,polyA_dmg_wt:BS_NB9XXBW6,polyA_dmg_wt:BS_R7NTZR4C,polyA_dmg_wt:BS_W7MFJZ5A,polyA_dmg_h3k28:BS_0ZA67BBC,polyA_dmg_h3k28:BS_1N7MQZGR,polyA_dmg_h3k28:BS_21ET39G7,polyA_dmg_h3k28:BS_2JP7RBMB,polyA_dmg_h3k28:BS_49CJNZ06,polyA_dmg_h3k28:BS_5VPM0F36,polyA_dmg_h3k28:BS_68KX6A42,polyA_dmg_h3k28:BS_6M2053M0,polyA_dmg_h3k28:BS_7WM3MNZ0,polyA_dmg_h3k28:BS_8ZY4GST0,polyA_dmg_h3k28:BS_C83TK159,polyA_dmg_h3k28:BS_EZ3147MX,polyA_dmg_h3k28:BS_FCDAH728,polyA_dmg_h3k28:BS_GWSJ4Z9H,polyA_dmg_h3k28:BS_H97S5SQN,polyA_dmg_h3k28:BS_JB43XBCQ,polyA_dmg_h3k28:BS_JQVAWTTM,polyA_dmg_h3k28:BS_KBEX4RT2,polyA_dmg_h3k28:BS_NEVYM2FP,polyA_dmg_h3k28:BS_NGSG2KB6,polyA_dmg_h3k28:BS_QNNX91SM,polyA_dmg_h3k28:BS_R9B92M75,polyA_dmg_h3k28:BS_X4DD4KSZ,polyA_dmg_h3k28:BS_XGDPK33A,polyA_dmg_h3k28:BS_XM1AHBDJ,polyA_dmg_h3k28:BS_YDEVMD24,polyA_dmg_h3k28:BS_Z3RCA1T9,polyA_dmg_h3k28:BS_ZF6BSFNF,nonpolyA_dmg_wt:BS_30VC3R2Q,nonpolyA_dmg_wt:BS_5965SFPZ,nonpolyA_dmg_wt:BS_6DCSD5Y6,nonpolyA_dmg_wt:BS_B1C6GZ84,nonpolyA_dmg_wt:BS_PZVHMSYN,nonpolyA_dmg_wt:BS_W4DB5RP1,nonpolyA_dmg_h3k28:BS_1D6PZNKN,nonpolyA_dmg_h3k28:BS_52MPKTJH,nonpolyA_dmg_h3k28:BS_7GJT4A6A,nonpolyA_dmg_h3k28:BS_97M1E2DW,nonpolyA_dmg_h3k28:BS_A7Q8G0Y1,nonpolyA_dmg_h3k28:BS_MKM0EEN1,nonpolyA_dmg_h3k28:BS_MN5XVP5A,nonpolyA_dmg_h3k28:BS_PQVTB6BB,nonpolyA_dmg_h3k28:BS_QB84TBA3,nonpolyA_dmg_h3k28:BS_RCCN63K6,nonpolyA_dmg_h3k28:BS_RXN5T5YT,nonpolyA_dmg_h3k28:BS_T1YHK6B7,nonpolyA_dmg_h3k28:BS_TGS109T3,nonpolyA_dmg_h3k28:BS_TXZ7KHTQ,nonpolyA_dmg_h3k28:BS_WH8G4VFB,nonpolyA_dmg_h3k28:BS_YEW22W35
/Users/naqvia/Desktop/neoepitope_HGG/HGG-DIPG-neoepitope-non-denovo/control-BS/control-BS_vs_BS_052PZFMK.94a6faf8-322a-48f0-b153-05f50e975b56.neoepitope_re-analysis.non_denovo.voila_tsv.tsv	TAF1	0.007;0.993
