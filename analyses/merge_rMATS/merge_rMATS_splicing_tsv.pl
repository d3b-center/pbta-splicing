#!/usr/bin/perl

#######################################################
# merge_rMATS_splicing_tsv.pl
#
# descr: merge all of rMATS tsv col_type
# usage: ./merge_rMATS_splicing_tsv.pl <config_file>
######################################################
my $config_file = $ARGV[0]; ## later iteration

# my ($SE_dir,$MXE_dir,$A5SS_dir,$A3SS_dir, $A3SS_dir,$RI_dir);

## directories of rMATS run (splicing types need to be in different dirs)
my @SE_dir=</Users/naqvia/Desktop/rmats_run/SE/*SE.MATS.JC.txt>;
my @MXE_dir=</Users/naqvia/Desktop/rmats_run/MXE/*.MATS.JC.txt>;
my @A5SS_dir=</Users/naqvia/Desktop/rmats_run/A5SS/*.MATS.JC.txt>;
my @A3SS_dir=</Users/naqvia/Desktop/rmats_run/A3SS/*.MATS.JC.txt>;
my @RI_dir=</Users/naqvia/Desktop/rmats_run/RI/*MATS.JC.txt>;


#my @rmats_tsv = </Users/naqvia/Desktop/rmats_run/SE/*control-BS*JC.txt>;

## full header for combining all splicing cases
my @full_header =("GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","1stExonStart_0base","1stExonEnd","2ndExonStart_0base","2ndExonEnd","riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE",	"longExonStart_0base","longExonEnd","shortES","shortEE",
                  "flankingES",	"flankingEE", "IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2","IncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference");

##print full header/col names/fields
print "splicing_case\tsample\tctrl\t";
print join "\t", @full_header,"\n";

## go through SE rMATS run files // can prob have each splicing type a diff script where you append to file
my @header;
foreach my $SE (@SE_dir)
{
  my $sample;
  my $ctrl_sample;

  if($SE=~/\_(BS\_\w+)\./)
  {
    $sample = $1;
  }
  if($SE=~/\.(.+)\_vs/)
  {
    $ctrl_sample = $1;
  }

  open(FIL,$SE) || die("Cannot Open File");
  while(<FIL>)
  {
    chomp;

    if($_=~/^ID/)
    {
      @header = split "\t";
    }
    else{
      my @cols = split "\t";

      ## ignore non-signficant results
      my $pval = $cols[18];
      next unless $pval <= 0.05;
      #print join " ", @header,":HEADER\n";
      #print "\n";
      print "SE\t";
      print $sample,"\t",$ctrl_sample,"\t";
      foreach my $full_header(@full_header)
      {
        my $i = 0;
        ## skip columns where it is not relevant for a given splicing case
        if(   ($full_header eq "1stExonStart_0base") ||($full_header eq "1stExonEnd") || ($full_header eq "2ndExonStart_0base") || ($full_header eq "2ndExonEnd")
            || ($full_header eq "riExonStart_0base") || ($full_header eq "riExonEnd") || ($full_header eq "longExonStart_0base")
            || ($full_header eq "longExonEnd") || ($full_header eq "shortES") || ($full_header eq "shortEE") || ($full_header eq "flankingES")
            || ($full_header eq "flankingEE") )
            {
              print "NA\t";
              next;
            }
        else{
          foreach my $item (@cols){
            my $current_header_name = $header[$i];
            #print "header: ".$current_header_name,":",$full_header,"\n";
            if($current_header_name =~/$full_header/)
            {
              print $cols[$i],"\t";
            }
            else
            {
              #print "NA";
            }
            #print "\t";
            $i++;
          }
        }
      }
      print "\n";
    }
  }
}

## go through MXE rMATS run files
my @header;
foreach my $MXE (@MXE_dir)
{
  my $sample;
  my $ctrl_sample;

  if($MXE=~/\_(BS\_\w+)\./)
  {
    $sample = $1;
  }
  if($MXE=~/\.(.+)\_vs/)
  {
    $ctrl_sample = $1;
  }

  open(FIL,$MXE) || die("Cannot Open File");
  while(<FIL>)
  {
    chomp;

    if($_=~/^ID/)
    {
      @header = split "\t";
    }
    else{
      my @cols = split "\t";
      my $pval = $cols[20];
      next unless $pval <= 0.05;
      print "MXE\t";
      print $sample,"\t",$ctrl_sample,"\t";

      foreach my $full_header(@full_header)
      {

        if( ($full_header eq "exonStart_0base") || ($full_header eq "exonEnd")
            || ($full_header eq "riExonStart_0base") || ($full_header eq "riExonEnd") || ($full_header eq "longExonStart_0base")
            || ($full_header eq "longExonEnd") || ($full_header eq "shortES") || ($full_header eq "shortEE") || ($full_header eq "flankingES")
            || ($full_header eq "flankingEE") )
            {
              #print "i:",$i," ";
              #print $full_header,":";

              print "NA\t";

              next;
            }
            my $i = 0; #

        foreach my $item (@cols){
          my $current_header_name = $header[$i];
          #print "header: ".$current_header_name,":",$full_header,"\n";
          if($current_header_name =~/$full_header/)
          {
            print $cols[$i],"\t";
          }
          else
          {
            #print "NA";
          }
          #print "\t";
          $i++;
        }
      }
      print "\n";
    }
  }
}

## go through A5SS rMATS run files
foreach my $A5SS (@A5SS_dir)
{
  #print $SE,"\n";
  my $sample;
  my $ctrl_sample;

  if($A5SS=~/\_(BS\_\w+)\./)
  {
    $sample = $1;
  }
  if($A5SS=~/\.(.+)\_vs/)
  {
    $ctrl_sample = $1;
  }

  open(FIL,$A5SS) || die("Cannot Open File");
  while(<FIL>)
  {
    chomp;

    if($_=~/^ID/)
    {
      @header = split "\t";
    }
    else{
      my @cols = split "\t";
      my $pval = $cols[18];
      next unless $pval <= 0.05;

      print "A5SS\t";
      print $sample,"\t",$ctrl_sample,"\t";
      foreach my $full_header(@full_header)
      {
        my $i = 0; #
        if(($full_header eq "exonStart_0base") || ($full_header eq "exonEnd")
        || ($full_header eq "upstreamES") || ($full_header eq "upstreamEE") || ($full_header eq "downstreamES")	||($full_header eq "downstreamEE")
        || ($full_header eq "1stExonStart_0base") || ($full_header eq "1stExonEnd") || ($full_header eq "2ndExonStart_0base") || ($full_header eq "2ndExonEnd")
        || ($full_header eq "riExonStart_0base") || ($full_header eq "riExonEnd") )
            {
              print "NA\t";
              next;
            }
        foreach my $item (@cols){
          my $current_header_name = $header[$i];
          #print "header: ".$current_header_name,":",$full_header,"\n";
          if($current_header_name =~/$full_header/)
          {
            print $cols[$i],"\t";
          }
          else
          {
            #print "NA";
          }
          #print "\t";
          $i++;
        }
      }
      print "\n";
    }
  }
}

## go through A3SS rMATS run files
foreach my $A3SS (@A3SS_dir)
{
  #print $SE,"\n";
  my $sample;
  my $ctrl_sample;

  if($A3SS=~/\_(BS\_\w+)\./)
  {
    $sample = $1;
  }
  if($A3SS=~/\.(.+)\_vs/)
  {
    $ctrl_sample = $1;
  }

  open(FIL,$A3SS) || die("Cannot Open File");
  while(<FIL>)
  {
    chomp;

    if($_=~/^ID/)
    {
      @header = split "\t";
    }
    else{
      my @cols = split "\t";
      my $pval = $cols[18];
      next unless $pval <= 0.05;

      print "A3SS\t";
      print $sample,"\t",$ctrl_sample,"\t";
      foreach my $full_header(@full_header)
      {
        my $i = 0; #
        if(($full_header eq "exonStart_0base") || ($full_header eq "exonEnd")
        || ($full_header eq "upstreamES") || ($full_header eq "upstreamEE") || ($full_header eq "downstreamES")	||($full_header eq "downstreamEE")
        || ($full_header eq "1stExonStart_0base") || ($full_header eq "1stExonEnd") || ($full_header eq "2ndExonStart_0base") || ($full_header eq "2ndExonEnd")
        || ($full_header eq "riExonStart_0base") || ($full_header eq "riExonEnd") )
            {
              print "NA\t";
              next;
            }
        foreach my $item (@cols){
          my $current_header_name = $header[$i];
          #print "header: ".$current_header_name,":",$full_header,"\n";
          if($current_header_name =~/$full_header/)
          {
            print $cols[$i],"\t";
          }
          else
          {
            #print "NA";
          }
          #print "\t";
          $i++;
        }
      }
      print "\n";
    }
  }
}

## go through RI rMATS run files
foreach my $RI (@RI_dir)
{
  #print $SE,"\n";
  my $sample;
  my $ctrl_sample;

  if($RI=~/\_(BS\_\w+)\./)
  {
    $sample = $1;
  }
  if($RI=~/\.(.+)\_vs/)
  {
    $ctrl_sample = $1;
  }

  open(FIL,$RI) || die("Cannot Open File");
  while(<FIL>)
  {
    chomp;

    if($_=~/^ID/)
    {
      @header = split "\t";
    }
    else{
      my @cols = split "\t";
      my $pval = $cols[18];
      next unless $pval <= 0.05;

      print "RI\t";
      print $sample,"\t",$ctrl_sample,"\t";
      foreach my $full_header(@full_header)
      {
        my $i = 0; #
        if(($full_header eq "exonStart_0base") || ($full_header eq "exonEnd")
        || ($full_header eq "1stExonStart_0base") || ($full_header eq "1stExonEnd") || ($full_header eq "2ndExonStart_0base") || ($full_header eq "2ndExonEnd")
        || ($full_header eq "longExonStart_0base")||($full_header eq "longExonEnd") || ($full_header eq "shortES") || ($full_header eq "shortEE") || ($full_header eq "flankingES") || ($full_header eq "flankingEE") )
            {
              print "NA\t";
              next;
            }
        foreach my $item (@cols){
          my $current_header_name = $header[$i];
          #print "header: ".$current_header_name,":",$full_header,"\n";
          if($current_header_name =~/$full_header/)
          {
            print $cols[$i],"\t";
          }
          else
          {
            #print "NA";
          }
          #print "\t";
          $i++;
        }
      }
      print "\n";
    }
  }
}

__DATA__
config file:
histology file = <file>
SE =<dir>
MXE =<dir>
A5SS =<dir>
A3SS =<dir>
RI =<dir>

new_cols --> <splice_case>
         --> <sample_name>
         --> <ctrl name>
