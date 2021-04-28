#!/usr/bin/perl



#
# extract_recurrent_dominant_lsvs.pl
#
# extract and append splicing ids from rMATS to be used to overlap with Uniprot
#

#my @tsv_files = <majiq_voila_tsv/*.tsv>;

## for denovo
#my @tsv_files = <majiq_denovo_voila_tsv/*.tsv>;

## for B-ALL
#my @tsv_files = <voila_deltapsi_disableIntron_disableDenovo_gcv31/ProB_*.tsv>;

#my @tsv_files = <voila_tsv/*non_denovo.lgg_wg_anno.voila_tsv.tsv>;

my $histology = $ARGV[0];
my @rmats_tsv = </Users/naqvia/Desktop/rmats_run/SE/*control-BS*JC.txt>;

my %dmg_samples;
my (@splicing_events, @splicing_event_deltapsis, @samples);
my %splicing_event_counts;
my %inc_levels;


## hash data structure to store values
my %gene_lsv;
my %lsv_ctrl_psi;
my %lsv_case_psi;
my (@samples, @splicing_ids);
my %dpsi;


my (%chr,%num_exons,%str);

my @tsv_files;
my %ids_dmg;

print "processing histology...\n";

## open histology file and annotate DMG samples
##add only dmg $samples
open(FIL,$histology) || die("Cannot Open File");
while(<FIL>)
{
  chomp;

  ## skip unless RNA-seq
  next unless ($_=~/RNA\-Seq/);

  my @cols = split "\t";
  my $molecular_subtype = $cols[30];
  next unless ($molecular_subtype=~/DMG/);
  my $bs_id = $cols[0];
  $dmg_samples{$bs_id} = $bs_id;
}
close(FIL);

print "processing rMAT files...\n";
foreach my $file(@rmats_tsv)
{
  my $bs_id = "";
  if($file=~/vs\_(BS\_\w+)\./)
  {
    #print "samp*".$1,"*\n";
    $bs_id = $1;
  }
  next unless $dmg_samples{$bs_id};
  push @samples, $bs_id;

  open(FIL,$file) || die("Cannot Open File");
  while(<FIL>)
  {
    chomp;
    my @cols = split "\t";
    my $gene         = $cols[2];

    ##exon coordinates
    my $chr          = $cols[3];
    my $str          = $cols[4];
    my $exonStart    = $cols[5]+1;
    my $exonEnd      = $cols[6];
    my $upstreamES   = $cols[7];
    my $upstreamEE   = $cols[8];
    my $downstreamES = $cols[9];
    my $downstreamEE = $cols[10];
    my $inc_level    = $cols[20];

    my $ctrl_IJC      = $cols[12];
    my $ctrl_SJC      = $cols[13];
    my $tumor_IJC     = $cols[14];
    my $tumor_SJC     = $cols[15];

    my $inc_len  = $cols[16];
    my $skip_len = $cols[17];
    my $pval = $cols[18];
    my $thr_diff = $cols[-1];

    ## only look at strong changes
    next unless (abs($thr_diff) >= .10);
    next unless ($pval<=0.05);
    next unless ( ($tumor_IJC >=10) || ($tumor_SJC >=10) );

    $gene=~s/\"//g;

    #print $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE,"\t",$inc_level,"*\n";
    my $splice_id= $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;
    my $splice_id_short= $gene."_".$exonStart."-".$exonEnd;

    ##keep one splicing event per exon (choose larger dPSI)
    my $related_event = 0;
    if($inc_levels{$splice_id_short}{$bs_id} > 0)
    {
      $related_event = $inc_levels{$splice_id_short}{$bs_id};
    }
    #print "splice_id_short:".$splice_id_short,"\tdPSI:",$thr_diff,"E\tRE:",$related_event,"\n";
    next unless(abs($thr_diff) > abs($related_event));
    #print "splice_id_short:".$splice_id_short,"\tdPSI:",$thr_diff,"E\tRE:",$related_event,"\n";

    $chr{$splice_id_short} = $chr;
    $str{$splice_id_short} = $str;
    $inc_levels{$splice_id_short}{$bs_id} = $thr_diff;
    push @splicing_events, $splice_id_short;
    push @{$splicing_event_deltapsis{$splice_id_short}}, $thr_diff;
    $splicing_event_counts{$splice_id_short}++;

  }
}
close(FIL);

## make arrays unique
my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

open(TAB, ">splicing_events.dpsi10.jc10.rec2.tsv");
open(BED, ">splicing_events.dpsi10.jc10.rec2.bed");

print "writing output...\n";
foreach my $event (@splicing_events_uniq)
{
  print "event: ",$event,"\t";
  print join "\t",@{$splicing_event_deltapsis{$event}},"\n";

  ## compute avg and stdevs of each lsv dPSI
  my $avg_dpsi = &average(\@{$splicing_event_deltapsis{$event}});
	my $std_dpsi = &stdev(\@{$splicing_event_deltapsis{$event}});

  my ($gene,$exon_coord) = split/\_/,$event;

  next unless $splicing_event_counts{$event} > 2;
  print TAB $gene,"\t",$event,"\t",$avg_dpsi,"\t",$std_dpsi,"\t";
  print TAB $splicing_event_counts{$event};

  print TAB "\t",$chr{$event},":",$exon_coord,"\t",$str{$event},"\n";

  print BED $chr{$event},"\t";
  if($event=~/(\w+)\_(\d+\-\d+)/)
  {
    #$lsv_coord = $1;
    $exon_coord = $2;
    #print BED $exon_coord." ";
    my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
    print BED $start_exon_coord,"\t",$end_exon_coord,"\t",$event,"\t",$avg_dpsi,"\t",$str{$event},"\n";
  }
}
close(TAB);
close(BED);

## chr11	132335046	132335286	NTM_ENSG00000182667.14:s:132314552-132315070_132335046-132335286	0.640240740740741	+
## calculate stdev given an array of PSI values
sub stdev{
        my($data) = @_;
        if(@$data == 1){

                return "1";
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

## calculate avg given an array of PSI values
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}


__DATA__

my @ctrl_bs_files = </Users/naqvia/Desktop/neoepitope_HGG/HGG-DIPG-neoepitope-non-denovo/control-BS/control-BS_vs_BS_*tsv>;
foreach my $file(@ctrl_bs_files)
{
  my @sep = split/\//,$file;
  my $biospecimen = $sep[-1];


  if($file=~/(BS\_\w+)/)
  {
    $biospecimen = $1;
    $biospecimen=~s/BS_vs_//;
  }
  else
  {
    #print "biosHGG:*",$tsv,"*\t*",$hgg_mappings{$tsv} ,"*\n";
    $biospecimen = $hgg_mappings{$file};
  }

  #print "biospec:".$biospecimen,"\n";
  next unless ($filter_wt{$biospecimen} || $filter_mt{$biospecimen});

  push @samples, $biospecimen;
  my $sample = $biospecimen;

  open(FIL,$file) || die("Cannot Open File $file");
  while(<FIL>)
  {
      next if ($_=~/^#/);
      my (@content)  = split "\t";
      my $gene_name  = $content[0];
      my $lsv_coord  = $content[2];
      my $psi_ctrl   = $content[6];
      my $psi_case   = $content[7];
      my $num_exons  = $content[13];
      my $chr        = $content[15];
      my $strand     = $content[16];
      my $exon_coords= $content[18];

      my ($ctrl_psi, $case_psi);

      my @psi_ctrl = split/;/,$psi_ctrl;
      my @psi_case = split/;/,$psi_case;
      my @exon_coord = split/;/,$exon_coords;


      ## remove complex LSVs
      next unless ( scalar(@psi_ctrl) == (scalar(@exon_coord)-1));

      print $gene_name,"\t",$lsv_coord,"\t",$psi_ctrl,"\t",$psi_case,"\t";
      print $exon_coords,"\t";
      #die;

      my $pos = 0;
      my $delta_psi = 0;
      my @relevant_psi_indices;
      my @relevant_dpsi;
      my %dpsi_index;

      ##look for PSI events that are 10% different
      while($pos < scalar(@psi_ctrl))
      {
        #print "enter\n";
        if( abs($psi_ctrl[$pos] - $psi_case[$pos]) >=.10 )
        {
          push @relevant_psi_indices, $pos;
          my $delta_psi = abs($psi_ctrl[$pos] - $psi_case[$pos]);
          push @relevant_dpsi, $delta_psi;
          $dpsi_index{$delta_psi} = $pos;
        }
        $pos++;
      }

      my @relevant_dpsi_sorted = sort { $a <=> $b } @relevant_dpsi;
      my $max_dpsi = $relevant_dpsi_sorted[-1];
      my $index_of_max_dpsi = $dpsi_index{$max_dpsi};

      print "max dpsi: ",$max_dpsi,"\tmax dpsi index: ",$index_of_max_dpsi,"\t";
      ## get relevant exon coordinates
#     foreach my $psi_index(@relevant_psi_indices)
#      {
        if( ($lsv_coord =~/\:s\:/) && ($strand=~/\+/) )
        {
          my $exon_coord_index = $index_of_max_dpsi + 1 ;
          print "enter max exon: ",$exon_coord[$exon_coord_index],"\n";
          my $new_lsv_id = $gene_name."_".$lsv_coord."_".$exon_coord[$exon_coord_index];
          $num_exons{$new_lsv_id} = $num_exons;
          $chr{$new_lsv_id} = $chr;
          $str{$new_lsv_id} = $strand;

          push @{$relevant_splicing_events{$new_lsv_id}{$sample}},$max_dpsi,
                 $exon_coord[$exon_coord_index];
          $splicing_event_count{$new_lsv_id}++;
          push @splicing_events, $new_lsv_id;
          push @{$splicing_event_deltapsis{$new_lsv_id}}, $max_dpsi;
        }

        # exons are in reverse order compared to psi vals
        elsif( ($lsv_coord =~/\:s\:/) && ($strand=~/\-/) ){
          my $exon_coord_index = $index_of_max_dpsi-1;
          my $num_exons = scalar(@exon_coord);
          my $num_psi   = scalar(@psi_ctrl);
          my $exon_coord_index = $num_exons - $num_psi;

          print "max exon: ",$exon_coord[$exon_coord_index],"\n";
          my $new_lsv_id = $gene_name."_".$lsv_coord."_".$exon_coord[$exon_coord_index];
          $num_exons{$new_lsv_id} = $num_exons;
          $chr{$new_lsv_id} = $chr;
          $str{$new_lsv_id} = $strand;

          push @{$relevant_splicing_events{$new_lsv_id}{$sample}},$max_dpsi,
                 $exon_coord[$exon_coord_index];
          $splicing_event_count{$new_lsv_id}++;
          push @splicing_events, $new_lsv_id;
          push @{$splicing_event_deltapsis{$new_lsv_id}}, $max_dpsi;

        }

        # exons are in reverse order compared to psi vals
        elsif( ($lsv_coord =~/\:t\:/) && ($strand=~/\-/) ){
          print scalar (@exon_coord)."&".scalar (@psi_ctrl),":";
          my $num_exons = scalar(@exon_coord);
          my $num_psi   = scalar(@psi_ctrl);
          my $exon_coord_index = $num_exons - $num_psi;
          print $exon_coord_index." ";
          #my $exon_coord_index = $index_of_max_dpsi;
          print "max exon: ",$exon_coord[$exon_coord_index],"\n";
          my $new_lsv_id = $gene_name."_".$lsv_coord."_".$exon_coord[$exon_coord_index];
          $num_exons{$new_lsv_id} = $num_exons;
          $chr{$new_lsv_id} = $chr;
          $str{$new_lsv_id} = $strand;

          push @{$relevant_splicing_events{$new_lsv_id}{$sample}},$max_dpsi,
                 $exon_coord[$exon_coord_index];
          $splicing_event_count{$new_lsv_id}++;
          push @splicing_events, $new_lsv_id;
          push @{$splicing_event_deltapsis{$new_lsv_id}}, $max_dpsi;

        }
        else{
          my $exon_coord_index = $index_of_max_dpsi;
          print "max exon: ",$exon_coord[$exon_coord_index],"\n";
          my $new_lsv_id = $gene_name."_".$lsv_coord."_".$exon_coord[$exon_coord_index];
          $num_exons{$new_lsv_id} = $num_exons;
          $chr{$new_lsv_id} = $chr;
          $str{$new_lsv_id} = $strand;

          push @{$relevant_splicing_events{$new_lsv_id}{$sample}},$max_dpsi,
                 $exon_coord[$exon_coord_index];
          $splicing_event_count{$new_lsv_id}++;
          push @splicing_events, $new_lsv_id;
          push @{$splicing_event_deltapsis{$new_lsv_id}}, $max_dpsi;
        }
#      }
      #die;
      print "largest_dpsi: ".$delta_psi,"\t";
  }
}

## make arrays unique
my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

foreach my $event (@splicing_events_uniq)
{
  ##compute avg and stdevs of each lsv dPSI
  print "event: ",$event,"\t";
  print join "\t",@{$splicing_event_deltapsis{$event}},"\n";
  my $avg_dpsi = &average(\@{$splicing_event_deltapsis{$event}});
	my $std_dpsi = &stdev(\@{$splicing_event_deltapsis{$event}});

  #next unless $splicing_event_count{$event} < 9;
  print "*",$event,"\t",$avg_dpsi,"\t",$std_dpsi,"\t",$num_exons{$event},"\t";
  print $splicing_event_count{$event};
  my ($gene,$lsv,$exon_coord) = split/\_/,$event;
  print "\t",$chr{$event},":",$exon_coord,"\t",$str{$event},"\n";

  print "bed*".$chr{$event},"\t";
  if($event=~/(\d+\-\d+)\_(\d+\-\d+)/)
  {
    $lsv_coord = $1;
    $exon_coord = $2;
    #print "*".$exon_coord." ";
    my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
    print $start_exon_coord,"\t",$end_exon_coord,"\t",$event,"\t",$avg_dpsi,"\t",$str{$event},"\n";
  }
}
## chr11	132335046	132335286	NTM_ENSG00000182667.14:s:132314552-132315070_132335046-132335286	0.640240740740741	+
## calculate stdev given an array of PSI values
sub stdev{
        my($data) = @_;
        if(@$data == 1){

                return "1";
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

## calculate avg given an array of PSI values
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
__DATA__
perl scripts/extract_recurrent_dominant_lsvs.pl | grep "*" | awk -F "\t" '{if ($5 == 47){ print $0 }}' | awk '{print $6"\t"$1"\t"$2"\t"$7}'| perl -pe 's/(chr[\d+XY]+)\:/$1\t/g' | awk -F "-" '{print $1"\t"$2"-"$3"-"$4}' | awk '{if($6 == "+"){ print $0 } else {print $0"-" }}' | perl -pe 's/\*//' >
