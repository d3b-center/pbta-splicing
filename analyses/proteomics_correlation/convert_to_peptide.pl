#!/usr/bin/perl

use strict;
use LWP::UserAgent;
use warnings;

$ENV{PERL_LWP_SSL_VERIFY_HOSTNAME} = 0;


my %exon2seq = ();
my $id = '';

my $fasta_file = $ARGV[0];
my $exon2seq;
my %strand;
my %event_splice_id;

open(FIL,"results/mixed_events.tsv") || die("Cannot Open File");
while(<FIL>)
{
	chomp;
	my @cols = split "\t";
	my $splice_id = $cols[5];
	my $chr = $cols[1];
	my $start = $cols[2];
	my $end = $cols[3];
	my $str = $cols[4];

	my $event = $chr.":".$start."-".$end;
	$strand{$event} = $str;
	$event_splice_id{$event} = $splice_id;
	#print $event,"X\tX",$str,"X\n";

}
open (FASTA,$fasta_file)|| die("Cannot Open File $fasta_file");
while(<FASTA>){
		chomp;
		#print $_,"\n";
		if($_ =~ /^>(.+)/){
				$id = $1;

		}else{
				$exon2seq{$id} .= $_;
				#print $id,"\t",$_,"\n";
		}
}
close (FASTA);

open(OUT,">results/splice_events.peptides.test.fa");

foreach my $event(keys %exon2seq)
{

	my $browser  = LWP::UserAgent->new;
	my $rna_sequence = $exon2seq{$event};
	my $response = $browser->post(
		'https://web.expasy.org/cgi-bin/translate/dna2aa.cgi',
		[
			'dna_sequence'    => $rna_sequence,
			'output_format'   => 'fasta'
		]);

	my $res = $response->content;
	#print $event,"\t",$res,"\n";
	#my @res;
	my @res = split/\n/,$res;
	my @frame_seqs;
	foreach my $entry(@res)
	{
		next if ($entry=~/\>/);
		next unless ($entry=~/\w+/);
		push @frame_seqs,$entry;
	}

	if($strand{$event}=~/\-/){ # 3-5'
		print OUT ">".$event_splice_id{$event},"\n",$frame_seqs[0],"\n";
		print OUT ">".$event_splice_id{$event},"\n",$frame_seqs[1],"\n";
		print OUT ">".$event_splice_id{$event},"\n",$frame_seqs[2],"\n";
	}

	if($strand{$event}=~/\+/){ # 5-3'
		print OUT ">".$event_splice_id{$event},"\n",$frame_seqs[3],"\n";
		print OUT ">".$event_splice_id{$event},"\n",$frame_seqs[4],"\n";
		print OUT ">".$event_splice_id{$event},"\n",$frame_seqs[5],"\n";
	}

}
