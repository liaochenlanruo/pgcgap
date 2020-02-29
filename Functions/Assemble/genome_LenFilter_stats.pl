#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use File::Basename;
use List::Util qw(sum min max);
use Getopt::Long;
use Pod::Usage;
use Getopt::Std;

my %options;

=head1 USAGE

  $perl genome_LenFilter_stats.pl --Scaf_suffix <STRING> --filter_length <INT>

=head1 OPTIONS

=over 30

=item B<[--help]>

Print the help message and exit

=back

=cut

$options{'help|h|?'} = \( my $opt_help );

=over 30

=item B<[--Scaf_suffix (STRING)]>

The suffix of scaffolds or genomes [Required by "All", "Annotate", "MASH", "ANI" and "AntiRes"] Here, "-8.fa" for Illumina data, ".contigs.fasta" for PacBio data and Oxford data. Users can also fill in other suffixes according to the actual situation ( Default -8.fa )

=back

=cut

$options{'Scaf_suffix=s'} = \( my $opt_Scaf_suffix = "-8.fa" );

=over 30

=item B<[--filter_length (INT)]>

I<[Required]> Sequences shorter than the 'filter_length' will be deleted from the assembled genomes. ( Default 200 )

=back

=cut

$options{'filter_length=i'} = \( my $opt_filter_length = 200 );


GetOptions(%options) or pod2usage(1);

pod2usage(1) if ($opt_help);


my @fas = glob("*$opt_Scaf_suffix");
		foreach my $fas (@fas) {
			$fas=~/(\S+)$opt_Scaf_suffix/;
			my $outpre = $1 . ".prefilter.stats";
			my $outfilter = $1 . ".filtered.stats";
			my $filterscaf = $1. ".filtered.fas";
			getstats($fas, $outpre);
			lenfilter($fas, $filterscaf, $opt_filter_length);
			getstats($filterscaf, $outfilter);
		}

{
	my $As = 0;
	my $Ts = 0;
	my $Gs = 0;
	my $Cs = 0;
	my $Ns = 0;
	sub getstats {
		my $file = shift;
		my $outFile = shift;

		my ($fileName, $filePath) = fileparse($file);
		$outFile = $file . "_n50_stat" if($outFile eq "");

		open(I, "<$file") or die "Can not open file: $file\n";
		open(O, ">$outFile") or die "Can not open file: $outFile\n";

		my @len = ();

		my $prevFastaSeqId = "";
		my $fastaSeqId = "";
		my $fastaSeq = "";

		while(my $line = <I>) {
			chomp $line;
			if($line =~ /^>/) {
				$prevFastaSeqId = $fastaSeqId;
				$fastaSeqId = $line;
				if($fastaSeq ne "") {
					push(@len, length $fastaSeq);
					baseCount($fastaSeq);
				}
				$fastaSeq = "";
			}
			else {
				$fastaSeq .= $line;
			}
		}
		if($fastaSeq ne "") {
			$prevFastaSeqId = $fastaSeqId;
			push(@len, length $fastaSeq);
			baseCount($fastaSeq);
		}

		my $totalReads = scalar @len;
		my $bases = sum(@len);
		my $minReadLen = min(@len);
		my $maxReadLen = max(@len);
		my $avgReadLen = sprintf "%0.2f", $bases/$totalReads;
		my $medianLen = calcMedian(@len);
		my $n25 = calcN50(\@len, 25);
		my $n50 = calcN50(\@len, 50);
		my $n75 = calcN50(\@len, 75);
		my $n90 = calcN50(\@len, 90);
		my $n95 = calcN50(\@len, 95);

		printf O "%-25s %d\n" , "Total sequences", $totalReads;
		printf O "%-25s %d\n" , "Total bases", $bases;
		printf O "%-25s %d\n" , "Min sequence length", $minReadLen;
		printf O "%-25s %d\n" , "Max sequence length", $maxReadLen;
		printf O "%-25s %0.2f\n", "Average sequence length", $avgReadLen;
		printf O "%-25s %0.2f\n", "Median sequence length", $medianLen;
		printf O "%-25s %d\n", "N25 length", $n25;
		printf O "%-25s %d\n", "N50 length", $n50;
		printf O "%-25s %d\n", "N75 length", $n75;
		printf O "%-25s %d\n", "N90 length", $n90;
		printf O "%-25s %d\n", "N95 length", $n95;
		printf O "%-25s %0.2f %s\n", "As", $As/$bases*100, "%";
		printf O "%-25s %0.2f %s\n", "Ts", $Ts/$bases*100, "%";
		printf O "%-25s %0.2f %s\n", "Gs", $Gs/$bases*100, "%";
		printf O "%-25s %0.2f %s\n", "Cs", $Cs/$bases*100, "%";
		printf O "%-25s %0.2f %s\n", "(A + T)s", ($As+$Ts)/$bases*100, "%";
		printf O "%-25s %0.2f %s\n", "(G + C)s", ($Gs+$Cs)/$bases*100, "%";
		printf O "%-25s %0.2f %s\n", "Ns", $Ns/$bases*100, "%";

		print "N50 Statisitcs file: $outFile\n";

		$As = 0;
		$Ts = 0;
		$Gs = 0;
		$Cs = 0;
		$Ns = 0;

		sub calcN50 {
			my @x = @{$_[0]};
			my $n = $_[1];
			@x=sort{$b<=>$a} @x;
			my $total = sum(@x);
			my ($count, $n50)=(0,0);
			for (my $j=0; $j<@x; $j++){
				$count+=$x[$j];
				if(($count>=$total*$n/100)){
					$n50=$x[$j];
					last;
				}
			}
			return $n50;
		}

		sub calcMedian {
			my @arr = @_;
			my @sArr = sort{$a<=>$b} @arr;
			my $arrLen = @arr;
			my $median;
			if($arrLen % 2 == 0) {
				$median = ($sArr[$arrLen/2-1] + $sArr[$arrLen/2])/2;
			}
			else {
				$median = $sArr[$arrLen/2];
			}
			return $median;
		}

		sub baseCount {
			my $seq = $_[0];
			my $tAs += $seq =~ s/A/A/gi;
			my $tTs += $seq =~ s/T/T/gi;
			my $tGs += $seq =~ s/G/G/gi;
			my $tCs += $seq =~ s/C/C/gi;
			$Ns += (length $seq) - $tAs - $tTs - $tGs - $tCs;
			$As += $tAs;
			$Ts += $tTs;
			$Gs += $tGs;
			$Cs += $tCs;
		}
	}
}

sub lenfilter{
	my $scaf = shift;
	my $scaf_filter = shift;
	my $length = shift;
	my $in = Bio::SeqIO->new(-file=>"$scaf", -format=>"fasta");
	my $out = Bio::SeqIO->new(-file=>">$scaf_filter", -format=>"fasta");
	while (my $seq = $in->next_seq) {
		my $len = $seq->length;
		if ($len >= $length) {
			$out->write_seq($seq);
		}
	}
}
