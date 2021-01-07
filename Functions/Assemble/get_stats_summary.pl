#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Getopt::Std;

my %options;

=head1 USAGE

  $perl get_stats_summary.pl --Scaf_suffix <STRING>

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


GetOptions(%options) or pod2usage(1);

pod2usage(1) if ($opt_help);


open OUT, ">Summary.txt" || die;
print OUT "Strain name\tTotal sequences filtered\tTotal bases filtered\tN50 filtered\tGC % filtered\tTotal sequences prefilter\tTotal bases prefilter\tN50 prefilter\tGC % prefilter\n";
my %hash;
my %good;
my $str;
my @stats = glob("*.stats");
foreach  (sort @stats) {
	if (/(.+)\.filtered.stats/) {
		$str = $1;
#		$hash{$str}++;
		print OUT $str;
		open IN, "$_" || die;
		while (<IN>) {
			chomp;
			if (/Total sequences\s+(\d+)/) {
				print OUT "\t$1";
			}
			if (/Total bases\s+(\d+)/) {
				print OUT "\t$1";
			}
			if (/N50 length\s+(\d+)/) {
				print OUT "\t$1";
				if ($1 < 50000) {
					$hash{$str}++;
				}else {
					$good{$str}++;
				}
			}
			if (/\(G \+ C\)s\s+(\S+)/) {
				print OUT "\t$1";
			}
		}
		close IN;
	}elsif (/(.+)\.prefilter.stats/) {
		if ($1 eq $str) {
			open IN, "$_" || die;
			while (<IN>) {
				chomp;
				if (/Total sequences\s+(\d+)/) {
					print OUT "\t$1";
				}
				if (/Total bases\s+(\d+)/) {
					print OUT "\t$1";
				}
				if (/N50 length\s+(\d+)/) {
					print OUT "\t$1";
				}
				if (/\(G \+ C\)s\s+(\S+)/) {
					print OUT "\t$1";
				}
			}
			print OUT "\n";
			close IN;
		}
	}
}
close OUT;

system("mkdir -p StatsFiles");
system("mv *.stats StatsFiles");

my @low = keys %hash;
if (@low) {
	system("mkdir -p N50lt50K");
	foreach  (@low) {
		my $genome = $_ . $opt_Scaf_suffix;
		my $filteredgenome = $_ . ".filtered.fas";
		system("mv $genome $filteredgenome N50lt50K");
	}
}

system("mkdir -p Filtered");
system("mkdir -p Prefilter");
system("mv *.filtered.fas Filtered");


my @perfect = keys %good;
if (@perfect) {
	foreach  (@perfect) {
		my $genome = $_ . $opt_Scaf_suffix;
		system("mv $genome Prefilter");
	}
}
