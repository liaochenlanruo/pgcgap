#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Getopt::Std;

my %options;

=head1 USAGE

  $ perl get_ANImatrix.pl

=head1 OPTIONS

=over 30

=item B<[--help]>

Print the help message and exit

=back

=cut

$options{'help|h|?'} = \( my $opt_help );

=over 30

=item B<[--Scaf_suffix (STRING)]>

The suffix of scaffolds or genomes ( Default "-8.fa" )

=back

=cut

$options{'Scaf_suffix=s'} = \( my $opt_Scaf_suffix = "-8.fa" );


GetOptions(%options) or pod2usage(1);

pod2usage(1) if ($opt_help);



my %hash;
my %scafl;
my %scafr;

open IN, "ANIs" || die;
open OUT,">ANIs.heatmap" || die;

while (<IN>) {
	chomp;
	my $l;
	my $r;
	my @lines = split /\t/;
	if ($lines[0]=~/.+\/(\S+)$opt_Scaf_suffix/) {
		$scafl{$1} = 1;
		$l = $1;
	}elsif ($lines[0]=~/.+\/(\S+)/) {
		$scafl{$1} = 1;
		$l = $1;
	}

	if ($lines[1]=~/.+\/(\S+)$opt_Scaf_suffix/) {
		$scafr{$1} = 1;
		$r = $1;
	}elsif ($lines[0]=~/.+\/(\S+)/) {
		$scafr{$1} = 1;
		$r = $1;
	}
	$hash{$l}{$r} = $lines[2];
}
close IN;

my @scafl = sort keys %scafl;
my @scafr = sort keys %scafr;
#my @scaflr = reverse @scafl;

print OUT "str";
foreach  (@scafl) {
	print OUT "\t" . $_;
}
print OUT "\n";


for (my $i=0; $i<@scafr; $i++) {
	print OUT $scafr[$i];
	for (my $j=0; $j<@scafl; $j++) {
		if (exists $hash{$scafr[$i]}{$scafl[$j]}) {
			print OUT "\t" . $hash{$scafr[$i]}{$scafl[$j]};
		}else{
			print OUT "\t" . "0";
		}
	}
	print OUT "\n";
}
close OUT;
