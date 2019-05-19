#!/usr/bin/env perl
use strict;
use warnings;


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
	if ($lines[0]=~/.+\/(\S+)-8.fa/) {
		$scafl{$1} = 1;
		$l = $1;
	}elsif ($lines[0]=~/.+\/(\S+)/) {
		$scafl{$1} = 1;
		$l = $1;
	}

	if ($lines[1]=~/.+\/(\S+)-8.fa/) {
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
		}
	}
	print OUT "\n";
}
close OUT;
