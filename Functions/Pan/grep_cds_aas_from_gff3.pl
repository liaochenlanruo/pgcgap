#!/usr/bin/perl
use strict;
use warnings;

=head1 #######################################
	
=head1 name
	grep_cds_aas_from_gff3.pl

=head1 description
	Getting CDs and AAS sequences from gff3 file.

=head1 example
	perl $0 gff_file file_prefix

=head1 author
	Hualin Liu
		#2021-4-03

=head1 #######################################
=cut

die `pod2text $0` unless @ARGV==2;
my ($gff,$prefix)=@ARGV;

# deal gff file
open GFF, "$ARGV[0]" || die;
my (%cds,%fasta);
local $/=">";
my $line = 0;
while(<GFF>){
	chomp;
	$line++;
	if ($line == 1) {
		my @lines = split /\n/;
		for (my $i=0; $i<@lines; $i++) {
			next if($lines[$i]=~/^#/);
			my @p=split /\t/,$lines[$i];
			my @q=split/;/,$p[8];
			my $pep;
			my $chr=$p[0];
			if($p[2] eq "CDS"){
				($pep=$q[0])=~s/ID=//;
				$cds{$pep} = [$chr, $p[3], $p[4], $p[6]];
			}
		}
	}else {# process the scaffolds
		my ($id, $seq) = split (/\n/, $_, 2);
		$id =~ /^(\S+)/;
		my $index= $1;
		$seq =~ s/\s+|\n|\-//g;
		$fasta{$index} = $seq;
	}
}
close GFF;

$/="\n";
# get id list
my $id_gene_cds_pep=$prefix.".cds.id";
open ID, ">",$id_gene_cds_pep||die;
foreach my $i(sort keys %cds){
	if($cds{$i}){
		my $out=join "\t",$i,$cds{$i}[0],$cds{$i}[1],$cds{$i}[2],$cds{$i}[3];
		print ID $out . "\n";
	}
}
warn "create file:$id_gene_cds_pep\n";
close ID;

# Get CDs sequences
open LIST, "$id_gene_cds_pep" || die;
my %gene;
while (<LIST>) {
	chomp;
	my @lines = split;
	if (exists $fasta{$lines[1]}) {
		my $cut = substr($fasta{$lines[1]},$lines[2]-1,$lines[3]-$lines[2]+1);
		if ($lines[4] eq "+") {
			$gene{$lines[0]} = $cut;
		}elsif ($lines[4] eq "-") {
			$cut = &reverse_complement($cut);
			$gene{$lines[0]} = $cut;
		}
	}else {
		print "No scaffold was found for $lines[1]\n";
	}
}


# Print CDs sequences
my $cds=$prefix.".cds";
warn "create file:$cds\n";
open CDS, ">$cds" || die;
foreach  (sort keys %gene) {
	print CDS ">$_\n$gene{$_} \n\n";
}
close CDS;

# Print AAs sequences
my $raw_pep=$prefix.".pep";
open PEP1, ">$raw_pep" || die;
my @raw_pep= &cds2pep($cds);
foreach  (@raw_pep) {
	print PEP1 $_ . "\n";
}
#print PEP1 $_ for@raw_pep;
close PEP1;
warn "create file:$raw_pep\n";


# subroutine
sub reverse_complement{
	my ($seq)=shift;
	$seq=reverse$seq;
	$seq=~tr/AaGgCcTt/TtCcGgAa/;
	return $seq;
}

# subroutine
sub cds2pep{
	my %code = (
						"standard" =>
							{
								'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
								'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
								'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Aci
								'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Aci
								'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanin
								'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
								'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
								'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
								'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
								'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
								'ATG' => 'M',                                                                         # Methionine
								'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
								'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
								'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
								'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
								'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
								'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
								'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
								'TGG' => 'W',                                                                         # Tryptophan
								'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
								'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U'                                              # Stop
							},
						"code11" =>
							{
								'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
								'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
								'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Aci
								'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Aci
								'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanin
								'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
								'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
								'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
								'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
								'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
								'ATG' => 'M',                                                                         # Methionine
								'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
								'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
								'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
								'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
								'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
								'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
								'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
								'TGG' => 'W',                                                                         # Tryptophan
								'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
								'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U'                                              # Stop
							}
							## more translate table could be added here in future
	);
	local $/=">";
	my $file=shift;
	open IN, $file || die;
	<IN>;
	my @results_set;
	while(<IN>){
		chomp;
		my ($id, $seq) = split (/\n/, $_);
		my @a=split /\s+/, $id;
		$seq =~ s/\s+|\n|\-//g;
#		$seq=~s/\n|>//g;
		my $len=length $seq;# length of sequence
		my $info_out= $a[0];# sequence id
		my ($pep_out, $triplet);
		for(my $i=0;$i<$len;$i+=3){
			$triplet=substr($seq,$i,3);
			next if (length $triplet != 3);
			if (exists $code{code11}{$triplet}){
				$pep_out .= $code{code11}{$triplet};# protein sequence
			}else{
				$pep_out .= "X";
			}
		}
		$pep_out=~s/U$// if ($pep_out=~/U$/);
		my $pep_len = length $pep_out;
		$pep_out=~s/([A-Z]{50})/$1\n/g;# split to 50 characters per line
		chop($pep_out) unless($pep_len%50);
		my $results= ">$info_out\n$pep_out\n";
		push @results_set, $results;
	}
	return @results_set;
	$/="\n";
}

__END__

