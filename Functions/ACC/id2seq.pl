#!/usr/bin/perl
use strict;

# Function: Given a id file, extracting sequences in another file.
# Author  : liaochenlanruo
# Date	: 2014-08-02 21:08
# Update  : 2022-04-06 16:33
# Usage   : perl $0 <id file> <input aa file> <output aa file>

my $EXE = "id2seq.pl";
my(@Options, $help, $ids, $seqin, $seqout);
setOptions();

# Option setting routines

sub setOptions {
  use Getopt::Long;
  @Options = (
    "GENERAL",
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"ids=s",  VAR=>\$ids, DESC=>"Specify the file containing the IDs of sequences, if the file has multiple columns, space or tab should be used to separate the columns."},
    {OPT=>"seqin=s",  VAR=>\$seqin, DESC=>"Specify the FASTA format file that contains the sequence."},
    {OPT=>"seqout=s",  VAR=>\$seqout, DESC=>"Specify the name of the output file that will be used to save the extracted sequences according to the user-supplied IDs."},
  );

  @ARGV or usage(1);

  &GetOptions(map {$_->{OPT}, $_->{VAR}} (grep { ref } @Options)) || usage(1);
# Now setup default values.
  foreach (grep { ref } @Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  my($exitcode) = @_;
  $exitcode = 0 if $exitcode eq 'help'; # what gets passed by getopt func ref
  $exitcode ||= 0;
  select STDERR if $exitcode; # write to STDERR if exitcode is error

  print "SYNOPSIS\n  Extract the corresponding sequences from a file to another file according to a number of ids in one file.\n";
  print "AUTHOR\n  Hualin Liu\n";
  print "USAGE\n";
  print "    $EXE --ids <IDs file> --seqin <FASTA file> --seqout <OutPut file>\n";

  foreach (@Options) {
    if (!ref $_) { print $_,"\n"; next }  # print text lines
    my $opt = $_->{OPT};
    $opt =~ s/!$//;
    $opt =~ s/=s$/ [X]/;
    $opt =~ s/=i$/ [N]/;
    $opt =~ s/=f$/ [n.n]/;
    printf "  --%-13s %s%s.\n",$opt,$_->{DESC},
           defined($_->{DEFAULT}) ? " [$_->{DEFAULT}]" : "";
  }

  exit($exitcode);
}


my $count = 0;
my $hits  = 0;
local $|  = 1;
my %ids_hash;

#-----------------------[reads ids]--------------
open ID,"$ids" || die "Can't open ID file";
while(<ID>) {
	chomp;
	my @lines = split;
	foreach  (@lines) {
		$ids_hash{$_}++;
	}
}
close ID;

# ---------show number of ids---------
my @ids = keys %ids_hash;
my $num = @ids;
print "\nRead $num ids.\n\n";

#-------------[ searching ]-------------
open IN, "$seqin" || die "Failed to open IN file\n";
open OUT, ">$seqout" || die "Failed to open OUT file\n";

local $/ = ">";
<IN>;
my ( $head, $seq );
while (<IN>) {
	$count++;
	s/\r?\n>//;
	( $head, $seq ) = split "\n", $_, 2;
	$seq =~ s/\s+//g;
	next unless $head =~/(\S+)/;
	if (exists $ids_hash{$1}){
		print OUT ">$head\n$seq\n" x $ids_hash{$1};
		$hits++;
	}
	# record hit number of a id
	$ids_hash{$1}++;
}

print "\rProcessing ${count} th record. hits: $hits";

# Displays IDs that do not match any records
my @idsm = grep {$ids_hash{$_} == 0} keys %ids_hash;
my $num = @idsm;
print "\n\n$num ids did not match any record in $ids\n";
print "@idsm\n";
$/ = "\n";

close IN;
close OUT;
