#!/usr/bin/perl
use strict;
use warnings;

my $EXE = "get_repeats.pl";
my(@Options, $help, $filein, $column, $sep);
setOptions();

# Option setting routines
sub setOptions {
  use Getopt::Long;
  @Options = (
    "GENERAL",
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"filein=s",  VAR=>\$filein, DESC=>"Specify the input file."},
    {OPT=>"column=i",  VAR=>\$column, DEFAULT=>0, DESC=>"Specifies which column is used for the calculation (default: 0 for the whole line)"},
    {OPT=>"sep=s",  VAR=>\$sep, DEFAULT=>"tab", DESC=>"Specify the separator (space, tab, comma, semicolon) between columns (Default: tab)."},

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

  print "SYNOPSIS\n  Counts the number of repeats of the string for the specified column in a given file.\n";
  print "AUTHOR\n  Hualin Liu\n";
  print "USAGE\n";
  print "    $EXE --filein <Input file> --column <INT> --sep <space, tab, comma>\n";

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

my %hash;
open IN, $filein || die;
my $out = $filein . ".num";
open OUT, ">$out" || die;

while (<IN>) {
	chomp;
	$_=~ s/[\r\n]+//g;
	if ($column == 0) {
		$hash{$_}++;
	}else {
		my @lines;
		if ($sep eq "space") {
			@lines = split /\s+/;
		}elsif ($sep eq "tab") {
			@lines = split /\t/;
		}elsif ($sep eq "comma") {
			@lines = split /,/;
		}elsif ($sep eq "semicolon") {
			@lines = split /\;/;
		}
		$hash{$lines[$column - 1]}++;
	}
}
close IN;

foreach  (keys %hash) {
	print OUT "$_\t$hash{$_}\n";
}
close OUT;

print "The result is saved in file $out, whose last column shows the number of repetitions.\n";