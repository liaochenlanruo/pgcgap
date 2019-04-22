#!/usr/bin/env perl
use strict;
use warnings;

#Date     : 2019.04.17
#Author   : Hualin Liu
#Function : Calculate the relative abundance of the flags and form a table



#INPUT Files: *.2Scog.table

#OUTPUT: All_flag_abundances.table

#---------------------------Calculate the relative abundance of the flags-------------------

my @super_cog = glob("*.2Scog.table");

foreach (@super_cog){

	$_ =~ /(\S+).2Scog.table/;

    my $name = $1;

	my $out = $name . "_standard.table";

    my ($lineCount, $totalWords, %map);

    open IN,"$_" || die;

    open OUT,">$out" || die;
    while ( <IN> ) {
       chomp;
       $lineCount++;
	   my @row = split /\t/;
       $row[1] =~ s/\S+//g;
       my $wordPos;
	   
my ($id,$flag,$ano) = split /\t/,$_,3;
       $wordPos++;
       $totalWords++;
       $map{$flag} = [] unless ( defined $map{$flag} );
       my $this = $map{$flag};                
       push @{$this}, "$lineCount-$wordPos";    
     }
    print OUT "There are total $lineCount lines in this context\n";

    foreach my $word ( sort keys %map ) {
       my $count = scalar ( @{$map{$word}} );
       my $stand = $count/$lineCount;

	   print OUT $word."\t".$stand."\n";
     }

}
close IN;
close OUT;



#-----------------------------make all the "*_standard.table" files to single file-----------------

my @standard = glob("*_standard.table");
foreach  (@standard){

	$_ =~ /(\S+)_standard.table/;
	my $name = $1;
	open IN,"$_" || die;
	
readline IN;
	open OUT,">>All_cogs.txt" || die;
	while (<IN>) {
	  chomp;
	  print OUT $name."\t$_\n";
	}
}
close IN; 
close OUT;




#--------formation a table of the flags from a single file which contend ids and flags，show the relative abundance of the flags----------------------"

my %hash;
my %ids;
my %flag;
open IN, "All_cogs.txt" || die;
open OUT, ">All_flags_relative_abundances.table" || die;

while (<IN>) {
	if (/(\S+)\t(\S+)\t(\S+)/) {
        my $ids = $1;
		$ids{$1} = "1";
		my $flag = $2;
		$flag{$2} = "1";
		$hash{$1}{$2} = "$3";
	}
}

my @ids = sort keys %ids;
my @flag = sort keys %flag;
print OUT "strains";
foreach  (@flag) {
	print OUT "\t" . $_;
}

print OUT "\n";
for (my $i = 0;$i < @ids;$i++) {
	print OUT $ids[$i];
	for (my $j = 0;$j < @flag;$j++) {
		if (exists $hash{$ids[$i]}{$flag[$j]}) {
		  print OUT "\t" . "$hash{$ids[$i]}{$flag[$j]}";
		}else{
			print OUT "\t" . "0";
		}
	}
	print OUT "\n";
}
close IN;
close OUT;
system("rm *_standard.table All_cogs.txt");

