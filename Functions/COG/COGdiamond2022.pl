#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Getopt::Std;
use Parallel::ForkManager;

my %options;

=head1 USAGE

  $perl COGdiamond2022.pl --threads <INT>

=head1 OPTIONS

=over 30

=item B<[--help]>

Print the help message and exit

=back

=cut

$options{'help|h|?'} = \( my $opt_help );

=over 30

=item B<[--threads (INT)]>

Number of threads to be used ( Default 4 )

=back

=cut

$options{'threads=i'} = \( my $opt_threads = 4 );

=over 30

=item B<[--AAsPath (PATH)]>

I<[Required]> Amino acids of all strains as fasta file paths, ( Default "./Results/Annotations/AAs" )

=back

=cut

$options{'AAsPath=s'} = \( my $opt_AAsPath = "./Results/Annotations/AAs" );

=over 30

=item B<[--evalue (FLOAT)]>

I<[Required]> maximum e-value to report alignments, ( Default 1e-3 )

=back

=cut

$options{'evalue=f'} = \( my $opt_evalue = "1e-3" );

=over 30

=item B<[--id (INT)]>

I<[Required]> minimum identity% to report an alignment, ( Default 40 )

=back

=cut

$options{'id=i'} = \( my $opt_id = 40 );

=over 30

=item B<[--query_cover (INT)]>

I<[Required]> minimum query cover% to report an alignment, ( Default 70 )

=back

=cut

$options{'query_cover=i'} = \( my $opt_query_cover = 70 );

=over 30

=item B<[--subject_cover (INT)]>

I<[Required]> minimum subject cover% to report an alignment, ( Default 50 )

=back

=cut

$options{'subject_cover=i'} = \( my $opt_subject_cover = 50 );

=over 30

=item B<[--strain_num (INT)]>

I<[Required]> The total number of strains used for analysis

=back

=cut

$options{'strain_num=i'} = \( my $opt_strain_num );

GetOptions(%options) or pod2usage(1);

pod2usage(1) if ($opt_help);


#=================Get the full path of COG database===================================
my $cogdb_dir;
my $path = `which COGdiamond_2020.dmnd`;
if ($path=~/(.+)\/COGdiamond_2020.dmnd/) {
	$cogdb_dir = $1;
}

#=============================== Get bin PATH ======================================================
my $COG_dir;
my $bin = `which COGdiamond2022.pl`;
if ($bin=~/(.+)\/COGdiamond2022.pl/) {
	$COG_dir = $1;
}

chdir $opt_AAsPath;
my @faa = glob("*.faa");
my $pm = new Parallel::ForkManager($opt_threads);
foreach (@faa) {
	$pm->start and next;
	&run_COG;
	$pm->finish;
}
$pm->wait_all_children;


#============================Get the relative abundance table==========================
if ($opt_strain_num > 1) {
	system("perl $COG_dir/get_flag_relative_abundances_table.pl");
	system("Rscript $COG_dir/Plot_COG_Abundance.R");
}

sub rm_duplicate {
	my $array_ref = shift;
	my %hash;
	foreach  ( @{$array_ref}) {
		$hash{$_}++;
	}
	return keys %hash;
}


sub run_COG {
	#my $faa = shift;
	$_ =~ /(.+).faa/;
	my $name = $1;
	my $blastout = $name . ".diamond.COG.xml"; 
	my $cog_gi = $name . ".2gi.table";
	my $cog_ids = $name . ".2id.table";
	my $super_id = $name . ".2Sid.table";
	my $super_cog = $name . ".2Scog.table";

#============================Run blastp===============================================

	my $threads;
	if ($opt_strain_num < $opt_threads) {
		$threads = int($opt_threads/$opt_strain_num);
	}else {
		$threads = 1;
	}
	system("diamond blastp --db $cogdb_dir/COGdiamond_2020.dmnd --query $_ --out $blastout --evalue $opt_evalue --outfmt 5 --max-target-seqs 1 --threads $threads --id $opt_id --query-cover $opt_query_cover --subject-cover $opt_subject_cover");

#======================================blast2gi========================================

	open CGI, ">$cog_gi" || die;

	my $start_tag = '^<Iteration>$';
	my $ref = &Tag_Reader($blastout,$start_tag);
	my $content;
	while($content = &$ref){
		my ($label,@content) = @$content;
		my ($query_def) = map {/Iteration_query-def>(.*)<\/Iteration_query-def/} @content;
		my $nohit_label = grep {/No hits found/} @content;
		if ($nohit_label){
			print CGI "$query_def\n";
			next;
		}else{
			my @gi = map {/.*>(\S+).*?<.*/g} @content;
			my @gi2;
			my $gi2;
			foreach my $gi (@gi) {
				if ($gi=~/(\S+)_(\d+)/) {
					$gi2 = $1 . "." . $2;
					push @gi2, $gi2;
				}else {
					push @gi2, $gi;
				}
			}
			print CGI "$query_def\t",(join "\t",@gi2),"\n";
		}
	}

	sub Tag_Reader {
		my $file = shift;
		my $start_tag = shift;
		die "Error:$file unreachable.\n" unless (-s $file);
		open IN, $file || die;
		my $temp1;
		while(<IN>){
			next unless /$start_tag/;
			$temp1 = $_;
			last;
		}

		my $temp2 = $temp1;
		my @temp;
		return sub {
			if (eof(IN)){
				close IN;
				return '';
			}
			while(!eof(IN)){
				$temp1 = $temp2;
				@temp=();
				while(<IN>){
					if (/$start_tag/){
						$temp2 = $_;
						last;
					}
				push @temp,$_;
				}
				return [$temp1,@temp];
			}
		}
	}
	close CGI;

#===============================================gi2id==================================
	open GI, $cog_gi || die;
	open CID, ">$cog_ids" || die;

	my %genes;
	while (<GI>){
		chomp;
		my ($gene,@cog_gi) = split /\s+/,$_;
		next unless (@cog_gi);
		foreach (@cog_gi){
			push @{$genes{$_}},$gene;
		}
	}
	close GI;

	open CSV, "$cogdb_dir/cog-20.cog.csv" || die;

	my %gene2cog_id;
	while(<CSV>){
		chomp;
		my @F = split /,/,$_;
		if($genes{$F[2]}){
			foreach (@{$genes{$F[2]}}){
				push @{$gene2cog_id{$_}},$F[6];
			}
		}
	}
	close CSV;

	foreach (sort keys %gene2cog_id){
		my $gene = $_;
		my @cog_ids = &rm_duplicate(\@{$gene2cog_id{$gene}});
		print CID ("$gene\t",join ("\t",@cog_ids),"\n");
	}
	close CID;

#=============================================id2Sid===================================
	my %cog2query;
	open QUERY, $cog_ids || die;
	open SID, ">$super_id" || die;
	while(<QUERY>){
		chomp;
		my ($query,@cog_id) = split /\s+/, $_;
		foreach (@cog_id){
			push @{$cog2query{$_}},$query;
		}
	}
	close QUERY;

	my %gi2final_cog;
	open NAME, "$cogdb_dir/cog-20.def.tab" || die;
	while (<NAME>){
		chomp;
		my ($cog_id,$flag,$name) = split /\t/,$_,3;
		if ($cog2query{$cog_id}){
			foreach (@{$cog2query{$cog_id}}){
				push @{$gi2final_cog{$_}},"$flag\t$name";
			}
		}
	}
	close NAME;

	for (keys %gi2final_cog){
		my $query = $_;
		foreach (&rm_duplicate(\@{$gi2final_cog{$query}})){
			print SID "$query\t$_\n";
		}
	}
	close SID;
	

#==============================================Sid2Scog================================
	
	open TABLE, $super_id || die;
	open SOG, ">$super_cog" || die;

	my %flag2queries;
	while (<TABLE>){
		chomp;
		my ($query,$flags)=(split /\s+/,$_)[0,1];
		my @flags = split //,$flags;
		for (@flags){
			push @{$flag2queries{$_}},$query;
		}
	}
	close TABLE;

	#open FUN, "$cogdb_dir/fun2003-2014.tab" || die;
	open FUN, "$cogdb_dir/fun-20.tab" || die;
	my %query2FLAG_name;
	while(<FUN>){
		chomp;
		my @F = (split /\t/,$_,3);
		if($flag2queries{$F[0]}){
			for (@{$flag2queries{$F[0]}}){
				push @{$query2FLAG_name{$_}},"$F[0]\t[$F[0]] $F[2]";
			}
		}
	}
	close CSV;

	my %FLAG_name2query;
	for (keys %query2FLAG_name){
		my $query = $_;
		for(&rm_duplicate(\@{$query2FLAG_name{$query}})){
			push @{$FLAG_name2query{$_}},$query;
		}
	}
	for (sort {$a cmp $b} keys %FLAG_name2query){
		my $FLAG_name = $_;
		for (&rm_duplicate(\@{$FLAG_name2query{$FLAG_name}})){
			print SOG "$_\t$FLAG_name\n";
		}
	}
	close SOG;


#=================================Plot=================================================
	system("Rscript $COG_dir/Plot_COG.R $super_cog");
}