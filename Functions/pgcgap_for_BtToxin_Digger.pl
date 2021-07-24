#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;#basic
use Pod::Usage;#acc
use Getopt::Std;#basic
use Bio::SeqIO;#bioperl
#use Data::Dumper;
use File::Tee qw(tee);#acc
use Cwd;#basic
use List::Util qw(sum min max);#acc
use File::Basename;#basic
#use POSIX;#acc
##use Sys::Info;#acc
##use Sys::Info::Constants qw( :device_cpu );
##use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);#basic

my %options;


=head1 NAME

PGCGAP

=head1 DESCRIPTION

The prokaryotic genomics and comparative genomics analysis pipeline

=head1 AUTHOR

Hualin Liu

=head1 CONTACT

liaochenlanruo@webmail.hzau.edu.cn

=head1 USAGE

  General usage: pgcgap [Modules] [Options]

  Show parameters for each module: pgcgap [Assemble|ACC]

  Show examples of each module: pgcgap Examples

=head1 OPTIONS

=over 30

=item B<[--help]>

Print the help message and exit

=back

=cut

$options{'help|h|?'} = \( my $opt_help );

=over 30

=item B<[--version]>

Show version number of PGCGAP and exit

=back

=cut

$options{'version'} = \( my $opt_version );

=over 30

=item B<[--check-external-programs]>

Check if all of the required external programs can be found and are executable, then exit

=back

=cut

$options{'check-external-programs'} = \( my $opt_check_external_programs = 0 );


=head2 *********************************************** Modules ************************************************

=for text



=over 30

=item B<[--Assemble]>

Assemble reads (short, long or hybrid) into contigs

=back

=cut

$options{'Assemble'} = \(my $opt_Assemble);

=over 30

=item B<[--ACC]>

Other useful gadgets

=back

=cut

$options{'ACC'} = \(my $opt_ACC);

=head2 *********************************************** Global Options *****************************************

=for text


=over 30

=item B<[--filter_length (INT)]>

I<[Required]> Sequences shorter than the 'filter_length' will be deleted from the assembled genomes [Required by "Assemble" and "Assess"]. ( Default 200 )

=back

=cut

$options{'filter_length=i'} = \(my $opt_filter_length = 200);

=over 30

=item B<[--logs (STRING)]>

Name of the log file ( Default Logs.txt )

=back

=cut

$options{'logs=s'} = \( my $opt_logs = "Logs.txt" );

=over 30

=item B<[--threads (INT)]>

Number of threads to be used ( Default 4 )

=back

=cut

$options{'threads=i'} = \( my $opt_threads = 4 );

=head2 *********************************************** Local Options ******************************************

=for text



=head3 =========================== Options of "--Assemble" for reads assembly ================================

=for text



=begin html

If you use the results of "--Assemble" function in your work, please also cite one of the following:

</br>

</br>Shaun D Jackman, Benjamin P Vandervalk, Hamid Mohamadi, Justin Chu, Sarah Yeo, S Austin Hammond, Golnaz Jahesh, Hamza Khan, Lauren Coombe, Ren&eacute; L Warren, and Inanc Birol (2017). ABySS 2.0: Resource-efficient assembly of large genomes using a Bloom filter. Genome research, 27(5), 768-777. doi:10.1101/gr.214346.116

</br>

</br>Koren S, Walenz BP, Berlin K, Miller JR, Phillippy AM. Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. Genome Research. (2017).

=end html

=over 30

=item B<[--ReadsPath (PATH)]>

I<[Required]> Reads of all strains as file paths ( Default ./Reads/Illumina )

=back

=cut

$options{'ReadsPath=s'} = \( my $opt_ReadsPath = "./Reads/Illumina" );

=over 30

=item B<[--reads1 (STRING)]>

I<[Required]> The suffix name of reads 1 ( for example: if the name of reads 1 is "YBT-1520_L1_I050.R1.clean.fastq.gz", "YBT-1520" is the strain same, so the suffix name should be ".R1.clean.fastq.gz" )

=back

=cut

$options{'reads1=s'} = \(my $opt_reads1);

=over 30

=item B<[--reads2 (STRING)]>

I<[Required]> The suffix name of reads 2( for example: if the name of reads 2 is "YBT-1520_2.fq", the suffix name should be _2.fq" )

=back

=cut

$options{'reads2=s'} = \(my $opt_reads2);

=over 30

=item B<[--suffix_len (INT)]>

I<[Required]> B<(Strongly recommended)> The suffix length of the reads file, that is the length of the reads name minus the length of the strain name. For example the --suffix_len of "YBT-1520_L1_I050.R1.clean.fastq.gz" is 26 ( "YBT-1520" is the strain name ) ( Default 0 )

=back

=cut

$options{'suffix_len=i'} = \(my $opt_suffix_len = 0);

=over 30

=item B<[--platform (STRING)]>

I<[Required]> Sequencing Platform, "illumina", "pacbio", "oxford" and "hybrid" available ( Default illumina )

=back

=cut

$options{'platform=s'} = \(my $opt_platform = "illumina");

=over 30

=item B<[--assembler (STRING)]>

I<[Required]> Software used for illumina reads assembly, "abyss", "spades" and "auto" available ( Default auto )

=back

=cut

$options{'assembler=s'} = \(my $opt_assembler = "auto");

=over 30

=item B<[--kmmer (INT)]>

I<[Required]> k-mer size for genome assembly of Illumina data with abyss( Default 81 )

=back

=cut

$options{'kmmer=i'} = \(my $opt_kmmer = 81);

=over 30

=item B<[--genomeSize (STRING)]>

I<[Required]> An estimate of the size of the genome. Common suffixes are allowed, for example, 3.7m or 2.8g. Needed by PacBio data and Oxford data ( Default Unset )

=back

=cut

$options{'genomeSize=s'} = \(my $opt_genomeSize);

=over 30

=item B<[--short1 (STRING)]>

I<[Required]> FASTQ file of first short reads in each pair. Needed by hybrid assembly ( Default Unset )

=back

=cut

$options{'short1=s'} = \(my $opt_short1);

=over 30

=item B<[--short2 (STRING)]>

I<[Required]> FASTQ file of second short reads in each pair. Needed by hybrid assembly ( Default Unset )

=back

=cut

$options{'short2=s'} = \(my $opt_short2);

=over 30

=item B<[--long (STRING)]>

I<[Required]> FASTQ or FASTA file of long reads. Needed by hybrid assembly ( Default Unset )

=back

=cut

$options{'long=s'} = \(my $opt_long);

=over 30

=item B<[--hout (STRING)]>

I<[Required]> Output directory for hybrid assembly ( Default ../../Results/Assembles/Hybrid )

=back

=cut

$options{'hout=s'} = \(my $opt_hout = '../../Results/Assembles/Hybrid');

=head3 ========================== Options for "--ACC" ========================================================

=over 30

=item B<[--Assess (STRING)]>

Filter short sequences in the genome and assess the status of the genome.

=back

=cut

$options{'Assess'} = \( my $opt_Assess);

=over 30

=item B<[--scafPath (PATH)]>

I<[Required]> Path for contigs/scaffolds ( Default "Results/Assembles/Scaf/Illumina" )

=back

=cut

$options{'scafPath=s'} = \(my $opt_scafPath = "./Results/Assembles/Scaf/Illumina");



=over 30

=item B<[--Scaf_suffix (STRING)]>

The suffix of scaffolds or genome files [Required]. This is an important parameter that must be set ( Default -8.fa )

=back

=cut

$options{'Scaf_suffix=s'} = \( my $opt_Scaf_suffix = "-8.fa" );

=head2 *************************** Paths of external programs *************************************************

=for text



Not needed if they were in the environment variables path. Users can check with the "--check-external-programs" option for the essential programs

=over 30

=item B<[--abyss-bin (PATH)]>

Path to abyss binary file. Default tries if abyss is in PATH;

=back

=cut

$options{'abyss-bin=s'} = \( my $opt_abyss_bin = `which abyss-pe 2>/dev/null` );

=over 30

=item B<[--canu-bin (PATH)]>

Path to canu binary file. Default tries if canu is in PATH;

=back

=cut

$options{'canu-bin=s'} = \( my $opt_canu_bin = `which canu 2>/dev/null` );

=over 30

=item B<[--unicycler-bin (PATH)]>

Path to unicycler binary file. Default tries if unicycler is in PATH;

=back

=cut

$options{'unicycler-bin=s'} = \( my $opt_unicycler_bin = `which unicycler 2>/dev/null` );


=begin text

  ################################### About The Software ###################################

		  ____       ____      ____     ____       _        ____    
		U|  _"\ u U /"___|u U /"___| U /"___|u U  /"\  u  U|  _"\ u 
		\| |_) |/ \| |  _ / \| | u   \| |  _ /  \/ _ \/   \| |_) |/ 
		 |  __/    | |_| |   | |/__   | |_| |   / ___ \    |  __/   
		 |_|        \____|    \____|   \____|  /_/   \_\   |_|      
		 ||>>_      _)(|_    _// \\    _)(|_    \\    >>   ||>>_    
		(__)__)    (__)__)  (__)(__)  (__)__)  (__)  (__) (__)__)   


  Software: PGCGAP - The prokaryotic genomics and comparative genomics analysis pipeline

  Author: Hualin Liu

  Contact: liaochenlanruo@webmail.hzau.edu.cn

  Citation: Liu H, Xin B, Zheng J, Zhong H, Yu Y, Peng D, Sun M. Build a bioinformatic analysis platform and apply it to routine analysis of microbial genomics and comparative genomics. Protocol exchange, 2021. DOI: 10.21203/rs.2.21224/v5

=end text

=cut

if ($opt_Assemble or $opt_ACC) {
	tee STDOUT, ">>$opt_logs";
}

GetOptions(%options) or pod2usage("Try '$0 --help' for more information.");

if($opt_version){
	print "PGCGAP version: 1.0.29\n";
	exit 0;
}

#pod2usage( -verbose => 1 ) if $opt_help;
if ($opt_help) {
	pod2usage(1);
	exit(0);
}

chomp($opt_abyss_bin, $opt_canu_bin, $opt_unicycler_bin);
check_external_programs() if($opt_check_external_programs);

pod2usage( -msg => 'abyss not in $PATH and binary not specified use --abyss-bin', -verbose => 0, -exitval => 1 ) unless ($opt_abyss_bin);
pod2usage( -msg => 'canu not in $PATH and binary not specified use --canu-bin', -verbose => 0, -exitval => 1 ) unless ($opt_canu_bin);
pod2usage( -msg => 'unicycler not in $PATH and binary not specified use --unicycler-bin', -verbose => 0, -exitval => 1 ) unless ($opt_unicycler_bin);



sub check_external_programs{
	my %programs = ("abyss" => $opt_abyss_bin, "canu" => $opt_canu_bin, "unicycler" => $opt_unicycler_bin);
	my $fail = 0;
	foreach my $p (sort keys %programs){
		my $path = $programs{$p};
		my $result = 'ok';
		if(! -X $path){
			$result = '!fail!';
			$fail = 1;
		}
		printf "%-10s%6s\t%s\n", $p, $result, $path;
	}
	exit($fail);
}


#=============================== Get bin PATH ======================================================
my $pgcgap_dir;
my $bin = `which pgcgap`;
if ($bin=~/(.+)\/pgcgap/) {
	$pgcgap_dir = $1;
}


#===================================================================================================
my $working_dir = getcwd;
#system("mkdir -p Results");

my $cpu_count = `cat /proc/cpuinfo| grep "cpu cores"| uniq`;
$cpu_count =~ /.+?(\d+)/;
my $threads_half = $1;

#my $threads_half = CPU();

# Genome Assemble with"Abyss" or "Canu"
if ($opt_Assemble) {
	system("mkdir -p Results/Assembles/Scaf");
	system("mkdir -p Results/Assembles/FASTQ_Preprocessor");#2020/4/15
	my $unicycler_set;
	if ($bin=~/(.+)bin\/pgcgap/) {
		$unicycler_set = $1 . "lib/python3.6/site-packages/unicycler/settings.py";
		if (-e $unicycler_set) {
			system("sed -i 's/RACON_POLISH_LOOP_COUNT_HYBRID = .*/RACON_POLISH_LOOP_COUNT_HYBRID = 2/g' $unicycler_set");
			system("sed -i 's/RACON_POLISH_LOOP_COUNT_LONG_ONLY = .*/RACON_POLISH_LOOP_COUNT_LONG_ONLY = 4/g' $unicycler_set");
		}else {
			print "Can not find the unicycler setting file\n";
		}
	}
	if ($opt_platform eq "illumina" and $opt_assembler eq "abyss") {
		#print "Performing --Assemble function for Illunina data with abyss...\n\n";
		system("mkdir -p Results/Assembles/Illumina");
		system("mkdir -p Results/Assembles/Scaf/Illumina");
		if ($opt_threads > $threads_half) {
			$opt_threads = $threads_half;
		}
		chdir $opt_ReadsPath;
		my @files = glob("*$opt_reads1");
		my %lists;
		foreach (@files) {
			if (/(\S+)$opt_reads1/) {
				$lists{$1} = "1";
			}
		}

		my @lists = keys %lists;

		foreach my $name(@lists) {
			my $read1 = $name . $opt_reads1;
			my $read2 = $name . $opt_reads2;
			my $str = substr($read1,0,(length($read1)-$opt_suffix_len));
			my $fastp_out1 = $name . ".fastp" . $opt_reads1;#2020/4/15
			my $fastp_out2 = $name . ".fastp" . $opt_reads2;#2020/4/15
			my $fastph = $str . ".fastp.html";#2020/4/15
			my $fastpj = $str . ".fastp.json";#2020/4/15
			print "Performing reads preprocessor with fastp\n\n";#2020/4/15
			system("fastp -i $read1 -I $read2 -o $fastp_out1 -O $fastp_out2 -j $fastpj -h $fastph -w $opt_threads -3");#2020/4/15
			print "Performing --Assemble function for Illunina data with abyss...\n\n";#2020/4/15
			system("abyss-pe name=$str k=$opt_kmmer in='$fastp_out1 $fastp_out2' np=$opt_threads");#2020/4/15
			#system("abyss-pe name=$str k=$opt_kmmer in='$read1 $read2' np=$opt_threads");
			print "Assemble complete !\n";
			my $assem = $str . "_assembly";
			system("mkdir -p $working_dir/Results/Assembles/Illumina/$assem");
			my $scaf = $str . "-8.fa";
#			system("mkdir Over");
			system("cp $scaf $working_dir/Results/Assembles/Scaf/Illumina/");
#			system("mv $read1 $read2 Over/");
			system("mv $str*.dot* $str*.fa $str*.path* $str*.dist $str*.fai $str*stats* $str*.hist coverage.hist $working_dir/Results/Assembles/Illumina/$assem/");
			system("mv $fastp_out1 $fastp_out2 $fastph $fastpj $working_dir/Results/Assembles/FASTQ_Preprocessor");#2020/4/15
		}
		chdir $working_dir;
		system("realpath $working_dir/Results/Assembles/Scaf/Illumina/* >> scaf.list");
		chdir "$working_dir/Results/Assembles/Scaf/Illumina/";
		my @fas = glob("*-8.fa");
		foreach  (@fas) {
			$_=~/(\S+)-8.fa/;
			my $outpre = $1 . ".prefilter.stats";
			my $outfilter = $1 . ".filtered.stats";
			my $filterscaf = $1. ".filtered.fas";
			getstats($_, $outpre);
			#system("perl N50Stat.pl -i $_ -o $outpre");
			lenfilter($_, $filterscaf, $opt_filter_length);
			getstats($filterscaf, $outfilter);
			#system("perl N50Stat.pl -i $filterscaf -o $outfilter");
		}
#		my $time_assemble = time();
#		my $time_assemblex = ($time_assemble - $time_start)/3600;
#		print "The 'Assemble' program runs for $time_assemblex hours.\n\n";
		chdir $working_dir;
	}elsif ($opt_platform eq "illumina" and $opt_assembler eq "spades") {
		print "Performing --Assemble function for Illunina data with spades...\n\n";
		system("mkdir -p Results/Assembles/Illumina");
		system("mkdir -p Results/Assembles/Scaf/Illumina");
		
		chdir $opt_ReadsPath;
		my @files = glob("*$opt_reads1");
		my %lists;
		foreach (@files) {
			if (/(\S+)$opt_reads1/) {
				$lists{$1} = "1";
			}
		}

		my @lists = keys %lists;

		foreach my $name(@lists) {
			my $read1 = $name . $opt_reads1;
			my $read2 = $name . $opt_reads2;
			my $str = substr($read1,0,(length($read1)-$opt_suffix_len));
			my $fastp_out1 = $name . ".fastp" . $opt_reads1;#2020/4/15
			my $fastp_out2 = $name . ".fastp" . $opt_reads2;#2020/4/15
			my $fastph = $str . ".fastp.html";#2020/4/15
			my $fastpj = $str . ".fastp.json";#2020/4/15
			print "Performing reads preprocessor with fastp\n\n";#2020/4/15
			system("fastp -i $read1 -I $read2 -o $fastp_out1 -O $fastp_out2 -j $fastpj -h $fastph -w $opt_threads -3");#2020/4/15
			print "Performing --Assemble function for Illunina data with unicycler...\n\n";#2020/4/15
			system("unicycler -1 $fastp_out1 -2 $fastp_out2 -t $opt_threads -o $str");#2020/4/15
			#print "Assembling...\n";
			#system("unicycler -1 $read1 -2 $read2 -t $opt_threads -o $str");
			print "Assemble complete !\n";
			my $scaf = $str . "-8.fa";
			system("mv $str/assembly.fasta $str/$scaf");
			system("cp $str/$scaf $working_dir/Results/Assembles/Scaf/Illumina/");
			system("cp -rf $str $working_dir/Results/Assembles/Illumina/");
			system("mv $fastp_out1 $fastp_out2 $fastph $fastpj $working_dir/Results/Assembles/FASTQ_Preprocessor");#2020/4/15
		}
		chdir $working_dir;
		system("realpath $working_dir/Results/Assembles/Scaf/Illumina/* >> scaf.list");
		chdir "$working_dir/Results/Assembles/Scaf/Illumina/";
		my @fas = glob("*-8.fa");
		foreach  (@fas) {
			$_=~/(\S+)-8.fa/;
			my $outpre = $1 . ".prefilter.stats";
			my $outfilter = $1 . ".filtered.stats";
			my $filterscaf = $1. ".filtered.fas";
			getstats($_, $outpre);
			#system("perl N50Stat.pl -i $_ -o $outpre");
			lenfilter($_, $filterscaf, $opt_filter_length);
			getstats($filterscaf, $outfilter);
			#system("perl N50Stat.pl -i $filterscaf -o $outfilter");
		}
#		my $time_assemble = time();
#		my $time_assemblex = ($time_assemble - $time_start)/3600;
#		print "The 'Assemble' program runs for $time_assemblex hours.\n\n";
		chdir $working_dir;
	}elsif ($opt_platform eq "illumina" and $opt_assembler eq "auto") {
		#print "Performing --Assemble function for Illunina data with abyss...\n\n";
		system("mkdir -p Results/Assembles/Illumina");
		system("mkdir -p Results/Assembles/Scaf/Illumina");
		
		chdir $opt_ReadsPath;
		my @files = glob("*$opt_reads1");
		my %lists;
		foreach (@files) {
			if (/(\S+)$opt_reads1/) {
				$lists{$1} = "1";
			}
		}

		my @lists = keys %lists;

		foreach my $name(@lists) {
			my $read1 = $name . $opt_reads1;
			my $read2 = $name . $opt_reads2;
			my $str = substr($read1,0,(length($read1)-$opt_suffix_len));
			my $fastp_out1 = $name . ".fastp" . $opt_reads1;#2020/4/15
			my $fastp_out2 = $name . ".fastp" . $opt_reads2;#2020/4/15
			my $fastph = $str . ".fastp.html";#2020/4/15
			my $fastpj = $str . ".fastp.json";#2020/4/15
			print "Performing reads preprocessor with fastp\n\n";#2020/4/15
			system("fastp -i $read1 -I $read2 -o $fastp_out1 -O $fastp_out2 -j $fastpj -h $fastph -w $opt_threads -3");#2020/4/15
			print "Performing --Assemble function for Illunina data with abyss...\n\n";#2020/4/15
			system("abyss-pe name=$str k=$opt_kmmer in='$fastp_out1 $fastp_out2' np=$threads_half");#2020/4/15
			#print "Assembling...\n";
			#system("abyss-pe name=$str k=$opt_kmmer in='$read1 $read2' np=$threads_half");
			print "Assemble complete !\n";
			my $assem = $str . "_assembly";
			system("mkdir -p $working_dir/Results/Assembles/Illumina/$assem");
			my $scaf = $str . "-8.fa";
			#system("cp $scaf $working_dir/Results/Assembles/Scaf/Illumina/");
			#system("mv $str*.dot* $str*.fa $str*.path* $str*.dist $str*.fai $str*stats* $str*.hist coverage.hist $working_dir/Results/Assembles/Illumina/$assem/");
			my $stats = $str . "-stats.tab";
			print "Checking the assembly stats...\n";
			open IN, "$stats" || die;
			my @array = <IN>;
			my $lastline = $array[-1];#get the last line of the file
			my @stats = split "\t", $lastline;
			if ($stats[5] < 50000) {
				system("mv $str*.dot* $str*.fa $str*.path* $str*.dist $str*.fai $str*stats* $str*.hist coverage.hist $working_dir/Results/Assembles/Illumina/$assem/");
				print "The N50 is less than 50k, now performing --Assemble function for Illunina data with unicycler to try to get a better assembly result...\n\n";
				system("unicycler -1 $fastp_out1 -2 $fastp_out2 -t $opt_threads -o $str");#2020/4/15
				#system("unicycler -1 $read1 -2 $read2 -t $opt_threads -o $str");
				print "Assemble completed!\n";
				#my $scaf = $str . "-8.fa";
				system("mv $str/assembly.fasta $str/$scaf");
				system("cp $str/$scaf $working_dir/Results/Assembles/Scaf/Illumina/");
				system("mv -f $str $working_dir/Results/Assembles/Illumina/");
			}else {
				system("cp $scaf $working_dir/Results/Assembles/Scaf/Illumina/");
				system("mv $str*.dot* $str*.fa $str*.path* $str*.dist $str*.fai $str*stats* $str*.hist coverage.hist $working_dir/Results/Assembles/Illumina/$assem/");
			}
			system("mv $fastp_out1 $fastp_out2 $fastph $fastpj $working_dir/Results/Assembles/FASTQ_Preprocessor");#2020/4/15
		}
		chdir $working_dir;
		system("realpath $working_dir/Results/Assembles/Scaf/Illumina/* >> scaf.list");
		chdir "$working_dir/Results/Assembles/Scaf/Illumina/";
		my @fas = glob("*-8.fa");
		foreach  (@fas) {
			$_=~/(\S+)-8.fa/;
			my $outpre = $1 . ".prefilter.stats";
			my $outfilter = $1 . ".filtered.stats";
			my $filterscaf = $1. ".filtered.fas";
			getstats($_, $outpre);
			#system("perl N50Stat.pl -i $_ -o $outpre");
			lenfilter($_, $filterscaf, $opt_filter_length);
			getstats($filterscaf, $outfilter);
			#system("perl N50Stat.pl -i $filterscaf -o $outfilter");
		}
#		my $time_assemble = time();
#		my $time_assemblex = ($time_assemble - $time_start)/3600;
#		print "The 'Assemble' program runs for $time_assemblex hours.\n\n";
		chdir $working_dir;
	}elsif ($opt_platform eq "pacbio") {
		#print "Performing --Assemble function for PacBio data...\n\n";
		system("mkdir -p Results/Assembles/PacBio");
		system("mkdir -p Results/Assembles/Scaf/PacBio");
		chdir $opt_ReadsPath;
		my $dir_ReadsPath = getcwd;
		my @files = glob("*$opt_reads1");
		foreach (@files) {
			my $name = substr($_,0,(length($_)-$opt_suffix_len));
#			if (/(\S+)$opt_reads1/) {
#				my $name = $1;
			my $outdir = "$working_dir/Results/Assembles/PacBio/" . $name;
			my $cir_outdir = $outdir . "Circlator";
			my $scaf = $name . ".contigs.fasta";
			my $correct_reads = $name . ".correctedReads.fasta.gz";
			my $cir_scaf = $name . ".fixstart.fasta";
			my $sam = $name . ".sam";
			my $polished_scaf = $name . ".polished.contigs.fasta";
			##my $fastp_out1 = $name . ".fastp" . $opt_reads1;#2020/4/15
			##my $fastph = $name . ".fastp.html";#2020/4/15
			##my $fastpj = $name . ".fastp.json";#2020/4/15
			##print "Performing reads preprocessor with fastp\n\n";#2020/4/15
			##system("fastp -i $_ -o $fastp_out1 -j $fastpj -h $fastph -w $opt_threads");#2020/4/15
			print "Performing --Assemble function for PacBio data...\n\n";
			##system("canu -p $name -d $outdir genomeSize=$opt_genomeSize maxThreads=$opt_threads useGrid=false -pacbio-raw $fastp_out1");#2020/4/15
			system("canu -p $name -d $outdir genomeSize=$opt_genomeSize maxThreads=$opt_threads useGrid=false -pacbio-raw $_");
			print "Begin polishing with racon\n\n";
			chdir $outdir;
			system("bwa index $scaf");
			system("bwa mem -x pacbio -t $opt_threads -o $sam $scaf $dir_ReadsPath/$_");
			system("racon -t $opt_threads $dir_ReadsPath/$_ $sam $scaf > $polished_scaf");
			##system("circlator all --assembler canu $outdir/$scaf $outdir/$correct_reads $cir_outdir");
			system("cp $polished_scaf $working_dir/Results/Assembles/Scaf/PacBio/");
			##system("mv $fastp_out1 $fastph $fastpj $working_dir/Results/Assembles/FASTQ_Preprocessor");#2020/4/15
			##system("cp $cir_outdir/06.fixstart.fasta $working_dir/Results/Assembles/Scaf/PacBio/$cir_scaf");
#			}
		}
		chdir $working_dir;
		system("realpath Results/Assembles/Scaf/PacBio/* >> scaf.list");
		chdir "$working_dir/Results/Assembles/Scaf/PacBio/";
		my @fas = glob("*.fasta");
		foreach  (@fas) {
			$_=~/(\S+).fasta/;
			my $outpre = $1 . ".prefilter.stats";
			my $outfilter = $1 . ".filtered.stats";
			my $filterscaf = $1. ".filtered.fas";
			getstats($_, $outpre);
			#system("perl N50Stat.pl -i $_ -o $outpre");
			lenfilter($_, $filterscaf, $opt_filter_length);
			getstats($filterscaf, $outfilter);
			#system("perl N50Stat.pl -i $filterscaf -o $outfilter");
		}
		chdir $working_dir;
	}elsif ($opt_platform eq "oxford") {
		#print "Performing --Assemble function for Oxford Nanopore data...\n\n";
		system("mkdir -p Results/Assembles/Oxford");
		system("mkdir -p Results/Assembles/Scaf/Oxford");
		chdir $opt_ReadsPath;
		my $dir_ReadsPath = getcwd;
		my @files = glob("*$opt_reads1");
		foreach (@files) {
			my $name = substr($_,0,(length($_)-$opt_suffix_len));
#			if (/(\S+)$opt_reads1/) {
#				my $name = $1;
			my $outdir = "$working_dir/Results/Assembles/Oxford/" . $name;
			my $cir_outdir = $outdir . "Circlator";
			my $scaf = $name . ".contigs.fasta";
			my $correct_reads = $name . ".correctedReads.fasta.gz";
			my $cir_scaf = $name . ".fixstart.fasta";
			my $sam = $name . ".sam";
			my $polished_scaf = $name . ".polished.contigs.fasta";
			##my $fastp_out1 = $name . ".fastp" . $opt_reads1;#2020/4/15
			##my $fastph = $name . ".fastp.html";#2020/4/15
			##my $fastpj = $name . ".fastp.json";#2020/4/15
			##print "Performing reads preprocessor with fastp\n\n";#2020/4/15
			##system("fastp -i $_ -o $fastp_out1 -j $fastpj -h $fastph -w $opt_threads");#2020/4/15
			print "Performing --Assemble function for Oxford Nanopore data...\n\n";
			##system("canu -p $name -d $outdir genomeSize=$opt_genomeSize maxThreads=$opt_threads useGrid=false -nanopore-raw $fastp_out1");#2020/4/15
			system("canu -p $name -d $outdir genomeSize=$opt_genomeSize maxThreads=$opt_threads useGrid=false -nanopore $_");
			print "Begin polishing with racon\n\n";
			chdir $outdir;
			system("bwa index $scaf");
			system("bwa mem -x ont2d -t $opt_threads -o $sam $scaf $dir_ReadsPath/$_");
			system("racon -t $opt_threads $dir_ReadsPath/$_ $sam $scaf > $polished_scaf");
			##system("circlator all --assembler canu --merge_min_id 85 --merge_breaklen 1000 $outdir/$scaf $outdir/$correct_reads $cir_outdir");
			system("cp $polished_scaf $working_dir/Results/Assembles/Scaf/Oxford/");
			##system("mv $fastp_out1 $fastph $fastpj $working_dir/Results/Assembles/FASTQ_Preprocessor");#2020/4/15
			##system("cp $cir_outdir/06.fixstart.fasta $working_dir/Results/Assembles/Scaf/Oxford/$cir_scaf");
#			}
		}
		chdir $working_dir;
		system("realpath Results/Assembles/Scaf/Oxford/* >> scaf.list");
		chdir "$working_dir/Results/Assembles/Scaf/Oxford/";
		my @fas = glob("*.fasta");
		foreach  (@fas) {
			$_=~/(\S+).fasta/;
			my $outpre = $1 . ".prefilter.stats";
			my $outfilter = $1 . ".filtered.stats";
			my $filterscaf = $1. ".filtered.fas";
			getstats($_, $outpre);
			#system("perl N50Stat.pl -i $_ -o $outpre");
			lenfilter($_, $filterscaf, $opt_filter_length);
			getstats($filterscaf, $outfilter);
			#system("perl N50Stat.pl -i $filterscaf -o $outfilter");
		}
		chdir $working_dir;
	}elsif ($opt_platform eq "hybrid") {
		print "Performing --Assemble function for hybrid data...\n\n";
		chdir $opt_ReadsPath;
		system("unicycler -1 $opt_short1 -2 $opt_short2 -l $opt_long -o $opt_hout -t $opt_threads");
		chdir $working_dir;
	}
}


if ($opt_ACC) {
	if ($opt_Assess) {
		chdir $opt_scafPath;
		system("genome_LenFilter_stats.pl --Scaf_suffix $opt_Scaf_suffix --filter_length $opt_filter_length");
		system("get_stats_summary.pl --Scaf_suffix $opt_Scaf_suffix");
		chdir $working_dir;
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

{
	my $As = 0;
	my $Ts = 0;
	my $Gs = 0;
	my $Cs = 0;
	my $Ns = 0;
	sub getstats {
		# Parameter variables
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

sub CPU{
	my %options;
	my $info = Sys::Info->new;
	my $cpu  = $info->device( CPU => %options );
	my $threads_num = $cpu->ht || 2;
	my $threads_half = $threads_num/2;
	return $threads_half;
}
sub printAssemble{
	print "[--platform (STRING)] Sequencing Platform, 'illumina', 'pacbio', 'oxford' and 'hybrid' available ( Default illumina )\n";
	print "[--assembler (STRING)] Software used for illumina reads assembly, 'abyss', 'spades' and 'auto' available ( Default auto )\n";
	print "[--ReadsPath (PATH)] Reads of all strains as file paths ( Default ./Reads/Illumina )\n";
	print "[--reads1 (STRING)] The suffix name of reads 1 ( for example: if the name of reads 1 is 'YBT-1520_L1_I050.R1.clean.fastq.gz', 'YBT-1520' is the strain same, the suffix name should be '.R1.clean.fastq.gz' )\n";
	print "[--reads2 (STRING)] The suffix name of reads 2( for example: if the name of reads 2 is 'YBT-1520_2.fq', the suffix name should be '_2.fq' )\n";
	print "[--suffix_len (INT)] (Strongly recommended) The suffix length of the reads file, that is the length of the reads name minus the length of the strain name. For example the --suffix_len of 'YBT-1520_L1_I050.R1.clean.fastq.gz' is 26 ( 'YBT-1520' is the strain name ) ( Default 0 )\n";
	print "[--kmmer (INT)] k-mer size for genome assembly of Illumina data with abyss( Default 81 )\n";
	print "[--filter_length (INT)] Sequences shorter than the 'filter_length' will be deleted from the assembled genomes. ( Default 200 )";
	print "[--genomeSize (STRING)] An estimate of the size of the genome. Common suffixes are allowed, for example, 3.7m or 2.8g. Needed by PacBio data and Oxford data ( Default Unset )\n";
	print "[--short1 (STRING)] FASTQ file of first short reads in each pair. Needed by hybrid assembly ( Default Unset )\n";
	print "[--short2 (STRING)] FASTQ file of second short reads in each pair. Needed by hybrid assembly ( Default Unset )\n";
	print "[--long (STRING)] FASTQ or FASTA file of long reads. Needed by hybrid assembly ( Default Unset )\n";
	print "[--hout (STRING)] Output directory for hybrid assembly ( Default ../../Results/Assembles/Hybrid )\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
}

sub printACC{
	print "Applets in ACC include 'Assess' now\n";
	print "Parameters for Assess include the following:\n    [--scafPath (PATH)] Path for contigs/scaffolds ( Default 'Results/Assembles/Scaf/Illumina' )\n    [--Scaf_suffix (STRING)] The suffix of scaffolds or genome files ( Default -8.fa )\n    [--filter_length (INT)] Sequences shorter than the 'filter_length' will be deleted from the assembled genomes. ( Default 200 )\n\n";
}

sub printExamples{
	print "\nThe main usage is as follows, visit the official website for step by step examples: https://liaochenlanruo.github.io/pgcgap/\n\n";

	print "Example 1: Conduct pair-end reads assembly\n\n";

	print "         pgcgap --Assemble --platform illumina --assembler abyss --ReadsPath <PATH> --reads1 <reads1 suffix> --reads2 <reads2 suffix> --suffix_len <INT> --kmmer <INT> --threads <INT>\n";
	print "         pgcgap --Assemble --platform illumina --assembler spades --ReadsPath <PATH> --reads1 <reads1 suffix> --reads2 <reads2 suffix> --suffix_len <INT> --threads <INT>\n";
	print "         pgcgap --Assemble --platform illumina --assembler auto --ReadsPath <PATH> --reads1 <reads1 suffix> --reads2 <reads2 suffix> --suffix_len <INT> --kmmer <INT> --threads <INT>\n\n";

	print "Example 2: Conduct PacBio/Oxford reads assembly\n\n";

	print "         pgcgap --Assemble --platform [pacbio|oxford] --ReadsPath <PATH> --reads1 <reads suffix> --suffix_len <INT> --genomeSize <STRING> --threads <INT>\n\n";

	print "Example 3: Conduct hybrid assembly\n\n";

	print "         pgcgap --Assemble --platform hybrid --ReadsPath <PATH> --short1 <pair-end-reads1> --short2 <pair-end-reads2> --long <long-reads> --hout <output_dir> --threads <INT>\n\n";


	print "Example 4: Perform the short sequences filter from the assembled genome and get the genome status\n\n";

	print "         pgcgap --ACC --Assess --scafPath <PATH> --Scaf_suffix <STRING> --filter_length <INT>\n\n";
}

if ( grep {$_ eq "Assemble"} @ARGV ){
	printAssemble();
	exit 0;
}

if ( grep {$_ eq "ACC"} @ARGV ) {
	printACC();
	exit 0;
}

if ( grep {$_ eq "Examples"} @ARGV ){
	printExamples();
	exit 0;
}