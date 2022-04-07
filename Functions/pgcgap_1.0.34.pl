#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper;
use File::Tee qw(tee);
use Cwd;
use List::Util qw(sum min max);
use File::Basename;
use POSIX;#ceil function
#use Sys::Info;
#use Sys::Info::Constants qw( :device_cpu );
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use Term::ANSIColor qw(:constants);

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

  Show parameters for each module: pgcgap [Assemble|Annotate|ANI|AntiRes|CoreTree|MASH|OrthoF|Pan|pCOG|VAR|STREE|ACC]

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

=over 30

=item B<[--check-update]>

Check if there is a new version of PGCGAP that can be upgraded

=back

=cut

$options{'check-update'} = \( my $opt_check_update = 0 );

=over 30

=item B<[--setup-COGdb]>

Setup COG database. Users should execute "pgcgap --setup-COGdb" after the first installation of pgcgap

=back

=cut

$options{'setup-COGdb'} = \( my $opt_setup_COGdb );

=over 30

=item B<[--setup-COGdb2]>

Alternate method to setup COG database. This option can be used to download and setup the COG database when network access is not available with 'setup-COGdb'

=back

=cut

$options{'setup-COGdb2'} = \( my $opt_setup_COGdb2 );

=head2 *********************************************** Modules ************************************************

=for text



=over 30

=item B<[--All]>

Perform Assemble, Annotate, CoreTree, Pan, OrthoF, ANI, MASH, AntiRes and pCOG functions with one command

=back

=cut

$options{'All'} = \(my $opt_All);

=over 30

=item B<[--Assemble]>

Assemble reads (short, long or hybrid) into contigs

=back

=cut

$options{'Assemble'} = \(my $opt_Assemble);

=over 30

=item B<[--Annotate]>

Genome annotation

=back

=cut

$options{'Annotate'} = \(my $opt_Annotate);

=over 30

=item B<[--CoreTree]>

Construct single-copy core proteins tree and core SNPs tree

=back

=cut

$options{'CoreTree'} = \(my $opt_CoreTree);

=over 30

=item B<[--Pan]>

Run "roary" pan genome pipeline with gff3 files, and construct a phylogenetic tree with the sing-copy core proteins called by roary

=back

=cut

$options{'Pan'} = \(my $opt_Pan);

=over 30

=item B<[--OrthoF]>

Identify orthologous protein sequence families

=back

=cut

$options{'OrthoF'} = \(my $opt_OrthoF);

=over 30

=item B<[--ANI]>

Compute whole-genome Average Nucleotide Identity ( ANI )

=back

=cut

$options{'ANI'} = \(my $opt_ANI);

=over 30

=item B<[--MASH]>

Genome and metagenome similarity estimation using MinHash

=back

=cut

$options{'MASH'} = \(my $opt_MASH);

=over 30

=item B<[--pCOG]>

Run COG annotation for each strain (*.faa), and generate a table containing the relative abundance of each flag for all strains

=back

=cut

$options{'pCOG'} = \(my $opt_pCOG);

=over 30

=item B<[--VAR]>

Rapid haploid variant calling and core genome alignment

=back

=cut

$options{'VAR'} = \(my $opt_VAR);

=over 30

=item B<[--AntiRes]>

Screening for antimicrobial and virulence genes

=back

=cut

$options{'AntiRes'} = \(my $opt_AntiRes);

=over 30

=item B<[--STREE]>

Construct a phylogenetic tree based on multiple sequences in one file

=back

=cut

$options{'STREE'} = \(my $opt_STREE);

=over 30

=item B<[--ACC]>

Other useful gadgets

=back

=cut

$options{'ACC'} = \(my $opt_ACC);

=head2 *********************************************** Global Options *****************************************

=for text



=over 30

=item B<[--strain_num (INT)]>

I<[Required by "All", "CoreTree", "Pan", "VAR" and "pCOG"]> The total number of strains used for analysis, not including the reference genome

=back

=cut

$options{'strain_num=i'} = \( my $opt_strain_num );

=over 30

=item B<[--ReadsPath (PATH)]>

I<[Required by "All", "Assemble" and "VAR"]> Reads of all strains as file paths ( Default ./Reads/Illumina )

=back

=cut

$options{'ReadsPath=s'} = \( my $opt_ReadsPath = "./Reads/Illumina" );

=over 30

=item B<[--scafPath (PATH)]>

I<[Required by "All", "Assess", "Annotate", "MASH" and "AntiRes"]> Path for contigs/scaffolds ( Default "Results/Assembles/Scaf/Illumina" )

=back

=cut

$options{'scafPath=s'} = \(my $opt_scafPath = "./Results/Assembles/Scaf/Illumina");

=over 30

=item B<[--AAsPath (PATH)]>

I<[Required by "All", "CoreTree", "OrthoF" and "pCOG"]> Amino acids of all strains as fasta file paths, ( Default "./Results/Annotations/AAs" )

=back

=cut

$options{'AAsPath=s'} = \( my $opt_AAsPath = "./Results/Annotations/AAs" );

=over 30

=item B<[--reads1 (STRING)]>

I<[Required by "All", "Assemble" and "VAR"]> The suffix name of reads 1 ( for example: if the name of reads 1 is "YBT-1520_L1_I050.R1.clean.fastq.gz", "YBT-1520" is the strain same, so the suffix name should be ".R1.clean.fastq.gz" )

=back

=cut

$options{'reads1=s'} = \(my $opt_reads1);

=over 30

=item B<[--reads2 (STRING)]>

I<[Required by "All", "Assemble" and "VAR"]> The suffix name of reads 2( for example: if the name of reads 2 is "YBT-1520_2.fq", the suffix name should be _2.fq" )

=back

=cut

$options{'reads2=s'} = \(my $opt_reads2);

=over 30

=item B<[--Scaf_suffix (STRING)]>

The suffix of scaffolds or genome files [Required by "All", "Assess", "Annotate", "MASH", "ANI" and "AntiRes"]. This is an important parameter that must be set ( Default .filtered.fas )

=back

=cut

$options{'Scaf_suffix=s'} = \( my $opt_Scaf_suffix = ".filtered.fas" );

=over 30

=item B<[--filter_length (INT)]>

I<[Required]> Sequences shorter than the 'filter_length' will be deleted from the assembled genomes [Required by "All", "Assemble" and "Assess"]. ( Default 200 )

=back

=cut

$options{'filter_length=i'} = \(my $opt_filter_length = 200);

=over 30

=item B<[--codon (INT)]>

I<[Required by "All", "Annotate", "CoreTree" and "Pan"]> Translation table ( Default 11 )

=back

=cut

$options{'codon=i'} = \( my $opt_codon = 11 );

=begin text

                                  1   Universal code
                                  2   Vertebrate mitochondrial code
                                  3   Yeast mitochondrial code
                                  4   Mold, Protozoan, and Coelenterate Mitochondrial code and Mycoplasma/Spiroplasma code
                                  5   Invertebrate mitochondrial
                                  6   Ciliate, Dasycladacean and Hexamita nuclear code
                                  9   Echinoderm and Flatworm mitochondrial code
                                  10  Euplotid nuclear code
                                  11  Bacterial, archaeal and plant plastid code ( Default )
                                  12  Alternative yeast nuclear code
                                  13  Ascidian mitochondrial code
                                  14  Alternative flatworm mitochondrial code
                                  15  Blepharisma nuclear code
                                  16  Chlorophycean mitochondrial code
                                  21  Trematode mitochondrial code
                                  22  Scenedesmus obliquus mitochondrial code
                                  23  Thraustochytrium mitochondrial code
    

=end text

=over 30

=item B<[--suffix_len (INT)]>

I<[Required by "All", "Assemble" and "VAR"]> B<(Strongly recommended)> The suffix length of the reads file, that is the length of the reads name minus the length of the strain name. For example the --suffix_len of "YBT-1520_L1_I050.R1.clean.fastq.gz" is 26 ( "YBT-1520" is the strain name ) ( Default 0 )

=back

=cut

$options{'suffix_len=i'} = \(my $opt_suffix_len = 0);

=over 30

=item B<[--logs (STRING)]>

Name of the log file ( Default Logs.txt )

=back

=cut

$options{'logs=s'} = \( my $opt_logs = "Logs.txt" );

=over 30

=item B<[--fasttree]>

I<[Can be used with "CoreTree", "Pan" and "OrthoF"]> Use FastTree to construct phylogenetic tree quickly instead of IQ-TREE.

=back

=cut

$options{'fasttree'} = \(my $opt_fasttree);

=over 30

=item B<[--bsnum (INT)]>

I<[Required by "CoreTree", "Pan", "OrthoF", "STREE", and "VAR"]> Replicates for bootstrap of IQ-TREE. ( Default 500 )

=back

=cut

$options{'bsnum=i'} = \( my $opt_bsnum = "500");

=over 30

=item B<[--fastboot (INT)]>

I<[Required by "CoreTree", "Pan", "OrthoF", "STREE", and "VAR"]> Replicates for ultrafast bootstrap of IQ-TREE. ( must >= 1000, Default 1000 )

=back

=cut

$options{'fastboot=i'} = \( my $opt_fastboot = "1000");

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

=head3 ========================== Options of "--Annotate" for genome annotation ==============================

=for text



=begin html

If you use the results of "--Annotate" function in your work, please also cite:

</br>

</br>Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063

=end html

=over 30

=item B<[--genus (STRING)]>

Genus name of the strain ( Default "NA" )

=back

=cut

$options{'genus=s'} = \(my $opt_genus = "NA");

=over 30

=item B<[--species (STRING)]>

Species name of the strain ( Default "NA" )

=back

=cut

$options{'species=s'} = \(my $opt_species = "NA");

=head3 ========================== Options for "--CoreTree" constructing ======================================

=for text



=begin html

If you use the results of "--CoreTree" function in your work, please also cite:

</br>

</br>CD-HIT: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik. Bioinformatics, (2006) 22:1658-1659

</br>

</br>CD-HIT: accelerated for clustering the next generation sequencing data", Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu & Weizhong Li. Bioinformatics, (2012) 28:3150-3152

</br>

</br>Katoh K, Misawa K, Kuma K, Miyata T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res. 2002;30(14):3059-3066

</br>

</br>Mikita Suyama, David Torrents, and Peer Bork (2006) PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Res. 34, W609-W612

</br>

</br>"SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments", Andrew J. Page, Ben Taylor, Aidan J. Delaney, Jorge Soares, Torsten Seemann, Jacqueline A. Keane, Simon R. Harris, Microbial Genomics 2(4), (2016)

=end html

=over 30

=item B<[--CDsPath (PATH)]>

I<[Required]> CDs of all strains as fasta file paths, ( Default "./Results/Annotations/CDs" )

=back

=cut

$options{'CDsPath=s'} = \( my $opt_CDsPath = "./Results/Annotations/CDs" );

=over 30

=item B<[-c (FLOAT)]>

Sequence identity threshold, ( Default 0.5)

=back

=cut

$options{'c=f'} = \( my $opt_c = 0.5 );

=over 30

=item B<[-n (INT)]>

Word_length, -n 2 for thresholds 0.4-0.5, -n 3 for thresholds 0.5-0.6, -n 4 for thresholds 0.6-0.7, -n 5 for thresholds 0.7-1.0 ( Default 2 )

=back

=cut

$options{'n=i'} = \( my $opt_n = 2 );

=over 30

=item B<[-G (INT)]>

Use global (set to 1) or local (set to 0) sequence identity, ( Default 0 )

=back

=cut

$options{'G=i'} = \( my $opt_G = 0 );

=over 30

=item B<[-t (INT)]>

Tolerance for redundance ( Default 0 )

=back

=cut

$options{'t=i'} = \( my $opt_t = 0 );

=over 30

=item B<[-aL (FLOAT)]>

Alignment coverage for the longer sequence. If set to 0.9, the alignment must covers 90% of the sequence ( Default 0.5 )

=back

=cut

$options{'aL=f'} = \( my $opt_aL = 0.5 );

=over 30

=item B<[-aS (FLOAT)]>

Alignment coverage for the shorter sequence. If set to 0.9, the alignment must covers 90% of the sequence ( Default 0.7 )

=back

=cut

$options{'aS=f'} = \( my $opt_aS = 0.7 );

=over 30

=item B<[-g (INT)]>

If set to 0, a sequence is clustered to the first cluster that meets the threshold (fast cluster). If set to 1, the program will cluster it into the most similar cluster that meets the threshold (accurate but slow mode, Default 1)

=back

=cut

$options{'g=i'} = \( my $opt_g = 1 );

=over 30

=item B<[-d (INT)]>

length of description in .clstr file. if set to 0, it takes the fasta defline and stops at first space ( Default 0 )

=back

=cut

$options{'d=i'} = \( my $opt_d = 0 );

=head3 ========================== Options for "--Pan" analysis ===============================================

=for text



=begin html

If you use the results of "--Pan" function in your work, please also cite:

</br>

</br>"Roary: Rapid large-scale prokaryote pan genome analysis", Andrew J. Page, Carla A. Cummins, Martin Hunt, Vanessa K. Wong, Sandra Reuter, Matthew T. G. Holden, Maria Fookes, Daniel Falush, Jacqueline A. Keane, Julian Parkhill, Bioinformatics, (2015). doi: http://dx.doi.org/10.1093/bioinformatics/btv421

=end html

=over 30

=item B<[--GffPath (PATH)]>

I<[Required]> Gff files of all strains as paths ( Default "./Results/Annotations/GFF" )

=back

=cut

$options{'GffPath=s'} = \( my $opt_GffPath = "./Results/Annotations/GFF" );

=over 30

=item B<[--PanTree]>

Construct a phylogenetic tree of single-copy core proteins called by roary

=back

=cut

$options{'PanTree'} = \(my $opt_PanTree);

=over 30

=item B<[--identi (INT)]>

Minimum percentage identity for blastp ( Default 95 )

=back

=cut

$options{'identi=i'} = \(my $opt_identi = "95");

=head3 ========================== Options for "--OrthoF" analysis ============================================

=for text



=begin html

If you use the results of "--OrthoF" function in your work, please also cite:

</br>

</br>Emms, D.M. and Kelly, S. (2018) OrthoFinder2: fast and accurate phylogenomic orthology analysis from gene sequences. bioRxiv

=end html

=over 30

=item B<[--Sprogram (STRING)]>

Sequence search program, Options: blast, mmseqs, blast_gz, diamond ( Default diamond )

=back

=cut

$options{'Sprogram=s'} = \( my $opt_Sprogram = "diamond" );

=head3 ========================== Options for "--ANI" analysis ===============================================

=for text



=begin html

If you use the results of "--ANI" function in your work, please also cite:

</br>

</br>Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat Commun. 2018;9(1):5114. Published 2018 Nov 30. doi:10.1038/s41467-018-07641-9

=end html

=over 30

=item B<[--queryL (FILE)]>

I<[Required]> The file containing full paths to query genomes, one per line ( Default scaf.list )

=back

=cut

$options{'queryL=s'} = \( my $opt_queryL = "scaf.list" );

=over 30

=item B<[--refL (FILE)]>

I<[Required]> The file containing full paths to reference genomes, one per line. ( Default scaf.list )

=back

=cut

$options{'refL=s'} = \( my $opt_refL = "scaf.list" );

=head3 ========================== Options for "--VAR" analysis ===============================================

=for text



=begin html

If you use the results of "--VAR" function in your work, please also cite:

</br>

</br>Joshi NA, Fass JN. (2011). Sickle: A sliding-window, adaptive, quality-based trimming tool for FastQ files (Version 1.33) [Software].  Available at https://github.com/najoshi/sickle.

</br>

</br>Seemann T (2015) snippy: fast bacterial variant calling from NGS reads https://github.com/tseemann/snippy.

</br>

</br>Croucher N. J., Page A. J., Connor T. R., Delaney A. J., Keane J. A., Bentley S. D., Parkhill J., Harris S.R. "Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins". Nucleic Acids Res. 2015 Feb 18;43(3):e15. doi: 10.1093/nar/gku1196

</br>

=end html

=over 30

=item B<[--refgbk (FILE)]>

I<[Required]> The B<full path and name> of reference genome in GENBANK format ( B<recommended> ), fasta format is also OK. For example: "/mnt/g/test/ref.gbk"

=back

=cut

$options{'refgbk=s'} = \( my $opt_refgbk );

=over 30

=item B<[--qualtype (STRING)]>

I<[Required]> Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)). ( Default sanger )

=back

=cut

$options{'qualtype=s'} = \(my $opt_qualtype = "sanger");

=over 30

=item B<[--qual (INT)]>

Threshold for trimming based on average quality in a window. ( Default 20 )

=back

=cut

$options{'qual=i'} = \(my $opt_qual = "20");

=over 30

=item B<[--length (INT)]>

Threshold to keep a read based on length after trimming. ( Default 20 )

=back

=cut

$options{'length=i'} = \(my $opt_length = "20");

=over 30

=item B<[--mincov (INT)]>

The minimum number of reads covering a site to be considered ( Default 10 )

=back

=cut

$options{'mincov=i'} = \(my $opt_mincov = "10");

=over 30

=item B<[--minfrac (FLOAT)]>

The minimum proportion of those reads which must differ from the reference ( Default 0.9 )

=back

=cut

$options{'minfrac=f'} = \(my $opt_minfrac = "0.9");

=over 30

=item B<[--minqual (INT)]>

The minimum VCF variant call "quality" ( Default 100 )

=back

=cut

$options{'minqual=i'} = \(my $opt_minqual = "100");

=over 30

=item B<[--ram (INT)]>

Try and keep RAM under this many GB ( Default 8 )

=back

=cut

$options{'ram=i'} = \(my $opt_ram = "8");

=head3 ========================== Options for "--AntiRes" analysis ===========================================

=over 30

=item B<[--db (STRING)]>

I<[Required]> The database to use, options: all, argannot, card, ecoh, ecoli_vf, megares, ncbi, plasmidfinder, resfinder and vfdb. ( Default all )

=back

=cut

$options{'db=s'} = \( my $opt_db = "all");

=over 30

=item B<[--identity (INT)]>

I<[Required]> Minimum %identity to keep the result, should be a number between 1 to 100. ( Default 75 )

=back

=cut

$options{'identity=i'} = \(my $opt_identity = "75");

=over 30

=item B<[--coverage (INT)]>

I<[Required]> Minimum %coverage to keep the result, should be a number between 0 to 100. ( Default 50 )

=back

=cut

$options{'coverage=i'} = \(my $opt_coverage = "50");

=head3 ========================== Options for "--STREE" ======================================================

=over 30

=item B<[--seqfile (STRING)]>

Path of the sequence file for analysis.

=back

=cut

$options{'seqfile=s'} = \( my $opt_seqfile);

=over 30

=item B<[--seqtype (STRING)]>

Type Of Sequence (p, d, c for Protein, DNA, Codons, respectively). ( Default p )

=back

=cut

$options{'seqtype=s'} = \( my $opt_seqtype = "p");

=head3 ========================== Options for "--pCOG" ======================================================

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

=head3 ========================== Options for "--ACC" ========================================================

=over 30

=item B<[--Assess (STRING)]>

Filter short sequences in the genome and assess the status of the genome.

=back

=cut

$options{'Assess'} = \( my $opt_Assess);

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

=item B<[--prodigal-bin (PATH)]>

Path to prodigal binary file. Default tries if prodigal is in PATH;

=back

=cut

$options{'prodigal-bin=s'} = \( my $opt_prodigal_bin = `which prodigal 2>/dev/null` );

=over 30

=item B<[--prokka-bin (PATH)]>

Path to prokka binary file. Default tries if prokka is in PATH;

=back

=cut

$options{'prokka-bin=s'} = \( my $opt_prokka_bin = `which prokka 2>/dev/null` );

=over 30

=item B<[--cd-hit-bin (PATH)]>

Path to cd-hit binary file. Default tries if cd-hit is in PATH;

=back

=cut

$options{'cd-hit-bin=s'} = \( my $opt_cdhit_bin = `which cd-hit 2>/dev/null` );

=over 30

=item B<[--mafft-bin (PATH)]>

Path to mafft binary file. Default tries if mafft is in PATH;

=back

=cut

$options{'mafft-bin=s'} = \( my $opt_mafft_bin = `which mafft 2>/dev/null` );

=over 30

=item B<[--pal2nal-bin (PATH)]>

Path to the pal2nal.pl binary file. Default tries if pal2nal.pl is in PATH;

=back

=cut

$options{'pal2nal-bin=s'} = \( my $opt_pal2nal_bin = `which pal2nal.pl 2>/dev/null` );

=over 30

=item B<[--snp-sites-bin (PATH)]>

Path to the snp-sites binary file. Default tries if snp-sites is in PATH;

=back

=cut

$options{'snp-sites-bin=s'} = \( my $opt_snpsites_bin = `which snp-sites 2>/dev/null` );

=over 30

=item B<[--roary-bin (PATH)]>

Path to the roary binary file. Default tries if roary is in PATH;

=back

=cut

$options{'roary-bin=s'} = \( my $opt_roary_bin = `which roary 2>/dev/null` );

=over 30

=item B<[--orthofinder-bin (PATH)]>

Path to the orthofinder binary file. Default tries if orthofinder is in PATH;

=back

=cut

$options{'orthofinder-bin=s'} = \( my $opt_orthofinder_bin = `which orthofinder 2>/dev/null` );

=over 30

=item B<[--fastANI-bin (PATH)]>

Path to the fastANI binary file. Default tries if fastANI is in PATH;

=back

=cut

$options{'fastANI-bin=s'} = \( my $opt_fastANI_bin = `which fastANI 2>/dev/null` );

=over 30

=item B<[--snippy-bin (PATH)]>

Path to the snippy binary file. Default tries if snippy is in PATH;

=back

=cut

$options{'snippy-bin=s'} = \( my $opt_snippy_bin = `which snippy 2>/dev/null` );

=over 30

=item B<[--sickle-bin (PATH)]>

Path to the sickle-trim binary file. Default tries if sickle is in PATH;

=back

=cut

$options{'sickle-bin=s'} = \( my $opt_sickle_bin = `which sickle 2>/dev/null` );

=over 30

=item B<[--mash-bin (PATH)]>

Path to mash binary file. Default tries if mash is in PATH;

=back

=cut

$options{'mash-bin=s'} = \( my $opt_mash_bin = `which mash 2>/dev/null` );

=over 30

=item B<[--abricate-bin (PATH)]>

Path to abricate binary file. Default tries if abricate is in PATH;

=back

=cut

$options{'abricate-bin=s'} = \( my $opt_abricate_bin = `which abricate 2>/dev/null` );

=over 30

=item B<[--unicycler-bin (PATH)]>

Path to unicycler binary file. Default tries if unicycler is in PATH;

=back

=cut

$options{'unicycler-bin=s'} = \( my $opt_unicycler_bin = `which unicycler 2>/dev/null` );

=over 30

=item B<[--muscle-bin (PATH)]>

Path to nuscle binary file. Default tries if muscle in PATH;

=back

=cut

$options{'muscle-bin=s'} = \( my $opt_muscle_bin = `which muscle 2>/dev/null` );

=over 30

=item B<[--trimAL-bin (PATH)]>

Path to trimAL binary file. Default tries if trimAL is in PATH;

=back

=cut

$options{'trimAL-bin=s'} = \( my $opt_trimAL_bin = `which trimal 2>/dev/null` );

=over 30

=item B<[--iqtree-bin (PATH)]>

Path to iqtree binary file. Default tries if iqtree is in PATH;

=back

=cut

$options{'iqtree-bin=s'} = \( my $opt_iqtree_bin = `which iqtree 2>/dev/null` );

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

  Version 1.0.34  Documentation, support and updates available at https://liaochenlanruo.fun/pgcgap

  Author: Hualin Liu

  Contact: liaochenlanruo@webmail.hzau.edu.cn or raise issues at GitHub https://github.com/liaochenlanruo/pgcgap

  Citation: Liu H, Xin B, Zheng J, Zhong H, Yu Y, Peng D, Sun M. Build a bioinformatic analysis platform and apply it to routine analysis of microbial genomics and comparative genomics. Protocol exchange, 2022. DOI: 10.21203/rs.2.21224/v6

=end text

=cut

if ($opt_All or $opt_Assemble or $opt_Annotate or $opt_CoreTree or $opt_Pan or $opt_OrthoF or $opt_ANI or $opt_MASH or $opt_pCOG or $opt_VAR or $opt_AntiRes or $opt_STREE or $opt_ACC) {
	tee STDOUT, ">>$opt_logs";
}

GetOptions(%options) or pod2usage("Try '$0 --help' for more information.");

if($opt_version){
	print RED,"PGCGAP version: " . BOLD, YELLOW, "1.0.34", RESET . "\n";
	print "Enter the command " . BOLD, YELLOW, "pgcgap --check-update", RESET . " to check if there is a new version, and update to the new version if it exists.\n";
	exit 0;
}

if ($opt_help) {
	pod2usage(1);
	exit(0);
}


chomp($opt_sickle_bin, $opt_snippy_bin, $opt_abyss_bin, $opt_canu_bin, $opt_prodigal_bin, $opt_prokka_bin, $opt_cdhit_bin, $opt_mafft_bin, $opt_snpsites_bin, $opt_pal2nal_bin, $opt_roary_bin, $opt_orthofinder_bin, $opt_fastANI_bin, $opt_mash_bin, $opt_abricate_bin, $opt_unicycler_bin, $opt_muscle_bin, $opt_trimAL_bin, $opt_iqtree_bin);
check_external_programs() if($opt_check_external_programs);
check_update() if ($opt_check_update);
pod2usage( -msg => 'cd-hit not in $PATH and binary not specified use --cd-hit-bin', -verbose => 0, -exitval => 1 ) unless ($opt_cdhit_bin);
pod2usage( -msg => 'mafft not in $PATH and binary not specified use --mafft-bin', -verbose => 0, -exitval => 1 ) unless ($opt_mafft_bin);
pod2usage( -msg => 'snp-sites not in $PATH and binary not specified use --snp-sites-bin', -verbose => 0, -exitval => 1 ) unless ($opt_snpsites_bin);
pod2usage( -msg => 'abyss not in $PATH and binary not specified use --abyss-bin', -verbose => 0, -exitval => 1 ) unless ($opt_abyss_bin);
pod2usage( -msg => 'canu not in $PATH and binary not specified use --canu-bin', -verbose => 0, -exitval => 1 ) unless ($opt_canu_bin);
pod2usage( -msg => 'prodigal not in $PATH and binary not specified use --prodigal-bin', -verbose => 0, -exitval => 1 ) unless ($opt_prodigal_bin);
pod2usage( -msg => 'prokka not in $PATH and binary not specified use --prokka-bin', -verbose => 0, -exitval => 1 ) unless ($opt_prokka_bin);
pod2usage( -msg => 'pal2nal.pl not in $PATH and binary not specified use --pal2nal-bin', -verbose => 0, -exitval => 1 ) unless ($opt_pal2nal_bin);
pod2usage( -msg => 'roary not in $PATH and binary not specified use --roary-bin', -verbose => 0, -exitval => 1 ) unless ($opt_roary_bin);
pod2usage( -msg => 'orthofinder not in $PATH and binary not specified use --orthofinder-bin', -verbose => 0, -exitval => 1 ) unless ($opt_orthofinder_bin);
pod2usage( -msg => 'fastANI not in $PATH and binary not specified use --fastANI-bin', -verbose => 0, -exitval => 1 ) unless ($opt_fastANI_bin);
pod2usage( -msg => 'snippy not in $PATH and binary not specified use --snippy-bin', -verbose => 0, -exitval => 1 ) unless ($opt_snippy_bin);
pod2usage( -msg => 'sickle not in $PATH and binary not specified use --sickle-bin', -verbose => 0, -exitval => 1 ) unless ($opt_sickle_bin);
pod2usage( -msg => 'mash not in $PATH and binary not specified use --mash-bin', -verbose => 0, -exitval => 1 ) unless ($opt_mash_bin);
pod2usage( -msg => 'abricate not in $PATH and binary not specified use --abricate-bin', -verbose => 0, -exitval => 1 ) unless ($opt_abricate_bin);
pod2usage( -msg => 'unicycler not in $PATH and binary not specified use --unicycler-bin', -verbose => 0, -exitval => 1 ) unless ($opt_unicycler_bin);
pod2usage( -msg => 'muscle not in $PATH and binary not specified use --muscle-bin', -verbose => 0, -exitval => 1 ) unless ($opt_muscle_bin);
pod2usage( -msg => 'trimAL not in $PATH and binary not specified use --trimAL-bin', -verbose => 0, -exitval => 1 ) unless ($opt_trimAL_bin);
pod2usage( -msg => 'iqtree not in $PATH and binary not specified use --iqtree-bin', -verbose => 0, -exitval => 1 ) unless ($opt_iqtree_bin);


sub check_external_programs{
	my %programs = ("snippy" => $opt_snippy_bin, "abyss" => $opt_abyss_bin, "canu" => $opt_canu_bin, "prodigal" => $opt_prodigal_bin, "prokka" => $opt_prokka_bin, "cd-hit" => $opt_cdhit_bin, "mafft" => $opt_mafft_bin, "snp-sites" => $opt_snpsites_bin, "pal2nal" => $opt_pal2nal_bin, "roary" => $opt_roary_bin, "orthofinder" => $opt_orthofinder_bin, "fastANI" => $opt_fastANI_bin, "mash" => $opt_mash_bin, "abricate" => $opt_abricate_bin, "unicycler" => $opt_unicycler_bin, "muscle" => $opt_muscle_bin, "trimAL" => $opt_trimAL_bin, "iqtree" => $opt_iqtree_bin);
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


sub check_update{
	my $search = `pgcgap --version`;
	$search=~/PGCGAP version: (\S+)/;
	my $current_version = $1;
	if (my $version_conda=`conda search pgcgap`) {
		my @lines = split /\n/, $version_conda;
		$lines[-1]=~/pgcgap\s+(\S+).+/;
		my $latest_version = $1;
		if ($current_version ne $latest_version) {
			print BOLD, RED, "Oh, No! You are running an old version of PGCGAP $current_version, we are going to update to the latest version ",RESET . BOLD, YELLOW, "$latest_version",RESET . BOLD, RED, " now!",RESET . "\n\n";
			print "Please wait patiently, take a break and have a cup of tea or coffee!\n";
			my $installation = `mamba install -y pgcgap=$latest_version`;
			print "$installation\n";
			exit(0);
		}else {
			print BOLD, GREEN, "Congratulations, You are running the latest version of PGCGAP ",RESET . BOLD, YELLOW, "v $latest_version",RESET . ".\n";
			exit(0);
		}
	}
}
#=============================== Get bin PATH ======================================================
my $pgcgap_dir;
my $bin = `which pgcgap`;
if ($bin=~/(.+)\/pgcgap/) {
	$pgcgap_dir = $1;
}


#=============================== setup COG database ================================================
if ($opt_setup_COGdb) {
	#https://ftp.ncbi.nih.gov/pub/COG/COG2020/
#	system("wget -c -r -nH -np -nd -R index.html -P ./ ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/");
#	system("gunzip prot2003-2014.fa.gz");
	system("wget -c --no-check-certificate -P ./ https://bcam.hzau.edu.cn/COGdb2020/cog-20.cog.csv");
	system("wget -c --no-check-certificate -P ./ https://bcam.hzau.edu.cn/COGdb2020/cog-20.def.tab");
	system("wget -c --no-check-certificate -P ./ https://bcam.hzau.edu.cn/COGdb2020/cog-20.fa");
	system("wget -c --no-check-certificate -P ./ https://bcam.hzau.edu.cn/COGdb2020/fun-20.tab");
	system("diamond makedb --in cog-20.fa --db COGdiamond_2020");
	#system("makeblastdb -parse_seqids -in cog-20.fa -input_type fasta -dbtype prot -out COG_2020");
	system("rm cog-20.fa");
	system("mv COGdiamond* cog-20.* fun-20.tab $pgcgap_dir/");
	system("chmod a+x $pgcgap_dir/COGdiamond_2020.dmnd");
	system("chmod a+x $pgcgap_dir/cog-20.*");
	system("chmod a+x $pgcgap_dir/fun-20.tab");
}

if ($opt_setup_COGdb2) {
	#https://ftp.ncbi.nih.gov/pub/COG/COG2020/
#	system("wget -c -r -nH -np -nd -R index.html -P ./ ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/");
	system("wget -c --no-check-certificate -r -nH -np -nd -R index.html -P ./ ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv");
	system("wget -c --no-check-certificate -r -nH -np -nd -R index.html -P ./ ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab");
	system("wget -c --no-check-certificate -r -nH -np -nd -R index.html -P ./ ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz");
	system("wget -c --no-check-certificate -r -nH -np -nd -R index.html -P ./ ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab");
	system("gunzip cog-20.fa.gz");
	system("diamond makedb --in cog-20.fa --db COGdiamond_2020");
	#system("makeblastdb -parse_seqids -in cog-20.fa -input_type fasta -dbtype prot -out COG_2020");
	system("rm cog-20.fa");
	system("mv COGdiamond* cog-20.* fun-20.tab $pgcgap_dir/");
	system("chmod a+x $pgcgap_dir/COGdiamond_2020.dmnd");
	system("chmod a+x $pgcgap_dir/cog-20.*");
	system("chmod a+x $pgcgap_dir/fun-20.tab");
}
#===================================================================================================
my $time_start = $^T;
my $working_dir = getcwd;
#system("mkdir -p Results");

# Phylogenetic tree construction with sequences of single/concatenated dna/protein
if ($opt_STREE) {
	system("mkdir -p $working_dir/Results/STREE");
	my $seqfile = $opt_seqfile;
	$seqfile =~ /(.+\/)*(.+)/;
	my $align_seq = $2 . ".aln";
	my $gblocks_out = $align_seq . ".gb";
	my $seqnum = `grep -c '^>' $seqfile`;
	$seqnum =~ s/[\n\r]+//;
	print "There are $seqnum sequences in the input file\n\n";
	#my $b12 = ceil($seqnum/2) + 1;
	print BOLD, CYAN, "Running muscle for sequence alignment...\n\n", RESET;
	#system("muscle -in $seqfile -out $working_dir/Results/STREE/$align_seq -log $working_dir/Results/STREE/Muscle.LOG"); # muscle < 5.1
	if ($seqnum < 400) {
		system("muscle -align $seqfile -output $working_dir/Results/STREE/$align_seq -threads $opt_threads"); # muscle >= 5.1
	}else {
		system("muscle -super5 $seqfile -output $working_dir/Results/STREE/$align_seq -threads $opt_threads"); # muscle >= 5.1
	}
	print BOLD, CYAN, "Running trimAL for selection of conserved blocks...\n\n", RESET;
	chdir "$working_dir/Results/STREE/";
	system("trimal -in $align_seq -out $gblocks_out -automated1");
	#system("Gblocks $align_seq -t=$opt_seqtype -b1=$b12 -b2=$b12 -b4=5 -b5=h -e=.gb");
	#system("iqtree -s $gblocks_out -nt AUTO -m MFP -mtree -b $opt_bsnum");
	if ($opt_fastboot) {
		print BOLD, CYAN, "Running IQ-TREE for phylogenetic tree construction with the fastboot mode...\n\n", RESET;
		system("iqtree -s $gblocks_out -nt $opt_threads -m MFP -mtree -B $opt_fastboot --wbtl --bnni --safe --keep-ident");
	}else {
		print BOLD, CYAN, "Running IQ-TREE for phylogenetic tree construction with the traditional bootstrap mode...\n\n", RESET;
		system("iqtree -s $gblocks_out -nt $opt_threads -m MFP -mtree -b $opt_bsnum --safe --keep-ident");
	}
	chdir $working_dir;
}

my $cpu_count = `cat /proc/cpuinfo| grep "cpu cores"| uniq`;
$cpu_count =~ /.+?(\d+)/;
my $threads_half = $1;
#my $threads_half = CPU();
#print $threads_half . "\n";
# Genome Assemble with"Abyss" or "Canu"
if ($opt_All or $opt_Assemble) {
	system("mkdir -p Results/Assembles/Scaf");
	system("mkdir -p Results/Assembles/FASTQ_Preprocessor");#2020/4/15
=pod
	my $unicycler_set;
	my $unicycler_set37;
	if ($bin=~/(.+)bin\/pgcgap/) {
		$unicycler_set = $1 . "lib/python3.6/site-packages/unicycler/settings.py";
		$unicycler_set37 = $1 . "lib/python3.7/site-packages/unicycler/settings.py";
		if (-e $unicycler_set) {
			system("sed -i 's/RACON_POLISH_LOOP_COUNT_HYBRID = .*/RACON_POLISH_LOOP_COUNT_HYBRID = 2/g' $unicycler_set");
			system("sed -i 's/RACON_POLISH_LOOP_COUNT_LONG_ONLY = .*/RACON_POLISH_LOOP_COUNT_LONG_ONLY = 4/g' $unicycler_set");
		}elsif (-e $unicycler_set37) {
			system("sed -i 's/RACON_POLISH_LOOP_COUNT_HYBRID = .*/RACON_POLISH_LOOP_COUNT_HYBRID = 2/g' $unicycler_set37");
			system("sed -i 's/RACON_POLISH_LOOP_COUNT_LONG_ONLY = .*/RACON_POLISH_LOOP_COUNT_LONG_ONLY = 4/g' $unicycler_set37");
		}else {
			print "Can not find the unicycler setting file\n";
		}
	}
=cut
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
			system("abyss-pe name=$str k=$opt_kmmer in='$fastp_out1 $fastp_out2' B=2G"); # Bloom filter mode abyss >= 2.3.4
			#system("abyss-pe name=$str k=$opt_kmmer in='$fastp_out1 $fastp_out2' np=$opt_threads");# MPI mode (legacy) abyss < 2.3.4
			print "Assemble complete !\n";
			my $assem = $str . "_assembly";
			system("mkdir -p $working_dir/Results/Assembles/Illumina/$assem");
			my $scaf = $str . "-8.fa";
#			system("mkdir Over");
			system("cp $scaf $working_dir/Results/Assembles/Scaf/Illumina/");
#			system("mv $read1 $read2 Over/");
			system("mv $str*.dot* $str*.fa $str*.path* $str*.dist $str*.fai $str*stats* $str*.hist coverage.hist $str*.tsv $working_dir/Results/Assembles/Illumina/$assem/");
			system("mv $fastp_out1 $fastp_out2 $fastph $fastpj $working_dir/Results/Assembles/FASTQ_Preprocessor");#2020/4/15
		}
		chdir $working_dir;
		#system("realpath $working_dir/Results/Assembles/Scaf/Illumina/* >> scaf.list");
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
		my $time_assemble = time();
		my $time_assemblex = ($time_assemble - $time_start)/3600;
		print "The 'Assemble' program runs for $time_assemblex hours.\n\n";
		chdir $working_dir;
		system("realpath $working_dir/Results/Assembles/Scaf/Illumina/*.filtered.fas >> scaf.list");
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
		#system("realpath $working_dir/Results/Assembles/Scaf/Illumina/* >> scaf.list");
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
		my $time_assemble = time();
		my $time_assemblex = ($time_assemble - $time_start)/3600;
		print "The 'Assemble' program runs for $time_assemblex hours.\n\n";
		chdir $working_dir;
		system("realpath $working_dir/Results/Assembles/Scaf/Illumina/*.filtered.fas >> scaf.list");
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
			system("abyss-pe name=$str k=$opt_kmmer in='$fastp_out1 $fastp_out2' B=2G"); # Bloom filter mode abyss >= 2.3.4
			#system("abyss-pe name=$str k=$opt_kmmer in='$fastp_out1 $fastp_out2' np=$opt_threads");# MPI mode (legacy) abyss < 2.3.4
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
				system("mv $str*.dot* $str*.fa $str*.path* $str*.dist $str*.fai $str*stats* $str*.hist coverage.hist $str*.tsv $working_dir/Results/Assembles/Illumina/$assem/");
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
		#system("realpath $working_dir/Results/Assembles/Scaf/Illumina/* >> scaf.list");
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
		my $time_assemble = time();
		my $time_assemblex = ($time_assemble - $time_start)/3600;
		print "The 'Assemble' program runs for $time_assemblex hours.\n\n";
		chdir $working_dir;
		system("realpath $working_dir/Results/Assembles/Scaf/Illumina/*.filtered.fas >> scaf.list");
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

#=========================================================================================================================================

##Annotate
if ($opt_All or $opt_Annotate) {
	system("mkdir -p Results/Annotations");
	system("mkdir -p Results/Annotations/CDs");
	system("mkdir -p Results/Annotations/AAs");
	system("mkdir -p Results/Annotations/GFF");
	chdir $opt_scafPath;
	my $path = `pwd`;
	print $path . "\n";
	my @scaf = glob("*$opt_Scaf_suffix");
	foreach my $scaf (@scaf) {
		$scaf=~/(.+)$opt_Scaf_suffix/;
		my $str = $1;
		my $faa = $str . ".faa";
		my $fna = $str . ".ffn";
		my $gff = $str . ".gff";
		my $outdir = $str . "_annotation";
		print "Running ORFs finding and annotating...\n";
		system("prokka --force --outdir $outdir --prefix $str --locustag $str --genus $opt_genus --species $opt_species --strain $str --gcode $opt_codon --cpus $opt_threads $scaf");
		system("cp $outdir/$faa $working_dir/Results/Annotations/AAs");
		system("cp $outdir/$fna $working_dir/Results/Annotations/CDs");
		system("cp $outdir/$gff $working_dir/Results/Annotations/GFF");
		system("mv $outdir $working_dir/Results/Annotations/");
	}
	chdir $working_dir;
}

#=========================================================================================================================================

## Wrapper to produce phylogenetic tree from the single core proteins and SNPs tree from the single core genes
if ($opt_All or $opt_CoreTree) {
	my $time_coretrees = time();
	print "Performing --CoreTree function...\n\n";
	system("mkdir -p Results/CoreTrees");
	system("cat $opt_AAsPath/*.faa > All.pep");
	print "Running CD-hit...\n";
	system("cd-hit -i All.pep -o All.pep.nr -c $opt_c -n $opt_n -G $opt_G -T $opt_threads -t $opt_t -aS $opt_aS -aL $opt_aL -g $opt_g -M 0 -d $opt_d");

	print "Starting to extract the list of single copied core proteins...\n\n";
	open CLSTR, "All.pep.nr.clstr" || die;
	open LIST, ">core.pep.list" || die;

	my $count = 0;
	my @array = ();
	my %hashc = ();
	my $locustag = '';
	my $strainame = '';
	my $geneno = 0;

	while (<CLSTR>){
		chomp;
		if (/^>Cluster/){
			foreach my $flag (keys %hashc){
				$count++;
			}
			$geneno = scalar @array;
			if(($count == $opt_strain_num) && ($geneno == $opt_strain_num)){
				print LIST join("\t",@array), "\n";
			}
			$count = 0;
			$geneno = 0;
			@array = ();
			%hashc = ();
		}
			if(/\s+>(\S+.*?)\.\.\./){
				$locustag = $1;
				push @array, $locustag;
				if($locustag =~ /(\S+)_\d+/){
					$strainame = $1;
					$hashc{$strainame}++;
					#$count++;
				}
			}
	}
	close CLSTR;

	foreach my $flag (keys %hashc){
		$count++;
	}
	$geneno = scalar @array;
	if(($count == $opt_strain_num) && ($geneno = $opt_strain_num)){
		print LIST join("\t",@array), "\n";
	}
	close LIST;
	if (-z "core.pep.list") {
		print "No single-copy core proteins were found\n\n";
		exit 1;
	}else {
		print "Starting to extract sequences of single copied core proteins and genes...\n";
	}
	#============================extract ortholog cluster of protein=============================
	open SEQP, "All.pep" || die;
	open OUT, ">All.pep2" || die;
	while (<SEQP>) {
		chomp;
		if (/^(>\S+)/) {
			print OUT $1 . "\n";
		}else {
			print OUT $_ . "\n";
		}
	}
	close SEQP;
	close OUT;

	open IN, "All.pep2" || die;
	local $/ = '>';
	my %hashP = ();

	while(<IN>){
		chomp;
		my ($name, $sequence) = split (/\n/, $_, 2);
		next unless ($name && $sequence);
		my ($n) = $name =~ /^(\S+)/;
		$sequence =~ s/\s+|\n|\-//g;
		$hashP{$n} = $sequence;
	}
	close(IN);

	$/ = "\n";
	open LISTP, "core.pep.list" || die;

	my $dirP = "faa";

	system("mkdir faa");

	my $d = 0;
	my $new_d = 0;

	while(<LISTP>){
		chomp;

		my @array = split (/\n/, $_);
		for my $ele (@array){
			$d++;
			$new_d = sprintf("%04d",$d);
			my @cluster = split (/\s+/, $ele);
			my $fna_file = "OG$new_d".".fa";

			open (OUT, ">$dirP/$fna_file") || die "cannot open $fna_file\n";
			for my $ele (@cluster){
				if(exists $hashP{$ele}){
					print OUT ">$ele\n$hashP{$ele}\n";
				}else{
					warn "error! The gene id $ele is missing in the sequence file.\n";
				}
			}
		}
	}
	close(LISTP);

	
	print "Running mafft...\n\n";
	chdir "faa";
#	system("unset MAFFT_BINARIES");
	my @fa = glob("*.fa");
	foreach (@fa){
		my $name=substr($_,0,(length($_)-3));
		my $in=$name.".fa";
		my $out=$name.".aln";
		system("mafft --quiet --auto --thread $opt_threads $in > $out");
	}


	##==============CONSTRUCT SINGLE CORE PROTEIN TREE========================================================
	#print "Starting to construct single core protein tree...\n\n";
	open CON, ">ALL.core.protein.fasta" || die "open ALL.core.protein.fasta failed\n";
	my $nfilesp = 0; # count number of files
	my %nseqp_hashp = (); # key:infile, val:nseqp
	my %seqid_count_hashp = (); # key:seqid, val:count
	my %HoHp              = ();   #
	my %seqid_HoHp        = ();   #
	my $first_namep       = q{};  # First name in matrix.
	my $lwidthp           = 60;   # default line width for fasta
	my $spacep            = "\t"; # spacepr for aligned print
	my $ncharp            = 0;    # ncharp for phyml header.
	my $nseqp;                    # nseqp for phyml header. Do not initiate!
	my $termp             = $/;   # input record separator
	my @hash_refp_arrayp   = ();   # array with hash references


	my @fasp = glob("*.aln");
	foreach my $argp (@fasp) {
		my $infilep  = $argp;
		my %seq_hashp = parse_fastap($infilep); # key: seqid, value:sequence
		$nfilesp++;

		## Save sequences in array with hash references. Does this work for really large number of fasta files?
		my $hash_refp     = \%seq_hashp;
		push(@hash_refp_arrayp, $hash_refp);

		## Add nseqps to global nseqp_hashp:
		$nseqp_hashp{$infilep} = scalar(keys(%seq_hashp));

		## Get length of sequence for all tax labels. Put in hashes.
		foreach my $tax_keyp (keys %seq_hashp) {
			$seqid_count_hashp{$tax_keyp}++;
			$HoHp{$infilep}{$tax_keyp} = length($seq_hashp{$tax_keyp});
			$seqid_HoHp{$infilep}{$tax_keyp}++;
		}

		## Check all seqs are same length
		my $length;
		my $lnamep;
		foreach my $name (keys %seq_hashp) {
			my $l = length $seq_hashp{$name};
			if (defined $length) {
				if ($length != $l) {
					print STDERR "Error!\nseqpuences in $infilep not all same length ($lnamep is $length, $name is $l)\n";
					exit(1);
				}
			}else {
				$length = length $seq_hashp{$name};
				$lnamep  = $name;
			}
		}
	} # Done with file


	#---------------------------------------------------------------------------
	#  Check if the same number of sequences
	#---------------------------------------------------------------------------
	my $lnamep;
	foreach my $file (keys %nseqp_hashp) {
		my $l = $nseqp_hashp{$file}; # val is a length
		if (defined $nseqp) {
			if ($nseqp != $l) {
				print STDERR "Error!\nNumber of sequences in files differ ($lnamep has $nseqp, $file has $l)\n";
				exit(1);
			}
		}else {
			$nseqp = $nseqp_hashp{$file};
			$lnamep  = $file;
		}
	}


	#---------------------------------------------------------------------------
	#  Check sequence id's
	#---------------------------------------------------------------------------
	if (scalar((keys %seqid_count_hashp)) != $nseqp) { # number of unique seqid's not eq to nseqps
		foreach my $key (sort { $seqid_count_hashp{$b} <=> $seqid_count_hashp{$a} } (keys %seqid_count_hashp)) {
			print STDERR "$key --> $seqid_count_hashp{$key}\n";
		}
		print STDERR "\nError!\nSome sequence labels does not occur in all files.\n";
		print STDERR "That is, sequence id's needs to be identical for concatenation.\n\n";
		exit(1);
	}else {
		## Find the longest taxon name for aligned printing
		my @sorted_names = sort { length($b) <=> length($a) } keys %seqid_count_hashp;
		$spacep = length( shift(@sorted_names) ) + 2;
		$first_namep = $sorted_names[0];
	}


	#---------------------------------------------------------------------------
	#Get ncharp
	#---------------------------------------------------------------------------
	foreach my $h_ref (@hash_refp_arrayp) {
		$ncharp = $ncharp + length($h_ref->{$first_namep});
	}


	#---------------------------------------------------------------------------
	#Print everything to STDOUT
	#---------------------------------------------------------------------------
	print STDERR "\nChecked $nfilesp files -- sequence labels and lengths seems OK.\n";
	print STDERR "Concatenated $nseqp sequences, length $ncharp.\n";
	print STDERR "Printing concatenation to 'ALL.core.protein.fasta'.\n\n";

	##Print the array with hash references (does this work with really large number of files (hashes))?
	##First, concatenate all sequences from hashes
	my %print_hashp = (); # key:label, value:sequence
	foreach my $h_ref (@hash_refp_arrayp) {
		foreach my $seqid (sort keys %$h_ref) {
			$print_hashp{$seqid} .= $h_ref->{$seqid};
		}
	}
	##Then print, and add line breaks in sequences
	foreach my $label (sort keys  %print_hashp) {
		print CON ">$label\n";

		##Print sequence
		##TODO: phylip strict printing of sequence in blocks of 10
		$print_hashp{$label} =~ s/\S{$lwidthp}/$&\n/gs; ## replace word of size $lwidthp with itself and "\n"
		print CON $print_hashp{$label}, "\n";
	}

	print STDERR "Concatenate FASTA alignments to FASTA format completed.\n\n";


	sub parse_fastap {
		my ($infilep) = @_;
		my $termp     = $/; # input record separator;
		my %seq_hashp = (); # key:seqid, val:seq
		open my $INFILEP, "<", $infilep or die "could not open infile '$infilep' : $! \n";
		$/ = ">";
		while(<$INFILEP>) {
		    chomp;
		    next if($_ eq '');
		    my ($id, @sequencelines) = split /\n/;
		    if ($id=~/(^\S+)_\S+$/) {
				$id = $1;
				foreach my $line (@sequencelines) {
					$seq_hashp{$id} .= $line;
				}
		    }
		}
		$/ = $termp;
		return(%seq_hashp);
	} # end of parse_fasta



	#print "\n\n";
	#system("mkdir -p $working_dir/Results/STREE");
	my $seqfile = "ALL.core.protein.fasta";
	#$seqfile =~ /(.+\/)*(.+)/;
	#my $align_seq = $2 . ".aln";
	my $gblocks_out = "ALL.core.protein.fasta.gb";
	my $seqnum = `grep -c '^>' $seqfile`;
	print "There are $seqnum sequences in the input file\n\n";
	my $b12 = ceil($seqnum/2) + 1;
	#print "Running muscle for sequence alignment...\n\n";
	#system("muscle -in $seqfile -out $working_dir/Results/STREE/$align_seq -log $working_dir/Results/STREE/Muscle.LOG");

	print "Running trimAL for selection of conserved blocks...\n\n";
	system("trimal -in $seqfile -out $gblocks_out -automated1");

	#print "Running Gblocks for selection of conserved blocks...\n\n";
	#system("Gblocks $seqfile -t=p -b1=$b12 -b2=$b12 -b4=5 -b5=h -e=.gb");

	print "Constructing ML tree of the single copy core proteins...\n\n";
	#===============================================================================
	if ($opt_fasttree) {
		print "Running FastTree for phylogenetic tree construction...\n\n";
		system("fasttree -quiet $gblocks_out > ALL.core.protein.nwk");
	}else {
		print "Running IQ-TREE for phylogenetic tree construction...\n\n";
		if ($opt_fastboot) {
			system("iqtree -s $gblocks_out -nt $opt_threads -m MFP -mtree -B $opt_fastboot --wbtl --bnni --safe --keep-ident");
		}else {
			system("iqtree -s $gblocks_out -nt $opt_threads -m MFP -mtree -b $opt_bsnum --safe --keep-ident");
		}
	}
	#===================================================================================
	#system("fasttree ALL.core.protein.fasta > ALL.core.protein.nwk");
	print "Constructing single copy core protein tree completed\n\n";
	system("mv ALL.core.protein.* ../Results/CoreTrees/");
	my $time_coretreem = time();
	my $time_coretreep = ($time_coretreem - $time_coretrees)/3600;
	print "The 'CoreTree' program runs for $time_coretreep hours to build single-copy core proteins tree.\n\n";

	#===============================================================================

	chdir "../";

	if ($opt_CDsPath ne "NO") {
		system("cat $opt_CDsPath/*.ffn > All.nuc");
		system("mkdir faa2ffn");
		#============================extract ortholog cluster of genes=============================
		open SEQN, "All.nuc" || die;
		local $/ = '>';
		my %hashN = ();
		<SEQN>;
		while(<SEQN>){
			#local $/ = '>';
			chomp;
			my ($name, $sequence) = split (/\n/, $_, 2);
			next unless ($name && $sequence);
			my ($n) = $name =~ /^(\S+)/;
			$sequence =~ s/\s+|\n|\-//g;
			$hashN{$n} = $sequence;
		}
		close(SEQN);
		$/ = "\n";
		open LISTN, "core.pep.list" || die;

		my $dirN = "ffn";

		system("mkdir ffn");

		my $dn = 0;
		my $new_dn = 0;

		while(<LISTN>){
			chomp;

			my @array = split (/\n/, $_);
			for my $ele (@array){
				$dn++;
				$new_dn = sprintf("%04d",$dn);
				my @cluster = split (/\s+/, $ele);
				my $fna_file = "OG$new_dn".".fa";

				open (OUT, ">$dirN/$fna_file") || die "cannot open $fna_file\n";
				for my $ele (@cluster){
					if(exists $hashN{$ele}){
						print OUT ">$ele\n$hashN{$ele}\n";
					}else{
						warn "error! The gene id is missing in the sequence file.\n";
					}
				}
			}
		}
		close(LISTN);


		opendir(DIR, "faa") || die "Can't open directory\n";
		my @store_array = ();
		@store_array = readdir(DIR);
		my $name = '';

		foreach my $file (@store_array) {
			next unless ($file =~ /^\S+\.aln$/);
			if ($file =~ /^(\S+)\.aln$/){
				$name = $1;
			}

			system("pal2nal.pl faa/$file ffn/$name.fa -nogap -output fasta -codontable $opt_codon > faa2ffn/$name.codon.aln");
		}


		chdir "faa2ffn";

		my @aln=glob("*.aln");
		foreach (@aln){
			my $name=substr($_,0,(length($_)-4));
			my $in=$name.".aln";
			my $out=$name.".fas";
			open ALNIN,"$in" or die;
			open ALNOUT, ">$out" or die;
			while(<ALNIN>){
				chomp;
				if (/(\>\S+)_\S+/){
					print ALNOUT $1."\n";
				}else{
					print ALNOUT $_."\n";
				}
			}
		}
		close ALNIN;
		close ALNOUT;

		open CON, ">ALL.core.nucl.fasta" || die "open ALL.core.nucl.fasta failed\n";
		my $nfiles = 0; # count number of files
		my %nseq_hash = (); # key:infile, val:nseq
		my %seqid_count_hash = (); # key:seqid, val:count
		my %HoH              = ();   #
		my %seqid_HoH        = ();   #
		my $first_name       = q{};  # First name in matrix.
		my $lwidth           = 60;   # default line width for fasta
		my $space            = "\t"; # spacer for aligned print
		my $nchar            = 0;    # nchar for phyml header.
		my $nseq;                    # nseq for phyml header. Do not initiate!
		my $term             = $/;   # input record separator
		my @hash_ref_array   = ();   # array with hash references


		my @fas = glob("*.fas");
		foreach my $arg (@fas) {
			my $infile  = $arg;
			my %seq_hash = parse_fasta($infile); # key: seqid, value:sequence
			$nfiles++;

			## Save sequences in array with hash references. Does this work for really large number of fasta files?
			my $hash_ref     = \%seq_hash;
			push(@hash_ref_array, $hash_ref);

			## Add nseqs to global nseq_hash:
			$nseq_hash{$infile} = scalar(keys(%seq_hash));

			## Get length of sequence for all tax labels. Put in hashes.
			foreach my $tax_key (keys %seq_hash) {
				$seqid_count_hash{$tax_key}++;
				$HoH{$infile}{$tax_key} = length($seq_hash{$tax_key});
				$seqid_HoH{$infile}{$tax_key}++;
			}

			## Check all seqs are same length
			my $length;
			my $lname;
			foreach my $name (keys %seq_hash) {
				my $l = length $seq_hash{$name};
				if (defined $length) {
					if ($length != $l) {
						print STDERR "Error!\nSequences in $infile not all same length ($lname is $length, $name is $l)\n";
						exit(1);
					}
				}
				else {
					$length = length $seq_hash{$name};
					$lname  = $name;
				}
			}
		} # Done with file


		#---------------------------------------------------------------------------
		#  Check if the same number of sequences
		#---------------------------------------------------------------------------
		my $lname;
		foreach my $file (keys %nseq_hash) {
			my $l = $nseq_hash{$file}; # val is a length
			if (defined $nseq) {
				if ($nseq != $l) {
					print STDERR "Error!\nNumber of sequences in files differ ($lname has $nseq, $file has $l)\n";
					exit(1);
				}
			}
			else {
				$nseq = $nseq_hash{$file};
				$lname  = $file;
			}
		}


		#---------------------------------------------------------------------------
		#  Check sequence id's
		#---------------------------------------------------------------------------
		if (scalar((keys %seqid_count_hash)) != $nseq) { # number of unique seqid's not eq to nseqs
			foreach my $key (sort { $seqid_count_hash{$b} <=> $seqid_count_hash{$a} } (keys %seqid_count_hash)) {
				print STDERR "$key --> $seqid_count_hash{$key}\n";
			}
			print STDERR "\nError!\nSome sequence labels does not occur in all files.\n";
			print STDERR "That is, sequence id's needs to be identical for concatenation.\n\n";
			exit(1);
		}
		else {
			## Find the longest taxon name for aligned printing
			my @sorted_names = sort { length($b) <=> length($a) } keys %seqid_count_hash;
			$space = length( shift(@sorted_names) ) + 2;
			$first_name = $sorted_names[0];
		}


		#---------------------------------------------------------------------------
		#Get nchar
		#---------------------------------------------------------------------------
		foreach my $h_ref (@hash_ref_array) {
			$nchar = $nchar + length($h_ref->{$first_name});
		}


		#---------------------------------------------------------------------------
		#Print everything to STDOUT
		#---------------------------------------------------------------------------
		print STDERR "\nChecked $nfiles files -- sequence labels and lengths seems OK.\n";
		print STDERR "Concatenated $nseq sequences, length $nchar.\n";
		print STDERR "Printing concatenation to 'ALL.core.nucl.fasta'.\n\n";

		## Print the array with hash references (does this work with really large number of files (hashes))?
		## First, concatenate all sequences from hashes
		my %print_hash = (); # key:label, value:sequence
		foreach my $h_ref (@hash_ref_array) {
			foreach my $seqid (sort keys %$h_ref) {
				$print_hash{$seqid} .= $h_ref->{$seqid};
			}
		}
		## Then print, and add line breaks in sequences
		foreach my $label (sort keys  %print_hash) {
			print CON ">$label\n";

			## Print sequence
			## TODO: phylip strict printing of sequence in blocks of 10
			$print_hash{$label} =~ s/\S{$lwidth}/$&\n/gs; ## replace word of size $lwidth with itself and "\n"
			print CON $print_hash{$label}, "\n";
		}

		print STDERR "Concatenate FASTA alignments to FASTA format completed.\n\n";


		sub parse_fasta {

			my ($infile) = @_;

			my $term     = $/; # input record separator;
			my %seq_hash = (); # key:seqid, val:seq

			open my $INFILE, "<", $infile or die "could not open infile '$infile' : $! \n";
			$/ = ">";
			while(<$INFILE>) {
				chomp;
				next if($_ eq '');
				my ($id, @sequencelines) = split /\n/;
				foreach my $line (@sequencelines) {
					$seq_hash{$id} .= $line;
				}
			}
			$/ = $term;

			return(%seq_hash);

		} # end of parse_fasta


		print "Calling core SNPs, only output columns containing exclusively ACGT.\n";
		system("snp-sites -o ALL.core.snp.fasta -c ALL.core.nucl.fasta");
		print "Running IQ-TREE for phylogenetic tree construction...\n\n";
		my $fconst = `snp-sites -C ALL.core.nucl.fasta`;
		chomp($fconst);
		if ($opt_fastboot) {
			print "Running IQ-TREE with ultrafast bootstrap $opt_fastboot\n\n";
			system("iqtree -fconst $fconst -s ALL.core.snp.fasta -nt AUTO -m MFP -mtree -B $opt_fastboot --wbtl --bnni --safe --keep-ident");
		}else {
			print "Running IQ-TREE with bootstrap $opt_bsnum\n\n";
			system("iqtree -fconst $fconst -s ALL.core.snp.fasta -nt AUTO -m MFP -mtree -b $opt_bsnum --safe --keep-ident");
		}
		system("mv ALL.core.snp.* ../Results/CoreTrees/");
		#===================end==========================================================

		chdir "../";
		rmove("faa2ffn", "./Results/CoreTrees/faa2ffn");
		rmove("ffn", "./Results/CoreTrees/ffn");
	}
	rmove("faa", "./Results/CoreTrees/faa");
	#system("mv faa ./Results/CoreTrees/");
	system("mv All.* ./Results/CoreTrees/");
	system("mv core.pep.list ./Results/CoreTrees/");
	my $time_coretreed = time();
	my $time_coretree = ($time_coretreed - $time_coretrees)/3600;
	print "The 'CoreTree' program runs for $time_coretree hours.\n\n";
}


if ($opt_All or $opt_Pan) {
	my $time_pans = time();
	print "Performing --Pan function...\n\n";
	my $pangenome = "Results/PanGenome";
	system("roary -p $opt_threads -r -t $opt_codon -i $opt_identi -f $pangenome $opt_GffPath/*.gff");
	chdir $pangenome;
	system("create_pan_genome_plots.R");#create pan genome plots
	system("Rscript $pgcgap_dir/plot_3Dpie.R");#plot pangenome 3D-pie
	system("python $pgcgap_dir/fmplot.py --labels accessory_binary_genes.fa.newick gene_presence_absence.csv");
	chdir $working_dir;
	if ($opt_PanTree) {
		#Constructing Roary single-copy core proteins tree
		system("mkdir $working_dir/$pangenome/Core");
		chdir $opt_GffPath;
		my @gff = glob("*.gff");
		foreach my $gff (@gff) {
			$gff=~/(.+).gff/;
			my $name = $1;
			print $name . "\n";
			system("perl $pgcgap_dir/grep_cds_aas_from_gff3.pl $gff $name");
		}
		system("mv *.id *.cds *.pep $working_dir/$pangenome/Core");
		chdir "$working_dir/$pangenome/Core";
		#chdir $working_dir;
		
		my %hash;
		system("cat *.pep > All_aa.fa");
		#chdir "Results/PanGenome/Core";
		open SEQP, "All_aa.fa" || die;
		open OUT, ">All_aa.fa2" || die;
		while (<SEQP>) {
			chomp;
			if (/^(>\S+)/) {
				print OUT $1 . "\n";
			}else {
				print OUT $_ . "\n";
			}
		}
		close SEQP;
		close OUT;

		local $/ = ">";
		open AA, "All_aa.fa2" || die;
		<AA>;
		while (<AA>) {
			chomp;
			my ($head, $seq) = split "\n", $_, 2;
			$head=~/^(\S+)/;
			$hash{$1} = $seq;
		}

		close AA;

		$/ = "\n";
		open IN, "../gene_presence_absence.csv" || die;
		open TBL, ">gene_presence_absence.tbl" || die;
		while (<IN>) {
			chomp;
			$_=~s/,"/\t/g;
			$_=~s/"//g;
			print TBL $_ . "\n";
		}
		close IN;
		close TBL;


		open INF, "gene_presence_absence.tbl" || die;
		open OUT, ">IDs.txt" || die;
		my $count;
		<INF>;
		while(<INF>){
			chomp;
			my @lines = split /\t/;
			if ($lines[3] == $opt_strain_num && $lines[5] == 1) {
				$count++;
				my $group = "Group_" . $count;
				print OUT $group;
				for (my $i=14; $i<@lines; $i++) {
					$lines[$i]=~/(\S+)/;
					print OUT "\t$1";
				}
				print OUT "\n";
			}
		}
		close INF;
		close OUT;

		open ID, "IDs.txt" || die;
		while (<ID>) {
			chomp;
			my @line = split /\t/;
			my $gene = $line[0] . ".aa";
			open OUTF, ">$gene" || die;
			for (my $j=1; $j<@line; $j++) {
				if (exists $hash{$line[$j]}) {
					print OUTF ">$line[$j]\n$hash{$line[$j]}\n";
				}
			}
			close OUTF;
		}
		close ID;


		print "Running mafft...\n\n";
		my @fa = glob("*.aa");
		foreach (@fa){
			my $name=substr($_,0,(length($_)-3));
			my $in=$name.".aa";
			my $out=$name.".aln";
			system("mafft --quiet --auto --thread $opt_threads $in > $out");
		}

		my @aln = glob("*.aln");
		foreach  (@aln) {
			$_=~/(\S+).aln/;
			my $aa = $1 . ".aa";
			my $file_size = -s $_;
			if ($file_size == 0) {
				system("rm -f $_");
				system("rm -f $aa");
			}
		}
		##==============CONSTRUCT SINGLE CORE PROTEIN TREE========================================================
		#print "Starting to construct single core protein tree...\n\n";
		open CON, ">Roary.core.protein.fasta" || die "open Roary.core.protein.fasta failed\n";
		my $nfilesr = 0; # count number of files
		my %nseqr_hashr = (); # key:infile, val:nseqp
		my %seqid_count_hashr = (); # key:seqid, val:count
		my %HoHr              = ();   #
		my %seqid_HoHr        = ();   #
		my $first_namer       = q{};  # First name in matrix.
		my $lwidthr           = 60;   # default line width for fasta
		my $spacer            = "\t"; # spacepr for aligned print
		my $ncharr            = 0;    # ncharp for phyml header.
		my $nseqr;                    # nseqp for phyml header. Do not initiate!
		my $termr             = $/;   # input record separator
		my @hash_refr_arrayr   = ();   # array with hash references


		my @fasr = glob("*.aln");
		foreach my $argr (@fasr) {
			my $infiler  = $argr;
			my %seq_hashr = parse_fastar($infiler); # key: seqid, value:sequence
			$nfilesr++;

			## Save sequences in array with hash references. Does this work for really large number of fasta files?
			my $hash_refr     = \%seq_hashr;
			push(@hash_refr_arrayr, $hash_refr);

			## Add nseqps to global nseqp_hashp:
			$nseqr_hashr{$infiler} = scalar(keys(%seq_hashr));

			## Get length of sequence for all tax labels. Put in hashes.
			foreach my $tax_keyr (keys %seq_hashr) {
				$seqid_count_hashr{$tax_keyr}++;
				$HoHr{$infiler}{$tax_keyr} = length($seq_hashr{$tax_keyr});
				$seqid_HoHr{$infiler}{$tax_keyr}++;
			}

			## Check all seqs are same length
			my $length;
			my $lnamer;
			foreach my $name (keys %seq_hashr) {
				my $l = length $seq_hashr{$name};
				if (defined $length) {
					if ($length != $l) {
						print STDERR "Error!\nseqpuences in $infiler not all same length ($lnamer is $length, $name is $l)\n";
						exit(1);
					}
				}else {
					$length = length $seq_hashr{$name};
					$lnamer  = $name;
				}
			}
		} # Done with file


		#---------------------------------------------------------------------------
		#  Check if the same number of sequences
		#---------------------------------------------------------------------------
		my $lnamer;
		foreach my $file (keys %nseqr_hashr) {
			my $l = $nseqr_hashr{$file}; # val is a length
			if (defined $nseqr) {
				if ($nseqr != $l) {
					print STDERR "Error!\nNumber of sequences in files differ ($lnamer has $nseqr, $file has $l)\n";
					exit(1);
				}
			}else {
				$nseqr = $nseqr_hashr{$file};
				$lnamer  = $file;
			}
		}


		#---------------------------------------------------------------------------
		#  Check sequence id's
		#---------------------------------------------------------------------------
		if (scalar((keys %seqid_count_hashr)) != $nseqr) { # number of unique seqid's not eq to nseqrs
			foreach my $key (sort { $seqid_count_hashr{$b} <=> $seqid_count_hashr{$a} } (keys %seqid_count_hashr)) {
				print STDERR "$key --> $seqid_count_hashr{$key}\n";
			}
			print STDERR "\nError!\nSome sequence labels does not occur in all files.\n";
			print STDERR "That is, sequence id's needs to be identical for concatenation.\n\n";
			exit(1);
		}else {
			## Find the longest taxon name for aligned printing
			my @sorted_names = sort { length($b) <=> length($a) } keys %seqid_count_hashr;
			$spacer = length( shift(@sorted_names) ) + 2;
			$first_namer = $sorted_names[0];
		}


		#---------------------------------------------------------------------------
		#Get ncharp
		#---------------------------------------------------------------------------
		foreach my $h_ref (@hash_refr_arrayr) {
			$ncharr = $ncharr + length($h_ref->{$first_namer});
		}


		#---------------------------------------------------------------------------
		#Print everything to STDOUT
		#---------------------------------------------------------------------------
		print STDERR "\nChecked $nfilesr files -- sequence labels and lengths seems OK.\n";
		print STDERR "Concatenated $nseqr sequences, length $ncharr.\n";
		print STDERR "Printing concatenation to 'Roary.core.protein.fasta'.\n\n";

		##Print the array with hash references (does this work with really large number of files (hashes))?
		##First, concatenate all sequences from hashes
		my %print_hashr = (); # key:label, value:sequence
		foreach my $h_ref (@hash_refr_arrayr) {
			foreach my $seqid (sort keys %$h_ref) {
				$print_hashr{$seqid} .= $h_ref->{$seqid};
			}
		}
		##Then print, and add line breaks in sequences
		foreach my $label (sort keys  %print_hashr) {
			print CON ">$label\n";

			##Print sequence
			##TODO: phylip strict printing of sequence in blocks of 10
			$print_hashr{$label} =~ s/\S{$lwidthr}/$&\n/gs; ## replace word of size $lwidthr with itself and "\n"
			print CON $print_hashr{$label}, "\n";
		}

		print STDERR "Concatenate FASTA alignments to FASTA format completed.\n\n";


		sub parse_fastar {
			my ($infiler) = @_;
			my $termp     = $/; # input record separator;
			my %seq_hashr = (); # key:seqid, val:seq
			open my $INFILER, "<", $infiler or die "could not open infile '$infiler' : $! \n";
			$/ = ">";
			while(<$INFILER>) {
				chomp;
				next if($_ eq '');
				my ($id, @sequencelines) = split /\n/;
				if ($id=~/(^\S+)_\S+$/) {
					$id = $1;
					foreach my $line (@sequencelines) {
						$seq_hashr{$id} .= $line;
					}
				}
			}
			$/ = $termr;
			return(%seq_hashr);
		} # end of parse_fastar

		my $seqfilen = "Roary.core.protein.fasta";
		my $gblocks_outn = "Roary.core.protein.fasta.gb";
		my $seqnumn = `grep -c '^>' $seqfilen`;
		print "There are $seqnumn sequences in the input file\n\n";
		my $b12n = ceil($seqnumn/2) + 1;

		print "Running trimAL for selection of conserved blocks...\n\n";
		system("trimal -in $seqfilen -out $gblocks_outn -automated1");

		#print "Running Gblocks for selection of conserved blocks...\n\n";
		#system("Gblocks $seqfilen -t=p -b1=$b12n -b2=$b12n -b4=5 -b5=h -e=.gb");
		print "Constructing ML tree of the single-copy core proteins...\n\n";
		#===============================================================================
		if ($opt_fasttree) {
			print "Running FastTree for phylogenetic tree construction...\n\n";
			system("fasttree -quiet $gblocks_outn > Roary.core.protein.nwk");
		}else {
			print "Running IQ-TREE for phylogenetic tree construction...\n\n";
			if ($opt_fastboot) {
				system("iqtree -s $gblocks_outn -nt $opt_threads -m MFP -mtree -B $opt_fastboot --wbtl --bnni --safe --keep-ident");
			}else {
				system("iqtree -s $gblocks_outn -nt $opt_threads -m MFP -mtree -b $opt_bsnum --safe --keep-ident");
			}
		}
		#===================================================================================
		#system("fasttree Roary.core.protein.fasta > Roary.core.protein.nwk");
		print "Constructing single-copy core protein tree completed\n\n";
	}
	my $time_pand = time();
	my $time_pan = ($time_pand - $time_pans)/3600;
	print "The 'Pan' program runs for $time_pan hours.\n\n";
	chdir $working_dir;
}

if ($opt_All or $opt_OrthoF) {
	my $time_OrthoFs = time();
	print "Performing --OrthoF function...\n\n";
	#system("mkdir Results/OrthoF");
	my $orthoFprefix = "orthoF";
	system("orthofinder -a $opt_threads -t $opt_threads -S $opt_Sprogram -n $orthoFprefix -f $opt_AAsPath");
	if (-e (glob("$working_dir/$opt_AAsPath/OrthoFinder/Results_orthoF/Single_Copy_Orthologue_Sequences/*.fa"))[0]) {
		system("mkdir -p $working_dir/$opt_AAsPath/OrthoFinder/Results_orthoF/Single_Copy_Orthologue_Tree");
		chdir "$working_dir/$opt_AAsPath/OrthoFinder/Results_orthoF/Single_Copy_Orthologue_Sequences";
		print "Running mafft...\n\n";
	#	system("unset MAFFT_BINARIES");
		my @fa = glob("*.fa");
		foreach (@fa){
			my $name=substr($_,0,(length($_)-3));
			my $in=$name.".fa";
			my $out=$name.".aln";
			system("mafft --quiet --auto --thread $opt_threads $in > $out");
		}


		##==============CONSTRUCT Single Copy Orthologue TREE========================================================
		#print "Starting to construct Single Copy Orthologue tree...\n\n";
		open CON, ">Single.Copy.Orthologue.fasta" || die "Open Single.Copy.Orthologue.fasta failed\n";
		my $nfilesp = 0; # count number of files
		my %nseqp_hashp = (); # key:infile, val:nseqp
		my %seqid_count_hashp = (); # key:seqid, val:count
		my %HoHp              = ();   #
		my %seqid_HoHp        = ();   #
		my $first_namep       = q{};  # First name in matrix.
		my $lwidthp           = 60;   # default line width for fasta
		my $spacep            = "\t"; # spacepr for aligned print
		my $ncharp            = 0;    # ncharp for phyml header.
		my $nseqp;                    # nseqp for phyml header. Do not initiate!
		my $termp             = $/;   # input record separator
		my @hash_refp_arrayp   = ();   # array with hash references


		my @fasp = glob("*.aln");
		foreach my $argp (@fasp) {
			my $infilep  = $argp;
			my %seq_hashp = parse_fastao($infilep); # key: seqid, value:sequence
			$nfilesp++;

			## Save sequences in array with hash references. Does this work for really large number of fasta files?
			my $hash_refp     = \%seq_hashp;
			push(@hash_refp_arrayp, $hash_refp);

			## Add nseqps to global nseqp_hashp:
			$nseqp_hashp{$infilep} = scalar(keys(%seq_hashp));

			## Get length of sequence for all tax labels. Put in hashes.
			foreach my $tax_keyp (keys %seq_hashp) {
				$seqid_count_hashp{$tax_keyp}++;
				$HoHp{$infilep}{$tax_keyp} = length($seq_hashp{$tax_keyp});
				$seqid_HoHp{$infilep}{$tax_keyp}++;
			}

			## Check all seqs are same length
			my $length;
			my $lnamep;
			foreach my $name (keys %seq_hashp) {
				my $l = length $seq_hashp{$name};
				if (defined $length) {
					if ($length != $l) {
						print STDERR "Error!\nseqpuences in $infilep not all same length ($lnamep is $length, $name is $l)\n";
						exit(1);
					}
				}else {
					$length = length $seq_hashp{$name};
					$lnamep  = $name;
				}
			}
		} # Done with file


		#---------------------------------------------------------------------------
		#  Check if the same number of sequences
		#---------------------------------------------------------------------------
		my $lnamep;
		foreach my $file (keys %nseqp_hashp) {
			my $l = $nseqp_hashp{$file}; # val is a length
			if (defined $nseqp) {
				if ($nseqp != $l) {
					print STDERR "Error!\nNumber of sequences in files differ ($lnamep has $nseqp, $file has $l)\n";
					exit(1);
				}
			}else {
				$nseqp = $nseqp_hashp{$file};
				$lnamep  = $file;
			}
		}


		#---------------------------------------------------------------------------
		#  Check sequence id's
		#---------------------------------------------------------------------------
		if (scalar((keys %seqid_count_hashp)) != $nseqp) { # number of unique seqid's not eq to nseqps
			foreach my $key (sort { $seqid_count_hashp{$b} <=> $seqid_count_hashp{$a} } (keys %seqid_count_hashp)) {
				print STDERR "$key --> $seqid_count_hashp{$key}\n";
			}
			print STDERR "\nError!\nSome sequence labels does not occur in all files.\n";
			print STDERR "That is, sequence id's needs to be identical for concatenation.\n\n";
			exit(1);
		}else {
			## Find the longest taxon name for aligned printing
			my @sorted_names = sort { length($b) <=> length($a) } keys %seqid_count_hashp;
			$spacep = length( shift(@sorted_names) ) + 2;
			$first_namep = $sorted_names[0];
		}


		#---------------------------------------------------------------------------
		#Get ncharp
		#---------------------------------------------------------------------------
		foreach my $h_ref (@hash_refp_arrayp) {
			$ncharp = $ncharp + length($h_ref->{$first_namep});
		}


		#---------------------------------------------------------------------------
		#Print everything to STDOUT
		#---------------------------------------------------------------------------
		print STDERR "\nChecked $nfilesp files -- sequence labels and lengths seems OK.\n";
		print STDERR "Concatenated $nseqp sequences, length $ncharp.\n";
		print STDERR "Printing concatenation to 'Single.Copy.Orthologue.fasta'.\n\n";

		##Print the array with hash references (does this work with really large number of files (hashes))?
		##First, concatenate all sequences from hashes
		my %print_hashp = (); # key:label, value:sequence
		foreach my $h_ref (@hash_refp_arrayp) {
			foreach my $seqid (sort keys %$h_ref) {
				$print_hashp{$seqid} .= $h_ref->{$seqid};
			}
		}
		##Then print, and add line breaks in sequences
		foreach my $label (sort keys  %print_hashp) {
			print CON ">$label\n";

			##Print sequence
			##TODO: phylip strict printing of sequence in blocks of 10
			$print_hashp{$label} =~ s/\S{$lwidthp}/$&\n/gs; ## replace word of size $lwidthp with itself and "\n"
			print CON $print_hashp{$label}, "\n";
		}

		print STDERR "Concatenate FASTA alignments to FASTA format completed.\n\n";


		sub parse_fastao {
			my ($infilep) = @_;
			my $termp     = $/; # input record separator;
			my %seq_hashp = (); # key:seqid, val:seq
			open my $INFILEP, "<", $infilep or die "could not open infile '$infilep' : $! \n";
			$/ = ">";
			while(<$INFILEP>) {
				chomp;
				next if($_ eq '');
				my ($id, @sequencelines) = split /\n/;
				if ($id=~/(^\S+)_\S+$/) {
					$id = $1;
					foreach my $line (@sequencelines) {
						$seq_hashp{$id} .= $line;
					}
				}
			}
			$/ = $termp;
			return(%seq_hashp);
		} # end of parse_fasta



		my $seqfile = "Single.Copy.Orthologue.fasta";
		my $gblocks_out = "Single.Copy.Orthologue.fasta.gb";
		my $seqnum = `grep -c '^>' $seqfile`;
		print "There are $seqnum sequences in the input file\n\n";
		my $b12 = ceil($seqnum/2) + 1;
		print "Running trimAL for selection of conserved blocks...\n\n";
		system("trimal -in $seqfile -out $gblocks_out -automated1");
		print "Constructing ML tree of the Single Copy Orthologue proteins...\n\n";
		#===============================================================================
		if ($opt_fasttree) {
			print "Running FastTree for phylogenetic tree construction...\n\n";
			system("fasttree -quiet $gblocks_out > Single.Copy.Orthologue.nwk");
		}else {
			print "Running IQ-TREE for phylogenetic tree construction...\n\n";
			if ($opt_fastboot) {
				system("iqtree -s $gblocks_out -nt $opt_threads -m MFP -mtree -B $opt_fastboot --wbtl --bnni --safe --keep-ident");
			}else {
				system("iqtree -s $gblocks_out -nt $opt_threads -m MFP -mtree -b $opt_bsnum --safe --keep-ident");
			}
		}
		#===================================================================================
		print "Constructing single copy Orthologue protein tree completed\n\n";
		system("mv Single.Copy.Orthologue.* ../Single_Copy_Orthologue_Tree/");
		#===============================================================================
	}
	chdir $working_dir;
	rmove("$opt_AAsPath/OrthoFinder/", "Results/OrthoFinder");
	my $time_OrthoFd = time();
	my $time_OrthoF = ($time_OrthoFd - $time_OrthoFs)/3600;
	print "The 'OrthoF' program runs for $time_OrthoF hours.\n\n";
	chdir $working_dir;
}

if ($opt_All or $opt_MASH) {
	my $time_MASHs = time();
	print "Performing --MASH function...\n\n";
	system("mkdir -p Results/MASH");
	chdir $opt_scafPath;
	my @genome = glob("*$opt_Scaf_suffix");
	foreach  (@genome) {
		system("mash sketch $_");
	}

	my @msh = glob("*.msh");
	for (my $i=0; $i<@msh; $i++) {
		for (my $j=0; $j<@msh; $j++) {
			system("mash dist $msh[$i] $msh[$j] >> MASH");
		}
	}
	open IN, "MASH" || die;
	open OUT, ">MASH2" || die;
	while (<IN>) {
		chomp;
		my @lines = split /\t/;
		my $dis = (1-$lines[2])*100;
		print OUT "$lines[0]\t$lines[1]\t$dis\t$lines[3]\t$lines[4]\n";
	}
	system("perl $pgcgap_dir/get_Mash_Matrix.pl --Scaf_suffix $opt_Scaf_suffix");
	system("Rscript $pgcgap_dir/Plot_MashHeatmap.R");
	system("rm -f *.msh");
	system("mv MASH MASH2 MASH.heatmap MASH_matrix.pdf $working_dir/Results/MASH");
	chdir $working_dir;
	my $time_MASHd = time();
	my $time_MASH = ($time_MASHd - $time_MASHs)/3600;
	print "The 'MASH' program runs for $time_MASH hours.\n\n";
}

if ($opt_All or $opt_ANI) {
	my $time_ANIs = time();
	print "Performing --ANI function...\n\n";
	system("mkdir -p Results/ANI");
	my $ANIO = "Results/ANI/ANIs";
	system("fastANI --matrix -t $opt_threads --ql $opt_queryL --rl $opt_refL -o $ANIO");
	chdir "Results/ANI";
	system("perl $pgcgap_dir/triangle2list.pl");
	system("perl $pgcgap_dir/get_ANImatrix.pl --Scaf_suffix $opt_Scaf_suffix");
	system("Rscript $pgcgap_dir/Plot_ANIheatmap.R");
	chdir $working_dir;
	my $time_ANId = time();
	my $time_ANI = ($time_ANId - $time_ANIs)/3600;
	print "The 'ANI' program runs for $time_ANI hours.\n\n";
}

if ($opt_VAR) {
	my $time_VARs = time();
	print "Performing --VAR function...\n\n";
	system("mkdir -p Results/Variants");
	chdir "$opt_ReadsPath";
	system("mkdir Trimmed");
	system("cp $opt_refgbk ./");
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
		my $trif = $str . "_trimmed_1.fastq";
		my $trir = $str . "_trimmed_2.fastq";
		my $tris = $str . "_trimmed_s.fastq";
		system("sickle pe -f $read1 -r $read2 -t $opt_qualtype -o Trimmed/$trif -p Trimmed/$trir -s Trimmed/$tris -q $opt_qual -l $opt_length");#Quality trimming
		system("snippy --cpus $opt_threads --ram $opt_ram --prefix $str --mincov $opt_mincov --minfrac $opt_minfrac --minqual $opt_minqual --outdir $working_dir/Results/Variants/$str --ref $opt_refgbk --R1 Trimmed/$trif --R2 Trimmed/$trir --report");
	}
	chdir $working_dir;
	system("snippy-core --ref $opt_refgbk $working_dir/Results/Variants/*");
	#system("snp-sites -c core.full.aln -o core.full.ATGC.aln");
	system("mkdir -p Results/Variants/Core");
	my $fconst = `snp-sites -C core.full.aln`;
	chomp($fconst);
	if ($opt_strain_num > 2) {
		if ($opt_fastboot) {
			print "Running IQ-TREE with ultrafast bootstrap $opt_fastboot\n\n";
			system("iqtree -fconst $fconst -s core.aln -nt AUTO -m MFP -mtree -B $opt_fastboot --wbtl --bnni --safe --keep-ident");
		}else {
			print "Running IQ-TREE with bootstrap $opt_bsnum\n\n";
			system("iqtree -fconst $fconst -s core.aln -nt AUTO -m MFP -mtree -b $opt_bsnum --safe --keep-ident");
		}
	}
	chdir $working_dir;
	system("mv core.* Results/Variants/Core/");
	#===================end==========================================================
	my $time_VARd = time();
	my $time_VAR = ($time_VARd - $time_VARs)/3600;
	print "The 'VAR' program runs for $time_VAR hours.\n\n";
}

if ($opt_All or $opt_AntiRes) {
	my $time_Antis = time();
	print "Performing --AntiRes function...\n\n";
	system("mkdir -p Results/AntiRes");
	chdir $opt_scafPath;
	my @genome = glob("*$opt_Scaf_suffix");
	foreach my $genome (@genome) {
		$genome=~/(.+)$opt_Scaf_suffix/;
		my $str = $1;
		if ($opt_db eq "all") {
			my @db = ("argannot", "card", "ecoh", "ecoli_vf", "ncbi", "plasmidfinder", "resfinder", "vfdb", "megares");
			foreach my $db (@db) {
				my $out = $str . "_" . $db . ".tab";
				system("abricate --threads $opt_threads --db $db --minid $opt_identity --mincov $opt_coverage $genome > $working_dir/Results/AntiRes/$out");
			}
		}else {
			my $out = $str . "_" . $opt_db . ".tab";
			system("abricate --threads $opt_threads --db $opt_db --minid $opt_identity --mincov $opt_coverage $genome > $working_dir/Results/AntiRes/$out");
		}
	}
	chdir $working_dir;
	system("abricate --summary Results/AntiRes/*.tab > Results/AntiRes/summary.txt");
	my $time_Antid = time();
	my $time_Anti = ($time_Antid - $time_Antis)/3600;
	print "The 'AntiRes' program runs for $time_Anti hours.\n\n";
}

if ($opt_All or $opt_pCOG) {
	my $time_COGs = time();
	print "Performing --COG function...\n\n";
	system("mkdir -p Results/COG");
	system("COGdiamond2022.pl --threads $opt_threads --strain_num $opt_strain_num --evalue $opt_evalue --id $opt_id --query_cover $opt_query_cover --subject_cover $opt_subject_cover --AAsPath $opt_AAsPath");
	system("mv $opt_AAsPath/*.table $opt_AAsPath/*.pdf $opt_AAsPath/*.xml $working_dir/Results/COG");
	chdir $working_dir;
	my $time_COGd = time();
	my $time_COG = ($time_COGd - $time_COGs)/3600;
	print "The 'pCOG' program runs for $time_COG hours.\n\n";
}

my $time_end = time();
my $time_total = ($time_end - $time_start)/3600;
#print "Total $time_total hours used.\n\n";

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

sub printAnnotate{
	print "[--scafPath (PATH)] Path for contigs/scaffolds ( Default 'Results/Assembles/Scaf/Illumina' )\n";
	print "[--Scaf_suffix (STRING)] The suffix of scaffolds or genome files. Users should set the suffixes according to the actual situation ( Default .filtered.fas )\n";
	print "[--codon (INT)] Translation table ( Default 11 )\n  1   Universal code\n  2   Vertebrate mitochondrial code\n  3   Yeast mitochondrial code\n  4   Mold, Protozoan, and Coelenterate Mitochondrial code and Mycoplasma/Spiroplasma code\n  5   Invertebrate mitochondrial\n  6   Ciliate, Dasycladacean and Hexamita nuclear code\n  9   Echinoderm and Flatworm mitochondrial code\n  10  Euplotid nuclear code\n  11  Bacterial, archaeal and plant plastid code ( Default )\n  12  Alternative yeast nuclear code\n  13  Ascidian mitochondrial code\n  14  Alternative flatworm mitochondrial code\n  15  Blepharisma nuclear code\n  16  Chlorophycean mitochondrial code\n  21  Trematode mitochondrial code\n  22  Scenedesmus obliquus mitochondrial code\n  23  Thraustochytrium mitochondrial code\n";
	print "[--genus (STRING)] Genus name of the strain ( Default 'NA' )\n";
	print "[--species (STRING)] Species name of the strain ( Default 'NA' )\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
}

sub printCoreTree{
	print "[--AAsPath (PATH)] Amino acids of all strains as fasta file paths, ( Default './Results/Annotations/AAs' )\n";
	print "[--CDsPath (PATH)] CDs of all strains as fasta file paths, ( Default './Results/Annotations/CDs' )\n";
	print "[--strain_num (INT)] The total number of strains used for analysis, not including the reference genome\n";
	print "[--codon (INT)] Translation table ( Default 11 )\n  1   Universal code\n  2   Vertebrate mitochondrial code\n  3   Yeast mitochondrial code\n  4   Mold, Protozoan, and Coelenterate Mitochondrial code and Mycoplasma/Spiroplasma code\n  5   Invertebrate mitochondrial\n  6   Ciliate, Dasycladacean and Hexamita nuclear code\n  9   Echinoderm and Flatworm mitochondrial code\n  10  Euplotid nuclear code\n  11  Bacterial, archaeal and plant plastid code ( Default )\n  12  Alternative yeast nuclear code\n  13  Ascidian mitochondrial code\n  14  Alternative flatworm mitochondrial code\n  15  Blepharisma nuclear code\n  16  Chlorophycean mitochondrial code\n  21  Trematode mitochondrial code\n  22  Scenedesmus obliquus mitochondrial code\n  23  Thraustochytrium mitochondrial code\n";
	print "[-c (FLOAT)] Sequence identity threshold, ( Default 0.5)\n";
	print "[-n (INT)] Word_length, -n 2 for thresholds 0.4-0.5, -n 3 for thresholds 0.5-0.6, -n 4 for thresholds 0.6-0.7, -n 5 for thresholds 0.7-1.0 ( Default 2 )\n";
	print "[-G (INT)] Use global (set to 1) or local (set to 0) sequence identity, ( Default 0 )\n";
	print "[-t (INT)] Tolerance for redundance ( Default 0 )\n";
	print "[-aL (FLOAT)] Alignment coverage for the longer sequence. If set to 0.9, the alignment must covers 90% of the sequence ( Default 0.5 )\n";
	print "[-aS (FLOAT)] Alignment coverage for the shorter sequence. If set to 0.9, the alignment must covers 90% of the sequence ( Default 0.7 )\n";
	print "[-g (INT)] If set to 0, a sequence is clustered to the first cluster that meets the threshold (fast cluster). If set to 1, the program will cluster it into the most similar cluster that meets the threshold (accurate but slow mode, Default 1)\n";
	print "[-d (INT)] length of description in .clstr file. if set to 0, it takes the fasta defline and stops at first space ( Default 0 )\n";
	print "[--fasttree] Use FastTree to construct phylogenetic tree quickly instead of the combination of Modeltest-ng and Raxml-ng\n";
	print "[--bsnum (INT)] Replicates for bootstrap of IQ-TREE ( Default 500 )\n\n";
	print "[--fastboot (INT)] Replicates for ultrafast bootstrap of IQ-TREE ( must >=1000, Default 1000 )\n\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
}

sub printPan{
	print "[--GffPath (PATH)] Gff files of all strains as paths ( Default './Results/Annotations/GFF' )\n";
	print "[--codon (INT)] Translation table ( Default 11 )\n  1   Universal code\n  2   Vertebrate mitochondrial code\n  3   Yeast mitochondrial code\n  4   Mold, Protozoan, and Coelenterate Mitochondrial code and Mycoplasma/Spiroplasma code\n  5   Invertebrate mitochondrial\n  6   Ciliate, Dasycladacean and Hexamita nuclear code\n  9   Echinoderm and Flatworm mitochondrial code\n  10  Euplotid nuclear code\n  11  Bacterial, archaeal and plant plastid code ( Default )\n  12  Alternative yeast nuclear code\n  13  Ascidian mitochondrial code\n  14  Alternative flatworm mitochondrial code\n  15  Blepharisma nuclear code\n  16  Chlorophycean mitochondrial code\n  21  Trematode mitochondrial code\n  22  Scenedesmus obliquus mitochondrial code\n  23  Thraustochytrium mitochondrial code\n";
	print "[--strain_num (INT)] The total number of strains used for analysis, not including the reference genome\n";
	print "[--PanTree] Construct a phylogenetic tree of single-copy core proteins called by roary\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
	print "[--identi (INT)] Minimum percentage identity for blastp ( Default 95 )\n";
	print "[--fasttree] Use FastTree to construct phylogenetic tree quickly instead of the combination of Modeltest-ng and Raxml-ng\n";
	print "[--bsnum (INT)] Replicates for bootstrap of IQ-TREE ( Default 500 )\n\n";
	print "[--fastboot (INT)] Replicates for ultrafast bootstrap of IQ-TREE ( must >=1000, Default 1000 )\n\n";
}

sub printOrthoF{
	print "[--AAsPath (PATH)] Amino acids of all strains as fasta file paths, ( Default './Results/Annotations/AAs' )\n";
	print "[--Sprogram (STRING)] Sequence search program, Options: blast, mmseqs, blast_gz, diamond ( Default diamond )\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
	print "[--fasttree] Use FastTree to construct phylogenetic tree quickly instead of the combination of Modeltest-ng and Raxml-ng\n";
	print "[--bsnum (INT)] Replicates for bootstrap of IQ-TREE ( Default 500 )\n\n";
	print "[--fastboot (INT)] Replicates for ultrafast bootstrap of IQ-TREE ( must >=1000, Default 1000 )\n\n";
}

sub printANI{
	print "[--queryL (FILE)] The file containing full paths to query genomes, one per line ( Default scaf.list )\n";
	print "[--refL (FILE)] The file containing full paths to reference genomes, one per line. ( Default scaf.list )\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
}

sub printMASH{
	print "[--scafPath (PATH)] Path for contigs/scaffolds ( Default 'Results/Assembles/Scaf/Illumina' )\n";
	print "[--Scaf_suffix (STRING)] The suffix of scaffolds or genome files. Users should set the suffixes according to the actual situation ( Default .filtered.fas )\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
}

sub printVAR{
	print "[--ReadsPath (PATH)] Reads of all strains as file paths ( Default ./Reads/Illumina )\n";
	print "[--reads1 (STRING)] The suffix name of reads 1 ( for example: if the name of reads 1 is 'YBT-1520_L1_I050.R1.clean.fastq.gz', 'YBT-1520' is the strain same, the suffix name should be '.R1.clean.fastq.gz' )\n";
	print "[--reads2 (STRING)] The suffix name of reads 2( for example: if the name of reads 2 is 'YBT-1520_2.fq', the suffix name should be '_2.fq' )\n";
	print "[--refgbk (FILE)] The B<full path and name> of reference genome in GENBANK format ( B<recommended> ), fasta format is also OK. For example: '/mnt/g/test/ref.gbk'\n";
	print "[--suffix_len (INT)] (Strongly recommended) The suffix length of the reads file, that is the length of the reads name minus the length of the strain name. For example the --suffix_len of 'YBT-1520_L1_I050.R1.clean.fastq.gz' is 26 ( 'YBT-1520' is the strain name ) ( Default 0 )\n";
	print "[--strain_num (INT)] The total number of strains used for analysis, not including the reference genome\n";
	print "[--qualtype (STRING)] Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)). ( Default sanger )\n";
	print "[--qual (INT)] Threshold for trimming based on average quality in a window. ( Default 20 )\n";
	print "[--length (INT)] Threshold to keep a read based on length after trimming. ( Default 20 )\n";
	print "[--mincov (INT)] The minimum number of reads covering a site to be considered ( Default 10 )\n";
	print "[--minfrac (FLOAT)] The minimum proportion of those reads which must differ from the reference ( Default 0.9 )\n";
	print "[--minqual (INT)] The minimum VCF variant call 'quality' ( Default 100 )\n";
	print "[--ram (INT)] Try and keep RAM under this many GB ( Default 8 )\n";
#	print "[--tree_builder (STRING)] Application to use for tree building [raxml|fasttree|hybrid] ( Default fasttree )\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
#	print "[--iterations (INT)] Maximum No. of iterations for gubbins ( Default 5 )\n";
	print "[--bsnum (INT)] Replicates for bootstrap of IQ-TREE ( Default 500 )\n\n";
	print "[--fastboot (INT)] Replicates for ultrafast bootstrap of IQ-TREE ( must >=1000, Default 1000 )\n\n";
}

sub printAntiRes{
	print "[--scafPath (PATH)] Path for contigs/scaffolds ( Default 'Results/Assembles/Scaf/Illumina' )\n";
	print "[--Scaf_suffix (STRING)] The suffix of scaffolds or genome files. Users should set the suffixes according to the actual situation ( Default .filtered.fas )\n";
	print "[--db (STRING)]> The database to use, options: all, argannot, card, ecoh, ecoli_vf, megares, ncbi, plasmidfinder, resfinder and vfdb. ( Default all )\n";
	print "[--identity (INT)] Minimum %identity to keep the result, should be a number between 1 to 100. ( Default 75 )\n";
	print "[--coverage (INT)] Minimum %coverage to keep the result, should be a number between 0 to 100. ( Default 50 )\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
}

sub printpCOG{
	print "[--AAsPath (PATH)] Amino acids of all strains as fasta file paths, ( Default './Results/Annotations/AAs' )\n";
	print "[--strain_num (INT)] The total number of strains used for analysis, not including the reference genome\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
	print "[--evalue (FLOAT)] Maximum e-value to report alignments, ( Default 1e-3 )\n";
	print "[--id (INT)] Minimum identity% to report an alignment, ( Default 40 )\n";
	print "[--query_cover (INT)] Minimum query cover% to report an alignment, ( Default 70 )\n";
	print "[--subject_cover (INT)] Minimum subject cover% to report an alignment, ( Default 50 )\n";
}

sub printSTREE{
	print "[--seqfile (STRING)] Path of the sequence file for analysis.\n";
	print "[--seqtype (STRING)] Type Of Sequence (p, d, c for Protein, DNA, Codons, respectively). ( Default p )\n";
	print "[--threads (INT)] Number of threads to be used ( Default 4 )\n";
	print "[--bsnum (INT)] Replicates for bootstrap of IQ-TREE ( Default 500 )\n\n";
	print "[--fastboot (INT)] Replicates for ultrafast bootstrap of IQ-TREE ( must >=1000, Default 1000 )\n\n";
}

sub printACC{
	print "Applets in ACC include 'Assess' now\n";
	print "Parameters for Assess include the following:\n    [--scafPath (PATH)] Path for contigs/scaffolds ( Default 'Results/Assembles/Scaf/Illumina' )\n    [--Scaf_suffix (STRING)] The suffix of scaffolds or genome files." . RED, " User specified required",RESET . " ( Default .filtered.fas )\n    [--filter_length (INT)] Sequences shorter than the 'filter_length' will be deleted from the assembled genomes. ( Default 200 )\n\n";
}

sub printExamples{
	print "\nThe main usage is as follows, visit the official website for step by step examples: https://liaochenlanruo.github.io/pgcgap/\n\n";

	print ON_BLUE, "Example 1: Perform all functions for pair-end reads. For the sake of flexibility, the 'VAR' module needs to be added separately.", RESET . "\n\n";

	print YELLOW, "         pgcgap",RESET . MAGENTA, " --All ",RESET . RED, "--platform",RESET . " illumina " . RED, "--ReadsPath",RESET . " <PATH>" . RED, " --reads1",RESET . " <reads1 suffix>" . RED, " --reads2",RESET . " <reads2 suffix>" . RED, " --suffix_len",RESET . " <INT>" . RED, " --kmmer",RESET . " <INT> " . RED, "--PanTree",RESET . RED, " --genus",RESET . " <STRING>" . RED, " --species",RESET . " <STRING>" . RED, " --codon",RESET . " <INT>" . RED, " --strain_num",RESET . " <INT>" . RED " --threads",RESET . " <INT>" . MAGENTA, " --VAR",RESET . RED, " --refgbk",RESET . " <full path>" . RED, " --qualtype",RESET . " <STRING>" . "\n\n";

	print ON_BLUE, "Example 2: Conduct pair-end reads assembly.", RESET . "\n\n";

	print GREEN,"         # Assemble with AbySS:", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--Assemble",RESET . " " . RED, "--platform",RESET . " illumina " . RED, "--assembler",RESET . " abyss " . RED, "--ReadsPath",RESET . " <PATH> " . RED, "--reads1",RESET . " <reads1 suffix> " . RED, "--reads2",RESET . " <reads2 suffix> " . RED, "--suffix_len",RESET . " <INT> " . RED, "--kmmer",RESET . " <INT> " . RED, "--threads",RESET . " <INT>\n";

	print GREEN,"         # Assemble with SPAdes:", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--Assemble",RESET . " " . RED, "--platform",RESET . " illumina " . RED, "--assembler",RESET . " spades " . RED, "--ReadsPath",RESET . " <PATH> " . RED, "--reads1",RESET . " <reads1 suffix> " . RED, "--reads2",RESET . " <reads2 suffix> " . RED, "--suffix_len",RESET . " <INT> " . RED, "--threads",RESET . " <INT>\n";

	print GREEN,"         # Assemble with AUTO mode:", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--Assemble",RESET . " " . RED, "--platform",RESET . " illumina " . RED, "--assembler",RESET . " auto " . RED, "--ReadsPath",RESET . " <PATH> " . RED, "--reads1",RESET . " <reads1 suffix> " . RED, "--reads2",RESET . " <reads2 suffix> " . RED, "--suffix_len",RESET . " <INT> " . RED, "--kmmer",RESET . " <INT> " . RED, "--threads",RESET . " <INT>\n\n";

	print ON_BLUE, "Example 3: Conduct PacBio/Oxford reads assembly.", RESET . "\n\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--Assemble",RESET . " " . RED, "--platform",RESET . " [pacbio|oxford] " . RED, "--ReadsPath",RESET . " <PATH> " . RED, "--reads1",RESET . " <reads suffix> " . RED, "--suffix_len",RESET . " <INT> " . RED, "--genomeSize",RESET . " <STRING> " . RED, "--threads",RESET . " <INT>\n\n";

	print ON_BLUE, "Example 4: Conduct hybrid assembly.", RESET . "\n\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--Assemble",RESET . " " . RED, "--platform",RESET . " hybrid " . RED, "--ReadsPath",RESET . " <PATH> " . RED, "--short1",RESET . " <pair-end-reads1> " . RED, "--short2",RESET . " <pair-end-reads2> " . RED, "--long",RESET . " <long-reads> " . RED, "--hout",RESET . " <output_dir> " . RED, "--threads",RESET . " <INT>\n\n";

	print ON_BLUE, "Example 5: Conduct gene prediction and annotation.", RESET . "\n\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--Annotate",RESET . " " . RED, "--scafPath",RESET . " <PATH> " . RED, "--Scaf_suffix",RESET . " <STRING> " . RED, "--genus",RESET . " <STRING> " . RED, "--species",RESET . " <STRING> " . RED, "--codon",RESET . " <INT> " . RED, "--threads",RESET . " <INT>\n\n";

	print ON_BLUE, "Example 6: Constructing the phylogenetic trees of single-copy core proteins and core SNPs.", RESET . "\n\n";

	print GREEN,"         # Construct phylogenetic tree with FastTree (Quick without best fit model testing)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--CoreTree",RESET . " " . RED, "--CDsPath",RESET . " <PATH> " . RED, "--AAsPath",RESET . " <PATH> " . RED, "--codon",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--threads",RESET . " <INT> " . RED, "--fasttree\n\n",RESET;

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap, DEFAULT)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--CoreTree",RESET . " " . RED, "--CDsPath",RESET . " <PATH> " . RED, "--AAsPath",RESET . " <PATH> " . RED, "--codon",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--threads",RESET . " <INT> " . RED, "--bsnum",RESET . " <INT>\n\n";

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--CoreTree",RESET . " " . RED, "--CDsPath",RESET . " <PATH> " . RED, "--AAsPath",RESET . " <PATH> " . RED, "--codon",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--threads",RESET . " <INT> " . RED, "--fastboot",RESET . " <INT>\n\n";

	print ON_BLUE, "Example 7: Constructing a single-copy core protein tree only.", RESET . "\n\n";

	print GREEN,"         # Construct phylogenetic tree with FastTree (Quick without best fit model testing)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--CoreTree",RESET . " " . RED, "--CDsPath",RESET . " NO " . RED, "--AAsPath",RESET . " <PATH> " . RED, "--codon",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--threads",RESET . " <INT> " . RED, "--fasttree\n\n",RESET;

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap, DEFAULT)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--CoreTree",RESET . " " . RED, "--CDsPath",RESET . " NO " . RED, "--AAsPath",RESET . " <PATH> " . RED, "--codon",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--threads",RESET . " <INT> " . RED, "--bsnum",RESET . " <INT>\n\n";

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--CoreTree",RESET . " " . RED, "--CDsPath",RESET . " NO " . RED, "--AAsPath",RESET . " <PATH> " . RED, "--codon",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--threads",RESET . " <INT> " . RED, "--fastboot",RESET . " <INT>\n\n";

	print ON_BLUE, "Example 8: Conduct pan-genome analysis and construct a phylogenetic tree of single-copy core proteins called by roary.", RESET . "\n\n";

	print GREEN,"         # Construct phylogenetic tree with FastTree (Quick without best fit model testing)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--Pan",RESET . " " . RED, "--codon",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--threads",RESET . " <INT> " . RED, "--identi",RESET . " <INT> " . RED, "--GffPath",RESET . " <PATH> " . RED, "--PanTree",RESET . " " . RED, "--fasttree\n\n",RESET;

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap, DEFAULT)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--Pan",RESET . " " . RED, "--codon",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--threads",RESET . " <INT> " . RED, "--identi",RESET . " <INT> " . RED, "--GffPath",RESET . " <PATH> " . RED, "--PanTree",RESET . " " . RED, "--bsnum",RESET . " <INT>\n\n";

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--Pan",RESET . " " . RED, "--codon",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--threads",RESET . " <INT> " . RED, "--identi",RESET . " <INT> " . RED, "--GffPath",RESET . " <PATH> " . RED, "--PanTree",RESET . " " . RED, "--fastboot",RESET . " <INT>\n\n";

	print ON_BLUE, "Example 9: Inference of orthologous gene groups and construct a phylogenetic tree of single-copy Orthologue proteins.", RESET . "\n\n";

	print GREEN,"         # Construct phylogenetic tree with FastTree (Quick without best fit model testing)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--OrthoF",RESET . " " . RED, "--threads",RESET . " <INT> " . RED, "--AAsPath",RESET . " <PATH> " . RED, "--fasttree\n\n",RESET;

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap, DEFAULT)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--OrthoF",RESET . " " . RED, "--threads",RESET . " <INT> " . RED, "--AAsPath",RESET . " <PATH> " . RED, "--bsnum",RESET . " <INT>\n\n";

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)", RESET . "\n";

	print "         pgcgap ",RESET . MAGENTA, "--OrthoF",RESET . " " . RED, "--threads",RESET . " <INT> " . RED, "--AAsPath",RESET . " <PATH> " . RED, "--fastboot",RESET . " <INT>\n\n";

	print ON_BLUE, "Example 10: Compute whole-genome Average Nucleotide Identity (ANI).", RESET . "\n\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--ANI",RESET . " " . RED, "--threads",RESET . " <INT> " . RED, "--queryL",RESET . " <FILE> " . RED, "--refL",RESET . " <FILE> " . RED, "--Scaf_suffix",RESET . " <STRING>\n\n";

	print ON_BLUE, "Example 11: Genome and metagenome similarity estimation using MinHash.", RESET . "\n\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--MASH",RESET . " " . RED, "--scafPath",RESET . " <PATH> " . RED, "--Scaf_suffix",RESET . " <STRING>\n\n";

	print ON_BLUE, "Example 12: Run COG annotation for each strain.", RESET . "\n\n";

	print YELLOW, "          pgcgap ",RESET . MAGENTA, "--pCOG",RESET . " " . RED, "--strain_num",RESET . " <INT> " . RED, "--threads",RESET . " <INT> " . RED, "--evalue",RESET . " <FLOAT> " . RED, "--id",RESET . " <INT> " . RED, "--query_cover",RESET . " <INT> " .RED, "--subject_cover",RESET . " <INT> " . RED, "--AAsPath",RESET . " <PATH>\n\n";

	print ON_BLUE, "Example 13: Variants calling and phylogenetic tree construction based on a reference genome.", RESET . "\n\n";

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap, DEFAULT)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--VAR",RESET . " " . RED, "--threads",RESET . " <INT> " . RED, "--refgbk",RESET . " <FILE with full path> " . RED, "--ReadsPath",RESET . " <PATH> " . RED, "--reads1",RESET . " <STRING> " . RED, "--reads2",RESET . " <STRING> " . RED, "--suffix_len",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--qualtype",RESET . " <STRING> " . RED, "--bsnum",RESET . " <INT>\n\n";

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--VAR",RESET . " " . RED, "--threads",RESET . " <INT> " . RED, "--refgbk",RESET . " <FILE with full path> " . RED, "--ReadsPath",RESET . " <PATH> " . RED, "--reads1",RESET . " <STRING> " . RED, "--reads2",RESET . " <STRING> " . RED, "--suffix_len",RESET . " <INT> " . RED, "--strain_num",RESET . " <INT> " . RED, "--qualtype",RESET . " <STRING> " . RED, "--fastboot",RESET . " <INT>\n\n";

	print ON_BLUE, "Example 14: Screening of contigs for antimicrobial and virulence genes.", RESET . "\n\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--AntiRes",RESET . " " . RED, "--scafPath",RESET . " <PATH> " . RED, "--Scaf_suffix",RESET . " <STRING> " . RED, "--threads",RESET . " <INT> " . RED, "--db",RESET . " <STRING> " . RED, "--identity",RESET . " <INT> " . RED, "--coverage",RESET . " <INT>\n";

	print ON_BLUE, "Example 15: Construct a phylogenetic tree based on multiple sequences in one file.", RESET . "\n\n";

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap, DEFAULT)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--STREE",RESET . " " . RED, "--seqfile",RESET . " <PATH> " . RED, "--seqtype",RESET . " <p|d|c> " . RED, "--bsnum",RESET . " <INT> " . RED, "--threads",RESET . " <INT>\n\n";

	print GREEN,"         # Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)", RESET . "\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--STREE",RESET . " " . RED, "--seqfile",RESET . " <PATH> " . RED, "--seqtype",RESET . " <p|d|c> " . RED, "--fastboot",RESET . " <INT> " . RED, "--threads",RESET . " <INT>\n\n";

	print ON_BLUE, "Example 16: Perform the short sequences filter from the assembled genome and get the genome status.", RESET . "\n\n";

	print YELLOW, "         pgcgap ",RESET . MAGENTA, "--ACC",RESET . " " . RED, "--Assess",RESET . " " . RED, "--scafPath",RESET . " <PATH> " . RED, "--Scaf_suffix",RESET . " <STRING> " . RED, "--filter_length",RESET . " <INT>\n\n";
}

if ( grep {$_ eq "Assemble"} @ARGV ){
	printAssemble();
	exit 0;
}

if ( grep {$_ eq "Annotate"} @ARGV ){
	printAnnotate();
	exit 0;
}

if ( grep {$_ eq "CoreTree"} @ARGV ){
	printCoreTree();
	exit 0;
}

if ( grep {$_ eq "Pan"} @ARGV ){
	printPan();
	exit 0;
}

if ( grep {$_ eq "OrthoF"} @ARGV ){
	printOrthoF();
	exit 0;
}

if ( grep {$_ eq "ANI"} @ARGV ){
	printANI();
	exit 0;
}

if ( grep {$_ eq "MASH"} @ARGV ){
	printMASH();
	exit 0;
}

if ( grep {$_ eq "VAR"} @ARGV ){
	printVAR();
	exit 0;
}

if ( grep {$_ eq "pCOG"} @ARGV ){
	printpCOG();
	exit 0;
}

if ( grep {$_ eq "AntiRes"} @ARGV ){
	printAntiRes();
	exit 0;
}

if ( grep {$_ eq "Examples"} @ARGV ){
	printExamples();
	exit 0;
}

if ( grep {$_ eq "STREE"} @ARGV ) {
	printSTREE();
	exit 0;
}

if ( grep {$_ eq "ACC"} @ARGV ) {
	printACC();
	exit 0;
}
