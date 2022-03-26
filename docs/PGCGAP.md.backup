# PGCGAP - the <font color=red>P</font>rokaryotic <font color=red>G</font>enomics and <font color=red>C</font>omparative <font color=red>G</font>enomics <font color=red>A</font>nalysis <font color=red>P</font>ipeline
---
<center>
![Platform](https://img.shields.io/badge/Platform-WSL%2FLinux%2FmacOS-green)
[![GitHub license](https://img.shields.io/github/license/liaochenlanruo/pgcgap)](https://github.com/liaochenlanruo/pgcgap/blob/master/LICENSE)
[![GitHubversion](https://badge.fury.io/gh/liaochenlanruo%2Fpgcgap.svg)](https://badge.fury.io/gh/liaochenlanruo%2Fpgcgap)
![Downloads conda](https://img.shields.io/conda/dn/bioconda/pgcgap.svg?style=flat)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pgcgap/README.html)
</center>

<center>
<script type="text/javascript" src="//rf.revolvermaps.com/0/0/4.js?i=0ett3z77b0x&amp;m=0&amp;h=128&amp;c=ff0000&amp;r=0" async="async"></script>
&nbsp;&nbsp;&nbsp;&nbsp;
<script type="text/javascript" src="//rf.revolvermaps.com/0/0/0.js?i=0ett3z77b0x&amp;d=3&amp;p=1&amp;b=1&amp;w=293&amp;g=2&amp;f=arial&amp;fs=13&amp;r=0&amp;c0=ff8a00&amp;c1=0006ff&amp;c2=000000&amp;ic0=0&amp;ic1=0" async="async"></script>
</center>

<center>
[English Readme](https://liaochenlanruo.github.io/pgcgap) | [Chinese Readme](https://liaochenlanruo.github.io/2019/04/28/PGCGAP%E4%B8%AD%E6%96%87%E8%AF%B4%E6%98%8E/)


		  ____       ____      ____     ____       _        ____    
		U|  _"\ u U /"___|u U /"___| U /"___|u U  /"\  u  U|  _"\ u 
		\| |_) |/ \| |  _ / \| | u   \| |  _ /  \/ _ \/   \| |_) |/ 
		 |  __/    | |_| |   | |/__   | |_| |   / ___ \    |  __/   
		 |_|        \____|    \____|   \____|  /_/   \_\   |_|      
		 ||>>_      _)(|_    _// \\    _)(|_    \\    >>   ||>>_    
		(__)__)    (__)__)  (__)(__)  (__)__)  (__)  (__) (__)__)   

</center>

## Introduction

PGCGAP is a pipeline for prokaryotic comparative genomics analysis. It can take the pair-end reads, Oxford reads or PacBio reads as input. In addition to genome assembly, gene prediction and annotation, it can also get common comparative genomics analysis results such as phylogenetic trees of single-core proteins and core SNPs, pan-genome, whole-genome Average Nucleotide Identity (ANI), orthogroups and orthologs, COG annotations, substitutions (SNPs) and insertions/deletions (indels), and antimicrobial and virulence genes mining with only one line of commands. __To follow this document, please upgrade PGCGAP to version v1.0.31 or later.__

## Installation

The software was tested successfully on Windows WSL, Linux x64 platform, and macOS. <font color=#ff0000>Because this software relies on a large number of other software, so it is recommended to install with __[Bioconda](https://bioconda.github.io/index.html)__</font>.






__Step1: Install PGCGAP__

<pre>
$conda create -n pgcgap python=3
$conda activate pgcgap
$conda install pgcgap (Users in China can input "conda install -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda pgcgap" for instead)
</pre>

<p><strong>Notice: What should we do when the installation is slow? </strong>As more and more software is contained in CONDA and the index files become larger, the search space for the software that satisfies all the software dependencies in the environment becomes larger and larger when installing a new software, making "Solving Environment" slower and slower. Sometimes we can't even install the software through CONDA. In fact, we can do something instead of just waiting.<br/></p>


- Method 1: use mamba to deal with the slow development of "solving environment" when using CONDA.
    <pre>
    $conda activate pgcgap
    $conda install mamba -c conda-forge
    $mamba install pgcgap
    </pre>

- Method 2: use "environment.yaml" we provided to deal with the slow development of "solving environment" when using CONDA. Run the following commands to download the <a href="https://github.com/liaochenlanruo/pgcgap/blob/master/conda/pgcgap_latest_env.yml">latest environmental file</a> and install PGCGAP:
    <pre>
    # download pgcgap_latest_env.yml
    $wget https://github.com/liaochenlanruo/pgcgap/blob/master/conda/pgcgap_latest_env.yml
	
    # create a conda environment named as pgcgap and install the latest version of PGCGAP
    $conda env create -f pgcgap_latest_env.yml
    </pre>

__Step2: Setup COG database__ (Users should execute this after the first installation of pgcgap)

<pre>
$conda activate pgcgap
$pgcgap --setup-COGdb
$conda deactivate
</pre>

Users with [docker container](https://hub.docker.com/) installed have another choice to install PGCGAP.

<pre>
$docker pull quay.io/biocontainers/pgcgap:<tag>
</pre>

(see [pgcgap/tags](https://quay.io/repository/biocontainers/pgcgap?tab=tags) for valid values for &lt;tag&gt;)


## Required dependencies


- [Abricate](https://github.com/tseemann/abricate)
- [ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss/)
- [Canu](http://canu.readthedocs.org/)
- [CD-HIT](http://weizhongli-lab.org/cd-hit/)
- [Coreutils](https://www.gnu.org/software/coreutils/)
- [Diamond](https://github.com/bbuchfink/diamond)
- [FastANI](https://github.com/ParBLiSS/FastANI)
- [Fastme](http://www.atgc-montpellier.fr/fastme/binaries.php)
- [Fastp](https://github.com/OpenGene/fastp)
- [FastTree](http://www.microbesonline.org/fasttree)
- [Gubbins](https://github.com/sanger-pathogens/gubbins) >=2.3.4
- [Htslib](https://github.com/samtools/htslib)
- [IQ-TREE](http://www.iqtree.org/)
- [Mafft](https://mafft.cbrc.jp/alignment/software/)
- [Mash](https://github.com/marbl/Mash)
- [Mmseqs2](https://github.com/soedinglab/mmseqs2)
- [Muscle](https://www.ebi.ac.uk/Tools/msa/muscle/)
- [NCBI-blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [OrthoFinder](https://github.com/davidemms/OrthoFinder)
- [OpenJDK8](https://openjdk.java.net/)
- [PAL2NAL v14](http://www.bork.embl.de/pal2nal/)
- [trimAL](http://trimal.cgenomics.org)
- [Perl](http://www.perl.org/get.html) & the modules
  - [perl-bioperl](http://metacpan.org/pod/BioPerl)
  - [perl-data-dumper](http://metacpan.org/pod/Data::Dumper)
  - [perl-file-tee](http://metacpan.org/pod/File::Tee)
  - [perl-getopt-long](http://metacpan.org/pod/Getopt::Long)
  - [perl-pod-usage](http://search.cpan.org/~marekr/Pod-Usage-1.69/)
  - [perl-parallel-forkmanager](https://metacpan.org/pod/release/DLUX/Parallel-ForkManager-0.7.5/ForkManager.pm)
- [Prokka](https://github.com/tseemann/prokka)
- [Python](https://www.python.org/) & the modules
  - [biopython](https://biopython.org/)
  - [matplotlib](https://matplotlib.org/)
  - [numpy](http://www.numpy.org/)
  - [pandas](http://pandas.pydata.org/)
  - [seaborn](http://seaborn.pydata.org/)
- [R](https://www.r-project.org/) & the packages
  - [corrplot](https://cran.r-project.org/web/packages/corrplot/index.html)
  - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/)
  - [gplots](https://cran.r-project.org/web/packages/gplots/)
  - [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html)
  - [plotrix](https://cran.r-project.org/web/packages/plotrix/)
- [Roary](https://sanger-pathogens.github.io/Roary/)
- [Sickle-trim](https://github.com/najoshi/sickle)
- [Snippy](https://github.com/tseemann/snippy)
- [Snp-sites](https://github.com/sanger-pathogens/snp-sites)
- [unicycler](https://github.com/rrwick/Unicycler)
- [wget](https://www.gnu.org/software/wget/)

## Usage


- __Print the help messages:__
    <pre>$pgcgap --help</pre>
- __Check for update:__
    <pre>$pgcgap --check-update</pre>
- __General usage:__
    <pre>$pgcgap [modules] [options]</pre>

- __Show parameters for each module:__
    <pre>$pgcgap [Assemble|Annotate|ANI|AntiRes|CoreTree|MASH|OrthoF|Pan|pCOG|VAR|STREE|ACC]</pre>

- __Show examples of each module:__
    <pre>$pgcgap Examples</pre>

<br/>

- __Setup COG database:__ (Users should execute this after the first installation of pgcgap)
    <pre>$pgcgap --setup-COGdb</pre>
<br/>

- __Modules:__

  - __[\-\-All]__                          Perform Assemble, Annotate, CoreTree, Pan, OrthoF, ANI, MASH, AntiRes and pCOG functions with one command

  - __[\-\-Assemble]__                     Assemble reads (short, long or hybrid) into contigs

  - __[\-\-Annotate]__                     Genome annotation

  - __[\-\-CoreTree]__                     Construct single-core proteins tree and SNPs tree of single-copy core genes

  - __[\-\-Pan]__                          Run "roary" pan-genome pipeline with gff3 files, and construct a phylogenetic tree with the sing-copy core proteins called by roary

  - __[\-\-OrthoF]__                       Identify orthologous protein sequence
                                         families with "OrthoFinder"

  - __[\-\-ANI]__                          Compute whole-genome Average Nucleotide Identity ( ANI )

  - __[\-\-MASH]__                         Genome and metagenome similarity estimation using MinHash

  - __[\-\-pCOG]__                         Run COG annotation for each strain (*.faa),
				          and generate a table containing the relative
				          abundance of each flag for all strains

  - **[\-\-VAR]**                          Rapid haploid variant calling and core
                                         genome alignment with "Snippy"

  - __[\-\-AntiRes]__                      Screening of contigs for antimicrobial and virulence genes

  - __[\-\-STREE]__                        Construct a phylogenetic tree based on multiple sequences in one file

  - __[\-\-ACC]__                          Other useful gadgets (now includes 'Assess' for filtering short sequences in the genome and assessing the status of the genome only)

<br/>

- __Global Options:__

  - **[\-\-strain\_num (INT)]**             [Required by "\-\-All", "\-\-CoreTree", "\-\-Pan", "\-\-VAR" and "\-\-pCOG"]
                                         The total number of strains used for analysis, not including the reference genome

  - __[\-\-ReadsPath (PATH)]__             [Required by "\-\-All", "\-\-Assemble" and "\-\-VAR"]
                                         Reads of all strains as file paths ( Default ./Reads/Illumina )

  - __[\-\-scafPath (PATH)]__              [Required by "\-\-All", "\-\-Assess", "\-\-Annotate", "\-\-MASH" and "\-\-AntiRes"] Path for contigs/scaffolds (Default "Results/Assembles/Scaf/Illumina")

  - __[\-\-AAsPath (PATH)]__               [Required by "\-\-All", "\-\-CoreTree", "\-\-Pan", "\-\-OrthoF" and "\-\-pCOG"] Amino acids of all strains as fasta file paths,
                                         ( Default "./Results/Annotations/AAs" )

  - __[\-\-reads1 (STRING)]__              [Required by "\-\-All", "\-\-Assemble" and
                                         "\-\-VAR"] The suffix name of reads 1 ( for
                                         example: if the name of reads 1 is
                                         "YBT-1520\_L1\_I050.R1.clean.fastq.gz",
                                         "YBT-1520" is the strain same, so the
                                         suffix name should be ".R1.clean.fastq.gz")

  - __[\-\-reads2 (STRING)]__              [Required by "\-\-All", "\-\-Assemble" and
                                         "\-\-VAR"] The suffix name of reads 2( for
                                         example: if the name of reads 2 is
                                         "YBT-1520\_2.fq", the suffix name should be "\_2.fq" )

  - **[\-\-Scaf\_suffix (STRING)]**         [Required by "\-\-All", "\-\-Assess", "\-\-Annotate" "\-\-MASH", "\-\-ANI" and "\-\-AntiRes"] The suffix of scaffolds or genome files. This is an important parameter that must be set (Default -8.fa)

  - **[\-\-filter\_length (INT)]**          [Required by "\-\-All", "\-\-Assemble" and "\-\-Assess"]> Sequences shorter than the 'filter\_length' will be deleted from the assembled genomes. ( Default 200 )

  - __[\-\-codon (INT)]__                  [Required by "\-\-All", "\-\-Annotate", "\-\-CoreTree" and "\-\-Pan"] Translation table ( Default 11 )

      -    1                             Universal code
      -    2                             Vertebrate mitochondrial code
      -    3                             Yeast mitochondrial code
      -    4                             Mold, Protozoan, and Coelenterate Mitochondrial code and Mycoplasma/Spiroplasma code
      -    5                             Invertebrate mitochondrial
      -    6                             Ciliate, Dasycladacean and Hexamita nuclear code
      -    9                             Echinoderm and Flatworm mitochondrial code
      -    10                            Euplotid nuclear code
      -    11                            Bacterial, archaeal and plant plastid code ( Default )
      -    12                            Alternative yeast nuclear code
      -    13                            Ascidian mitochondrial code
      -    14                            Alternative flatworm mitochondrial code
      -    15                            Blepharisma nuclear code
      -    16                            Chlorophycean mitochondrial code
      -    21                            Trematode mitochondrial code
      -    22                            Scenedesmus obliquus mitochondrial code
      -    23                            Thraustochytrium mitochondrial code

  - __[\-\-suffix\_len (INT)]__             [Required by "\-\-All", "\-\-Assemble" and
                                         "\-\-VAR"] __(Strongly recommended)__ The suffix length of the reads,
                                         that is the length of your reads name
                                         minus the length of your strain name. For
                                         example the \-\-suffix\_len of
                                         "YBT-1520\_L1\_I050.R1.clean.fastq.gz" is 26
                                         ( "YBT-1520" is the strain name ) ( Default 0 )
  - __[--fasttree]__             [Can be used with "CoreTree", "Pan" and "OrthoF"] Use FastTree to construct phylogenetic tree quickly instead of IQ-TREE
  - __[--bsnum (INT)]__             [Required by "CoreTree", "Pan", "OrthoF", "STREE", and "VAR"] Replicates for bootstrap of IQ-TREE ( Default 500 )
  - __[--fastboot (INT)]__             [Required by "CoreTree", "Pan", "OrthoF", "STREE", and "VAR"] Replicates for ultrafast bootstrap of IQ-TREE. ( must >= 1000, Default 1000 )

  - __[\-\-logs (STRING)]__                Name of the log file ( Default Logs.txt )

  - __[\-\-threads (INT)]__                Number of threads to be used ( Default 4 )
<br/>

- __Local Options:__

  - __\-\-Assemble__

      - __[\-\-platform (STRING)]__         [Required] Sequencing Platform, "illumina", "pacbio", "oxford" and "hybrid" available ( Default illumina )

      - __[\-\-assembler (STRING)]__        [Required] Software used for illumina reads assembly, "abyss" and "spades" available ( Default abyss )

      - __[\-\-kmmer (INT)]__               [Required] k-mer size for genome assembly of Illumina data ( Default 81 )

      - __[\-\-genomeSize (STRING)]__       [Required] An estimate of the size of the genome. Common suffixes are allowed, for example, 3.7m or 2.8g. Needed by PacBio data and Oxford data ( Default Unset )

      - __[\-\-short1 (STRING)]__           [Required] FASTQ file of first short reads in each pair. Needed by hybrid assembly ( Default Unset )

      - __[\-\-short2 (STRING)]__           [Required] FASTQ file of second short reads in each pair. Needed by hybrid assembly ( Default Unset )

      - __[\-\-long (STRING)]__             [Required] FASTQ or FASTA file of long reads. Needed by hybrid assembly ( Default Unset )

      - __[\-\-hout (STRING)]__             [Required] Output directory for hybrid assembly ( Default ../../Results/Assembles/Hybrid )

  - __\-\-Annotate__

      - __[\-\-genus (STRING)]__           Genus name of your strain ( Default "NA" )

      - __[\-\-species (STRING)]__         Species name of your strain ( Default "NA")
  <br/>

  - **\-\-CoreTree**

      - __[\-\-CDsPath (PATH)]__           [Required] CDs of all strains as fasta file
                                         paths ( Default "./Results/Annotations/CDs" ),
					 if set to "NO", the SNPs of single-copy core
					 genes will not be called

      - __[-c (FLOAT)]__                 Sequence identity threshold, ( Default 0.5)

      - __[-n (INT)]__                   Word_length,  -n 2 for thresholds 0.4-0.5,
                                         -n 3 for thresholds 0.5-0.6, -n 4 for thresholds 0.6-0.7,
                                         -n 5 for thresholds 0.7-1.0 ( Default 2 )

      - __[-G (INT)]__                   Use global (set to 1) or local (set to 0)
                                         sequence identity, ( Default 0 )

      - __[-t (INT)]__                   Tolerance for redundance ( Default 0 )

      - __[-aL (FLOAT)]__                Alignment coverage for the longer
                                         sequence. If set to 0.9, the alignment
                                         must cover 90% of the sequence ( Default 0.5 )

      - __[-aS (FLOAT)]__                Alignment coverage for the shorter sequence.
                                         If set to 0.9, the alignment must covers
                                         90% of the sequence ( Default 0.7 )

      - __[-g (INT)]__                   If set to 0, a sequence is clustered to
                                         the first cluster that meets the threshold
                                         (fast cluster). If set to 1, the program
                                         will cluster it into the most similar
                                         cluster that meets the threshold (accurate
                                         but slow mode, Default 1)

      - __[-d (INT)]__                   length of description in .clstr file. if
                                         set to 0, it takes the fasta defline and
                                         stops at first space ( Default 0 )
  <br/>

  - __\-\-Pan__

      - __[\-\-identi (INT)]__                  Minimum percentage identity for blastp ( Default 95 )
  <br/>
      - __[\-\-PanTree]__                        Construct a phylogenetic tree of single-copy core proteins called by roary
 <br/>
      - __[\-\-GffPath (PATH)]__           [Required] Gff files of all strains as paths
                                         ( Default "./Results/Annotations/GFF" )
  <br/>

  - __\-\-OrthoF__

      - __[\-\-Sprogram (STRING)]__        Sequence search program, Options: blast,
                                         mmseqs, blast_gz, diamond ( Default diamond)

      - __[\-\-PanTree]__                  Construct a phylogenetic tree of single-copy core proteins called by roary
  <br/>

  - __\-\-ANI__

      - __[\-\-queryL (FILE)]__            [Required] The file containing paths to query
                                         genomes, one per line ( Default scaf.list )

      - __[\-\-refL (FILE)]__              [Required] The file containing paths to reference genomes,
                                         one per line. ( Default scaf.list )

      - __[\-\-ANIO (FILE)]__              The name of the output file ( Default "Results/ANI/ANIs" )
  <br/>

  - **\-\-VAR**

      - __[\-\-refgbk (FILE)]__            [Required] The full path and name of
                                         reference genome in GENBANK format (
                                         recommended ), fasta format is also OK.
                                         For example: "/mnt/g/test/ref.gbk"

      - __[\-\-qualtype (STRING)]__        [Required] Type of quality values (solexa
                                         (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7),
				         sanger (which is CASAVA >= 1.8)). ( Default sanger )

      - __[\-\-qual (INT)]__               Threshold for trimming based on average
                                         quality in a window. ( Default 20 )

      - __[\-\-length (INT)]__             Threshold to keep a read based on length
                                         after trimming. ( Default 20 )

      - __[\-\-mincov (INT)]__             The minimum number of reads covering a
                                         site to be considered ( Default 10 )

      - __[\-\-minfrac (FLOAT)]__          The minimum proportion of those reads
                                         which must differ from the reference ( Default 0.9 )

      - __[\-\-minqual (INT)]__            The minimum VCF variant call "quality" ( Default 100 )

      - __[\-\-ram (INT)]__                Try and keep RAM under this many GB ( Default 8 )

      - __[\-\-tree\_builder (STRING)]__    Application to use for tree building
                                         [raxml|fasttree|hybrid] ( Default fasttree)

      - __[\-\-iterations (INT)]__         Maximum No. of iterations for gubbins ( Default 5 )
<br/>

  - __\-\-AntiRes__

      - __[\-\-db (STRING)]__              [Required] The database to use, options: all, [argannot](https://www.ncbi.nlm.nih.gov/pubmed/24145532), 
	                                 [card](https://www.ncbi.nlm.nih.gov/pubmed/27789705), [ecoh](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5343136/), ecoli_vf, [megares](https://megares.meglab.org), [ncbi](https://www.biorxiv.org/content/10.1101/550707v1), [plasmidfinder](https://www.ncbi.nlm.nih.gov/pubmed/24777092), [resfinder](https://www.ncbi.nlm.nih.gov/pubmed/22782487) and [vfdb](https://www.ncbi.nlm.nih.gov/pubmed/26578559). ( Default all)

      - __[\-\-identity (INT)]__           [Required] Minimum %identity to keep the result, 
	                                 should be a number between 1 to 100. ( Default 75 )

      - __[\-\-coverage (INT)]__           [Required] Minimum %coverage to keep the result, 
	                                 should be a number between 0 to 100. ( Default 50 )

  - __\-\-STREE__

      - __[\-\-seqfile (STRING)]__        [Required] Path of the sequence file for analysis.

      - __[\-\-seqtype (INT)]__           [Required] Type Of Sequence (p, d, c for Protein, DNA, Codons, respectively). ( Default p )

      - __[\-\-bsnum (INT)]__             [Required] Times for bootstrap. ( Default 1000 )

  - __\-\-ACC__

      - __[\-\-Assess (STRING)]__           Filter short sequences in the genome and assess the status of the genome

- __Paths of external programs__

    Not needed if they were in the environment variables path. Users can check with the "\-\-check-external-programs" option for the essential programs.
  <br/>

  - __[\-\-abricate-bin (PATH)]__          Path to abyss binary file. Default tries if abyss is in PATH;

  - __[\-\-abyss-bin (PATH)]__             Path to abyss binary file. Default tries if abyss is in PATH;

  - __[\-\-canu-bin (PATH)]__              Path to canu binary file. Default tries if canu is in PATH;

  - __[\-\-cd-hit-bin (PATH)]__            Path to cd-hit binary file. Default tries if cd-hit is in PATH;

  - __[\-\-fastANI-bin (PATH)]__           Path to the fastANI binary file. Default tries if fastANI is in PATH;

  - __[\-\-Gblocks-bin (PATH)]__           Path to the Gblocks binary file. Default tries if Gblocks is in PATH;

  - __[\-\-gubbins-bin (PATH)]__           Path to the run\_gubbins.py binary file. Default tries if run\_gubbins.py is in PATH;

  - __[\-\-iqtree-bin (PATH)]__            Path to the iqtree binary file. Default tries if iqtree is in PATH;

  - __[\-\-mafft-bin (PATH)]__             Path to mafft binary file. Default tries if mafft is in PATH;

  - __[\-\-mash-bin (PATH)]__              Path to the mash binary file. Default tries if mash is in PATH.

  - __[\-\-modeltest-ng-bin (PATH)]__      Path to the modeltest-ng binary file. Default tries if modeltest-ng is in PATH.

  - __[\-\-muscle-bin (PATH)]__            Path to the muscle binary file. Default tries if muscle is in PATH.

  - __[\-\-orthofinder-bin (PATH)]__       Path to the orthofinder binary file. Default tries if orthofinder is in PATH;

  - __[\-\-pal2nal-bin (PATH)]__           Path to the pal2nal.pl binary file. Default tries if pal2nal.pl is in PATH;

  - __[\-\-prodigal-bin (PATH)]__          Path to prodigal binary file. Default tries if prodigal is in PATH;

  - __[\-\-prokka-bin (PATH)]__            Path to prokka binary file. Default tries if prokka is in PATH;

  - __[\-\-raxml-ng-bin (PATH)]__          Path to the raxml-ng binary file. Default tries if raxml-ng is in PATH;

  - __[\-\-roary-bin (PATH)]__             Path to the roary binary file. Default tries if roary is in PATH;

  - __[\-\-sickle-bin (PATH)]__            Path to the sickle-trim binary file. Default tries if sickle is in PATH.

  - __[\-\-snippy-bin (PATH)]__            Path to the snippy binary file. Default tries if snippy is in PATH;

  - __[\-\-snp-sites-bin (PATH)]__         Path to the snp-sites binary file. Default tries if snp-sites is in PATH;

  - __[\-\-unicycler-bin (PATH)]__         Path to the unicycler binary file. Default tries if unicycler is in PATH;
<br/>

- __Setup COG database__

  - __[\-\-setup-COGdb]__                  Users should execute this after first installation of pgcgap.
<br/>

- Check the required external programs (__It is strongly recommended that this step be performed after the installation of PGCGAP__):
    <pre>$pgcgap --check-external-programs</pre>

<br/>

## Examples

  - __Example 1:__ Perform all functions, take the *Escherichia coli* as an example, total 6 strains for analysis.<br/>

    __Notice__: For the sake of flexibility, The "VAR" function needs to be added additionally.<br/>

    <pre>
    $pgcgap --All --platform illumina --filter_length 200 --ReadsPath Reads/Illumina --reads1 _1.fastq.gz --reads2 _2.fastq.gz --suffix_len 11 --kmmer 81 --genus Escherichia --species “Escherichia coli” --codon 11 --strain_num 6 --threads 4 --VAR --refgbk /mnt/h/PGCGAP_Examples/Reads/MG1655.gbff --qualtype sanger
    </pre>

  - __Example 2:__ Genome assembly.

    - Illumina reads assembly

        In this dataset, the naming format of the genome is “strain\_1.fastq.gz” and “strain\_2.fastq.gz”. The string after the strain name is “\_1.fastq.gz”, and its length is 11, so "\-\-suffix\_len" was set to 11.

         <pre>
     $pgcgap --Assemble --platform illumina --assembler abyss --filter_length 200 --ReadsPath Reads/Illumina --reads1 _1.fastq.gz --reads2 _2.fastq.gz --kmmer 81 --threads 4 --suffix_len 11
     $pgcgap --Assemble --platform illumina --assembler spades --filter_length 200 --ReadsPath Reads/Illumina --reads1 _1.fastq.gz --reads2 _2.fastq.gz --threads 4 --suffix_len 11
     $pgcgap --Assemble --platform illumina --assembler auto --filter_length 200 --ReadsPath Reads/Illumina --reads1 _1.fastq.gz --reads2 _2.fastq.gz --kmmer 81 --threads 4 --suffix_len 11</pre>

    - Oxford reads assembly

          Oxford nanopore only produces one reads file, so only the parameter of "\-\-reads1" needs to be set, where the value is ".fasta". “\-\-genomeSize” is the estimated genome size, and users can check the genome size of similar strains in the NCBI database for reference. The parameter was set to "4.8m" here. The suffix of the reads file here is ".fasta" and its length is 6, so "\-\-suffix_len" was set to 6.

         <pre>$pgcgap --Assemble --platform oxford --filter_length 200 --ReadsPath Reads/Oxford --reads1 .fasta --genomeSize 4.8m --threads 4 --suffix_len 6</pre>

    - PacBio reads assembly

         PacBio also produces only one reads file "pacbio.fastq", the parameter settings are similar to Oxford. The strain name is "pacbio" with the suffix ".fastq" and the suffix length is 6, so "\-\-suffix_len" was set to 6.

         <pre>$pgcgap --Assemble --platform pacbio --filter_length 200 --ReadsPath Reads/PacBio --reads1 .fastq --genomeSize 4.8m --threads 4 --suffix_len 6</pre>

    - Hybrid assembly of short reads and long reads

         Paired-end short reads and long reads in the directory “Reads/Hybrid/” were used as inputs. Illumina reads and long reads must be from the same isolates.

         <pre>$pgcgap --Assemble --platform hybrid --ReadsPath Reads/Hybrid --short1 short_reads_1.fastq.gz --short2 short_reads_2.fastq.gz --long long_reads_high_depth.fastq.gz --threads 4</pre>

  - __Example 3__: Gene prediction and annotation

     <pre>$pgcgap --Annotate --scafPath Results/Assembles/Scaf/Illumina --Scaf_suffix -8.fa --genus Escherichia --species “Escherichia coli” --codon 11 --threads 4</pre>

  - __Example 4__: Constructing single-copy core protein tree and core SNPs tree

     <pre># Construct phylogenetic tree with FastTree (Quick without best fit model testing)<br/>$pgcgap --CoreTree --CDsPath Results/Annotations/CDs --AAsPath Results/Annotations/AAs --codon 11 --strain_num 6 --threads 4 --fasttree<br/># Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap)<br/>$pgcgap --CoreTree --CDsPath Results/Annotations/CDs --AAsPath Results/Annotations/AAs --codon 11 --strain_num 6 --threads 4 --bsnum 500<br/># Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)<br/>$pgcgap --CoreTree --CDsPath Results/Annotations/CDs --AAsPath Results/Annotations/AAs --codon 11 --strain_num 6 --threads 4 --fastboot 1000</pre>

  - __Example 5:__ Constructing single-copy core protein tree only.
    <pre># Construct phylogenetic tree with FastTree (Quick without best fit model testing)<br/>$pgcgap --CoreTree --CDsPath NO --AAsPath Results/Annotations/AAs --codon 11 --strain_num 6 --threads 4 --fasttree<br/># Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap)<br/>$pgcgap --CoreTree --CDsPath NO --AAsPath Results/Annotations/AAs --codon 11 --strain_num 6 --threads 4 --bsnum 500<br/># Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)<br/>$pgcgap --CoreTree --CDsPath NO --AAsPath Results/Annotations/AAs --codon 11 --strain_num 6 --threads 4 --fastboot 1000</pre>

  - __Example 6:__ Conduct pan-genome analysis and construct a phylogenetic tree of single-copy core proteins called by roary. **Applicable to v1.0.27 and later**.

    <pre># Construct phylogenetic tree with FastTree (Quick without best fit model testing)<br/>$pgcgap --Pan --codon 11 --identi 95 --strain_num 6 --threads 4 --GffPath Results/Annotations/GFF --PanTree --fasttree<br/># Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap)<br/>$pgcgap --Pan --codon 11 --identi 95 --strain_num 6 --threads 4 --GffPath Results/Annotations/GFF --PanTree --bsnum 500<br/># Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)<br/>$pgcgap --Pan --codon 11 --identi 95 --strain_num 6 --threads 4 --GffPath Results/Annotations/GFF --PanTree --fastboot 1000</pre>

  - __Example 7:__ Inference of orthologous gene groups and construct a phylogenetic tree of single-copy Orthologue proteins. __Applicable to v1.0.29 and later__

    <pre># Construct phylogenetic tree with FastTree (Quick without best fit model testing)<br/>$pgcgap --OrthoF --threads 4 --AAsPath Results/Annotations/AAs --fasttree<br/># Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap)<br/>$pgcgap --OrthoF --threads 4 --AAsPath Results/Annotations/AAs --bsnum 500<br/># Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)<br/>$pgcgap --OrthoF --threads 4 --AAsPath Results/Annotations/AAs --fastboot 1000</pre>

  - __Example 8:__ Compute whole-genome Average Nucleotide Identity (ANI).

    <pre>$pgcgap --ANI --threads 4 --queryL scaf.list --refL scaf.list --ANIO Results/ANI/ANIs --Scaf_suffix .fa</pre>

  - __Example 9:__ Genome and metagenome similarity estimation using MinHash
    <pre>$pgcgap --MASH --scafPath <PATH> --Scaf_suffix <STRING></pre>

  - __Example 10:__ Run COG annotation for each strain.

    <pre>$pgcgap --pCOG --threads 4 --strain_num 6 --AAsPath Results/Annotations/AAs</pre>

  - __Example 11:__ Variants calling and phylogenetic tree construction based on the reference genome.

    <pre># Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap)<br/>$pgcgap --VAR --threads 4 --refgbk /mnt/h/PGCGAP_Examples/Reads/MG1655.gbff --ReadsPath Reads/Illumina --reads1 _1.fastq.gz --reads2 _2.fastq.gz --suffix_len 11 --strain_num 6 --qualtype sanger --bsnum 500<br/># Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)<br/>$pgcgap --VAR --threads 4 --refgbk /mnt/h/PGCGAP_Examples/Reads/MG1655.gbff --ReadsPath Reads/Illumina --reads1 _1.fastq.gz --reads2 _2.fastq.gz --suffix_len 11 --strain_num 6 --qualtype sanger --fastboot 1000</pre>

  - __Example 12:__ Screening of contigs for antimicrobial and virulence genes

    <pre>$pgcgap --AntiRes --scafPath Results/Assembles/Scaf/Illumina --Scaf_suffix -8.fa --threads 6 --db all --identity 75 --coverage 50</pre>

  - __Example 13:__ Filter short sequences in the genome and assess the status of the genome

    <pre>$pgcgap --ACC --Assess --scafPath Results/Assembles/Scaf/Illumina --Scaf_suffix -8.fa --filter_length 200</pre>

  - __Example 14:__ Construct a phylogenetic tree based on multiple sequences in one file

    <pre># Construct phylogenetic tree with IQ-TREE (Very slow with best fit model testing, traditional bootstrap)<br/>$pgcgap --STREE --seqfile proteins.fas --seqtype p --bsnum 500 --threads 4<br/># Construct phylogenetic tree with IQ-TREE (Slow with best fit model testing, ultrafast bootstrap)<br/>$pgcgap --STREE --seqfile proteins.fas --seqtype p --fastboot 1000 --threads 4</pre>


## Generating Input files

### Working directory
The directory where the PGCGAP software runs.

### Assemble
Pair-end reads of all strains in a directory or PacBio reads or Oxford nanopore reads (Default: ./Reads/Illumina/ under the working directory).

### Annotate
Genomes files (complete or draft) in a directory (Default: Results/Assembles/Scaf/Illumina under the working directory).

### ANI

QUERY\_LIST and REFERENCE\_LIST files containing full paths to genomes, one per line (default: scaf.list under the working directory). If the “\-\-Assemble” function was run first, the list file will be generated automatically.

### MASH

Genomes files (complete or draft) in a directory (Default: Results/Assembles/Scaf/Illumina under the working directory).

### CoreTree
Amino acids file (With “.faa” as the suffix) and nucleotide (With “.ffn” as the suffix) file of each strain placed into two directories (default: “./Results/Annotations/AAs/” and “./Results/Annotations/CDs/”). The “.faa” and “.ffn” files of the same strain should have the same prefix name. The name of protein IDs and gene IDs should be started with the strain name. The “Prokka” software was suggested to generate the input files. If the “\-\-Annotate” function was run first, the files will be generated automatically. If the “\-\-CDsPath” was set to “NO”, the nucleotide files will not be needed.

### OrthoF
A set of protein sequence files (one per species) in FASTA format under a directory (default: “./Results/Annotations/AAs/”). If the “\-\-Annotate” function was run first, the files will be generated automatically.

### Pan
GFF3 files (With “.gff” as the suffix) of each strain placed into a directory. They must contain the nucleotide sequence at the end of the file. All GFF3 files created by Prokka are valid (default: ./Results/Annotations/GFF/). If the “\-\-Annotate” function was run first, the files will be generated automatically.

### pCOG
Amino acids file (With “.faa” as the suffix) of each strain placed into a directory (default: ./Results/Annotations/AAs/). If the “\-\-Annotate” function was run first, the files will be generated automatically.

### VAR
- Pair-end reads of all strains in a directory (default: ./Reads/Over/ under the working directory).
<br/>

- The full path of reference genome in fasta format or GenBank format (__must be provided__).

### AntiRes
Genomes files (complete or draft) in a directory (Default: Results/Assembles/Scaf/Illumina under the working directory).

### STREE
Multiple-FASTA sequences in a file, can be Protein, DNA and Codons.

## Output Files

### Assemble

- **Results/Assembles/Illumina/**<br/>
Directories contain Illumina assembly files and information of each strain.
<br/>

- **Results/Assembles/PacBio/**<br/>
Directories contain PacBio assembly files and information of each strain.
<br/>

- **Results/Assembles/Oxford/**<br/>
Directories contain ONT assembly files and information of each strain.
<br/>

- **Results/Assembles/Hybrid/**<br/>
Directory contains hybrid assembly files of the short reads and long reads of the same strain.
<br/>

- __Results/Assembles/Scaf/Illumina__<br/>
Directory contains Illumina contigs/scaffolds of all strains. "\*.filtered.fas" is the genome after excluding short sequences. "\*.prefilter.stats" describes the stats of the genome before filtering, and "\*.filtered.stats" describes the stats of the genome after filtering.
<br/>

- __Results/Assembles/Scaf/Oxford__<br/>
Directory contains Oxford nanopore contigs/scaffolds of all strains.
<br/>

- __Results/Assembles/Scaf/PacBio__<br/>
Directory contains PacBio contigs/scaffolds of all strains.
<br/>

### Annotate

- **Results/Annotations/\*\_annotation**<br/>
directories contain [annotation files](https://github.com/tseemann/prokka?_blank) of each strain.
<br/>

- __Results/Annotations/AAs__<br/>
Directory contain amino acids sequences of all strains.
<br/>

- __Results/Annotations/CDs__<br/>
Directory contain nucleotide sequences of all strains.
<br/>

- __Results/Annotations/GFF__<br/>
Directory contain the master annotation of all strains in GFF3 format.


### ANI

- __Results/ANI/ANIs__<br/>
The file contains comparation information of genome pairs. The document is composed of five columns, each of which represents query genome, reference genome, ANI value, count of bidirectional fragment mappings, total query fragments.
<br/>

- __Results/ANI/ANIs.matrix__<br/>
file with identity values arranged in a [phylip-formatted lower triangular matrix](https://www.mothur.org/wiki/Phylip-formatted_distance_matrix?_blank).
<br/>

- __Results/ANI/ANIs.heatmap__<br/>
An ANI matrix of all strains.
<br/>

- __Results/ANI/ANI\_matrix.pdf__<br/>
The heatmap plot of "ANIs.heatmap".

### MASH

- __Results/MASH/MASH__<br/>
The pairwise distance between pair genomes, each column represents Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes, respectively.

- __Results/MASH/MASH2__<br/>
The pairwise similarity between pair genomes, each column represents Reference-ID, Query-ID, similarity, P-value, and Matching-hashes, respectively.

- __Results/MASH/MASH.heatmap__<br/>
A similarity matrix of all genomes.

- __Results/MASH/MASH\_matrix.pdf__<br/>
A heat map plot of "MASH.heatmap".

### CoreTree

- __Results/CoreTrees/ALL.core.protein.fasta__<br/>
Concatenated and aligned sequences file of single-copy core proteins.
<br/>

- __Results/CoreTrees/ALL.core.protein.nwk__<br/>
The phylogenetic tree file of single-copy core proteins for all strains constructed by FastTree.

- __Results/CoreTrees/ALL.core.protein.fasta.gb.treefile__<br/>
The phylogenetic tree file of single-copy core proteins for all strains constructed by IQ-TREE.

- __Results/CoreTrees/faa2ffn/ALL.core.nucl.fasta__<br/>
Concatenated and aligned sequences file of single-copy core genes.
<br/>

- __Results/CoreTrees/ALL.core.snp.fasta__<br/>
Core SNPs of single-copy core genes in fasta format.
<br/>

- __Results/CoreTrees/ALL.core.snp.nwk__<br/>
The phylogenetic tree file of SNPs of single-copy core genes for all strains constructed by FastTree.

- __Results/CoreTrees/ALL.core.snp.fasta.gb.treefile__<br/>
The phylogenetic tree file of SNPs of single-copy core genes for all strains constructed by IQ-TREE
<br/>

- __Results/CoreTrees/"Other\_files"__<br/>
Intermediate directories and files.
<br/>

### OrthoF

- __Results/OrthoFinder/Results\_orthoF__<br/>
Same as [OrthoFinder](https://github.com/davidemms/OrthoFinder?_blank) outputs.
<br/>
- __Results/OrthoFinder/Results_orthoF/Single_Copy_Orthologue_Tree/__<br/>
Directory contains Phylogenetic tree files based on Single Copy Orthologue sequences.<br/>

- __Results/OrthoFinder/Results_orthoF/Single_Copy_Orthologue_Tree/Single.Copy.Orthologue.nwk__<br/>
Phylogenetic tree constructed by FastTree.<br/>

- __Results/OrthoFinder/Results_orthoF/Single_Copy_Orthologue_Tree/Single.Copy.Orthologue.fasta.gb.treefile__<br/>
Phylogenetic tree constructed by IQ-TREE.<br/>

### Pan

- __Results/PanGenome/Pangenome\_Pie.pdf__<br/>
A 3D pie chart and a fan chart of the breakdown of genes and the number of isolates they are present in.
<br/>

- __Results/PanGenome/pangenome\_frequency.pdf__<br/>
A graph with the frequency of genes versus the number of genomes.
<br/>

- __Results/PanGenome/Pangenome\_matrix.pdf__<br/>
A figure showing the tree compared to a matrix with the presence and absence of core and accessory genes.
<br/>

- __Results/PanGenome/Core/Roary.core.protein.fasta__<br/>
Alignments of single-copy core proteins called by roary software.
<br/>

- __Results/PanGenome/Core/Roary.core.protein.nwk__<br/>
A phylogenetic tree of Roary.core.protein.fasta constructed by FastTree.
<br/>

- __Results/PanGenome/Core/Roary.core.protein.fasta.gb.treefile__<br/>
A phylogenetic tree of Roary.core.protein.fasta constructed by IQ-TREE.
<br/>

- __Results/PanGenome/Other\_files__<br/>
see [roary](https://sanger-pathogens.github.io/Roary/?_blank) outputs.
<br/>

### pCOG

- __\*.COG.xml, \*.2gi.table, \*.2id.table, \*.2Sid.table__<br/>
Intermediate files.
<br/>

- __\*.2Scog.table__<br/>
The super COG table of each strain.
<br/>

- __\*.2Scog.table.pdf__<br/>
A plot of super COG table in pdf format.
<br/>

- __All\_flags\_relative\_abundances.table__
A table containing the relative abundance of each flag for all strains.

### VAR
- __Results/Variants/directory-named-in-strains__<br/>
directories containing substitutions (snps) and insertions/deletions (indels) of each strain. See [Snippy](https://github.com/tseemann/snippy?_blank) outputs for detail.

- __Results/Variants/Core__<br/>
  The directory containing SNP phylogeny files.

  - __core.aln__ : A core SNP alignment includes only SNP sites.
  - __core.full.aln__ : A whole genome SNP alignment (includes invariant sites).
  - __core.*.treefile__ : Phylogenetic tree of the core SNP alignment based on the best-fit model of evolution selected using IQ-TREE (ignoring possible recombination).
  - **gubbins.core.full.node\_labelled.final\_tree.tre** : Phylogenetic tree of the whole genome SNP alignment constructed with __gubbins__ (get rid of recombination).

### AntiRes
- __Results/AntiRes/*.tab__ : Screening results of each strain.
- __Results/AntiRes/summary.txt__ : A matrix of gene presence/absence for all strains.

### STREE
- __Results/STREE/*.aln__ : Aligned sequences.
- __Results/STREE/*.aln.gb__ : Conserved blocks of the sequences.
- __Results/STREE/*.aln.gb.treefile__ : The final phylogenetic tree.

## License

PGCGAP is free software, licensed under GPLv3.

## Feedback and Issues

Please report any issues to the [issues page](https://github.com/liaochenlanruo/pgcgap/issues?_blank) or email us at [liaochenlanruo@webmail.hzau.edu.cn](mailto:liaochenlanruo@webmail.hzau.edu.cn).

## Citation

- If you use this software please cite: Liu H, Xin B, Zheng J, Zhong H, Yu Y, Peng D, Sun M. Build a bioinformatics analysis platform and apply it to routine analysis of microbial genomics and comparative genomics. *Protocol exchange*, 2020. DOI: [10.21203/rs.2.21224/v5](https://doi.org/10.21203/rs.2.21224/v5)

- If you use "\-\-Assemble", please also cite one or two of [Fastp](https://github.com/OpenGene/fastp), [ABySS](https://doi.org/10.1101/gr.214346.116), [SPAdes](http://link.springer.com/chapter/10.1007%2F978-3-642-37195-0_13), [Canu](https://doi.org/10.1101/gr.215087.116), or [Unicycler](https://doi.org/10.1371/journal.pcbi.1005595).

- If you use "\-\-Annotate", please also cite [Prokka](https://www.pixiv.net/member_illust.php?mode=medium&illust_id=24642063).

- If you use "\-\-CoreTree", please also cite [CD-HIT](https://doi.org/10.1093/bioinformatics/btl158), [MAFFT](https://doi.org/10.1093/nar/gkf436), [PAL2NAL](https://doi.org/10.1093/nar/gkl315), [trimAL](https://doi.org/10.1093/bioinformatics/btp348), [FastTree](https://doi.org/10.1371/journal.pone.0009490) or [IQ-TREE](https://doi.org/10.1093/molbev/msaa015), and [SNP-sites](https://dx.doi.org/10.1099%2Fmgen.0.000056).

- If you use "\-\-Pan", please also cite [Roary](https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtv421), [MAFFT](https://doi.org/10.1093/nar/gkf436), [trimAL](https://doi.org/10.1093/bioinformatics/btp348), [FastTree](https://doi.org/10.1371/journal.pone.0009490) or [IQ-TREE](https://doi.org/10.1093/molbev/msaa015).

- If you use "\-\-OrthoF", please also cite [OrthoFinder](https://dx.doi.org/10.1186%2Fs13059-019-1832-y), [trimAL](https://doi.org/10.1093/bioinformatics/btp348), [FastTree](https://doi.org/10.1371/journal.pone.0009490) or [IQ-TREE](https://doi.org/10.1093/molbev/msaa015).

- If you use "\-\-ANI", please also cite [fastANI](https://dx.doi.org/10.1038%2Fs41467-018-07641-9).

- If you use "\-\-MASH", please also cite [Mash](https://dx.doi.org/10.1186%2Fs13059-016-0997-x).

- If you use "\-\-VAR", please also cite [Sickle](https://github.com/najoshi/sickle), [Snippy](https://github.com/tseemann/snippy), [Gubbins](https://dx.doi.org/10.1093%2Fnar%2Fgku1196),  [IQ-TREE](https://doi.org/10.1093/molbev/msaa015), and [SnpEff](https://dx.doi.org/10.4161%2Ffly.19695).

- If you use "\-\-AntiRes", please also cite [Abricate](https://github.com/tseemann/abricate) and the corresponding database you used: [NCBI AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6811410), [CARD](https://www.ncbi.nlm.nih.gov/pubmed/27789705), [Resfinder](https://www.ncbi.nlm.nih.gov/pubmed/22782487), [ARG-ANNOT](https://www.ncbi.nlm.nih.gov/pubmed/24145532), [VFDB](https://www.ncbi.nlm.nih.gov/pubmed/26578559), [PlasmidFinder](https://www.ncbi.nlm.nih.gov/pubmed/24777092), [EcOH](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5343136/), or [MEGARES 2.00](https://academic.oup.com/nar/article/48/D1/D561/5624973).

- If you use "\-\-STREE", please also cite [Muscle](http://europepmc.org/abstract/MED/30976793) , [trimAL](https://doi.org/10.1093/bioinformatics/btp348), and [IQ-TREE](https://doi.org/10.1093/molbev/msaa015).


## FAQ


### Q1 VAR function ran failed to get annotated VCFs and Core results


Check the log file named in "strain_name.log" under Results/Variants/<strain\_name\>/ directory. If you find a sentence like "WARNING: All frames are zero! This seems rather odd, please check that 'frame' information in your 'genes' file is accurate." This is a snpEff error. Users can install JDK8 to solve this problem.



<pre>$conda install java-jdk=8.0.112</pre>



Click [here](https://github.com/tseemann/snippy/issues/259?_blank) for more solutions.

### Q2 Could not determine version of minced please install version 2 or higher
When running the Annotate function, this error could happen, the error message shows as following:

<pre>
Error: A JNI error has occurred, please check your installation and try again
Exception in thread "main" java.lang.UnsupportedClassVersionError: minced has been compiled by a more recent version of the Java Runtime (class file version 55.0), this version of the Java Runtime only recognizes class file versions up to 52.0
	at java.lang.ClassLoader.defineClass1(Native Method)
	at java.lang.ClassLoader.defineClass(ClassLoader.java:763)
	at java.security.SecureClassLoader.defineClass(SecureClassLoader.java:142)
	at java.net.URLClassLoader.defineClass(URLClassLoader.java:468)
	at java.net.URLClassLoader.access$100(URLClassLoader.java:74)
	at java.net.URLClassLoader$1.run(URLClassLoader.java:369)
	at java.net.URLClassLoader$1.run(URLClassLoader.java:363)
	at java.security.AccessController.doPrivileged(Native Method)
	at java.net.URLClassLoader.findClass(URLClassLoader.java:362)
	at java.lang.ClassLoader.loadClass(ClassLoader.java:424)
	at sun.misc.Launcher$AppClassLoader.loadClass(Launcher.java:349)
	at java.lang.ClassLoader.loadClass(ClassLoader.java:357)
	at sun.launcher.LauncherHelper.checkAndLoadMain(LauncherHelper.java:495)
[01:09:40] Could not determine version of minced - please install version 2.0 or higher
</pre>
Users can downgrade the minced to version 0.3 to solve this problem.

<pre>$conda install minced=0.3</pre>

Click [here](https://github.com/bioconda/bioconda-recipes/pull/15407?_blank) for detail informations.

### Q3 dyld: Library not loaded: @rpath/libcrypto.1.0.0.dylib

This error may happen when running function "VAR" on macOS. It is an error of openssl. Users can solve this problem as the following:

<pre>
#Firstly, install brew if have not installed before
$ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

#Install openssl with brew
$brew install openssl

#Create the soft link for libraries
$ln -s /usr/local/opt/openssl/lib/libcrypto.1.0.0.dylib /usr/local/lib/

$ln -s /usr/local/opt/openssl/lib/libssl.1.0.0.dylib /usr/local/lib/
</pre>

Click [here](https://gist.github.com/aklap/e885721ef15c8668ed0a1dd64d2ea1a7) for more informations

### Q4 Use of uninitialized value in require at Encode.pm line 61

This warning may happen when running function "Pan". It is a warning of Roary software.
The content of line 61 is "require Encode::ConfigLocal;". Users can ignore the warning.
Click [here](https://github.com/sanger-pathogens/Roary/issues/323) for details.

## Updates

- V1.0.3

  - Updated ANI function.

- V1.0.4

  - Add parallel for function "pCOG".
  - Optimized drawing of ANI heat map.

- V1.0.5

  - Bug repair for the input of gubbins.

- V1.0.6

  - Modified CoreTree to split protein and SNPs tree constructing.

- V1.0.7

  - Split Assemble and Annotate into two functions.
  - Added third-generation genome assembly function.
  - Changed the default parameters of the CoreTree function (aS 0.8 to 0.7 and aL 0.8 to 0.5).
  - Changed the name of function "COG" to "pCOG".
  - Fixed the sorting bug for ANI heat map.

- V1.0.8

  - Add the "MASH" function to compute genome distance and similarity using MinHash.

- V1.0.9

  - The function of constructing a single-copy core protein phylogenetic tree was added to "Pan".
  - Fixed a bug of plot_3Dpie.R, Optimized image display, and a fan chart has been added.
  - Fixed a bug for plotting the ANI matrix.

- V1.0.10

  - Add the "AntiRes" function to screening of contigs for antimicrobial and virulence genes.

- V1.0.11

  - Users now can choose "abyss" or "spades" for illumina reads aseembly.
  - New support for hybrid assembly of paired-end short reads and long reads.
  - Add the selecting of best-fit model of evolution for DNA and protein alignments before constructing a phylogenetic tree.
  - Optimized display of help information. Users can check parameters for each modulewith command "pgcgap \[Assemble|Annotate|ANI|AntiRes|CoreTree|MASH|OrthoF|Pan|pCOG|VAR\]", and can look up the examples of each module with command "pgcgap Examples".

- V1.0.12

  - Added automatic mode for illumina genome assembly. First, PGCGAP calls "ABySS" for genome assembly. When the assembled N50 is less than 50,000, it automatically calls "SPAdes" to try multiple parameters for assembly.
  - Added ability to filter short sequences of assembled genomes.
  - Added function of genome assembly status assessment.
  - Modified the drawing script of ANI and MASH modules so that it can automatically adjust the font size according to the number of samples.

- V1.0.13

  - Fixed the "running error" bug of function "Assess" in module "ACC".
  - Added module "STREE" for constructing a phylogenetic tree based on multiple sequences in one file.

- V1.0.14
  - The relative_abundances of flags among strains will not be called while the strain number is less than two.
  - Fixed the error of function "Assess" in module "ACC".

- V1.0.15
  - When the number of threads set by the user exceeds the number of threads owned by the system, PGCGAP will automatically adjust the number of threads to avoid program crash.
  - Add FASTQ preprocessor before Illunima genome assembly: adapter trimming, polyG tail trimming of Illumina NextSeq/NovaSeq reads, quality filtering (Q value filtering, N base filtering, sliding window filtering), length filtering.

- V1.0.16
  - Reduced the number of Racon polishing rounds for better speed performance when peforming genome assembly.
  - Force overwriting existing output folder when running "Annotate" analysis to avoid program crash.

- V1.0.17
  - Fixed a bug that the program can not go back to the working directory after genome annotation.
  - Added scripts to check if there were single-copy core proteins found while running module "CoreTree".
  - Modified the help message.

- V1.0.18
  - Updated the downloading link of COG database.
  - Users can choose the number of threads used for running module "STREE".

- V1.0.19
  - Can resume from break-point when downloading the COG database.
  - Fixed a bug that failed to create multi-level directories.

- V1.0.20
  - Fixed a little bug (path error) of module "VAR".
  - Fixed a little bug of module "CoreTree" to avoid the interference of special characters in sequence ID to the program.

- V1.0.21
  - Change the default search program "blast" to "diamond" of module "OrthoF".
  - Fixed a bug of module "pCOG" to output the right figure.

- V1.0.22
  - The drawing function of module "ANI" and module "MASH" has been enhanced, including automatic adjustment of font size and legend according to the size of the picture.
  - Fixed a bug of module "ANI", that is no heatmap will be drawn when there is "NA" in the ANI matrix in the previous versions.
  - When the ANI value or genome similarity is greater than 95%, an asterisk (*) will be drawn in the corresponding cell of the heatmap.

- V1.0.23
  - The "--Assess" function of module "ACC" was enhanced to (1) generate a summary file containing the status of all genomes (before and after the short sequence filtering), (2) auto move the low-quality genomes (that is genomes with N50 length less than 50 k) to a directory, and others to another directory.

- V1.0.24
  - Fixed a little bug of module "Pan" to avoid the interference of special characters (>) in sequence ID to the program.

- V1.0.25
  - Gblocks was used to eliminate poorly aligned positions and divergent regions of an alignment of DNA or protein sequences in module "CoreTree" and "Pan".
  - The parameter "--identi" was added into module "Pan" to allow users to set the minimum percentage identity for blastp.
 
 - V1.0.26
   - Adjusted the font size with the variation of genome number and the string length of the genome name when plotting the heat map of module "ANI" and "MASH".
   - Two heat map are provided, one of which with a star (means the similarity of the two genomes is larger than 95%) and another without a star, when performing the "ANI" and "MASH" analysis.

 - V1.0.27
   - The Amino Acid files are no longer needed when performing the Pan-genome analysis with module Pan.

 - V1.0.28
   - Users can check and install the latest version of PGCGAP by the command "pgcgap --check-update".
   - Update module Assemble to allow polish after the assembly of PacBio and ONT data.
   - Update module pCOG to adjust the latest database of [COG 2020](https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/COG).
   - Optimized the drawing and color scheme of the module pCOG.
   - Fixed the parameter "CoreTree" in the module Pan to avoid program termination caused by the ">" in non-sequence lines.

- V1.0.29
  - Function added to module OrthoF: Phylogenetic tree can be constructed automatically with the Single Copy Orthologue Sequences called by module OrthoF.
  - Fixed the "permission denied" error when moving directories on the WSL platform.


- V1.0.30
 - Replace Gblocks with trimAL to trim MSA (module CoreTree, Pan, STREE, and OrthoF).
 - Replaced Modeltest-ng and Raxml-ng with IQ-TREE (module CoreTree, Pan OrthoF, and VAR).
 - Added the option of using FastTree to build phylogenetic tree (module CoreTree, Pan, and OrthoF).

- V1.0.31
  - The default replicates for bootstrap testing was set to 500.
  - Add the method for phylogenetic tree constructing with ultrafast bootstrap of IQ-TREE.
  - Prevent the log from being written to the tree file generated by FastTree.

- V1.0.32
  - A more colorful version, try "pgcgap Examples" to have a look.
  - Updated module AntiRes: the parameter --db had been modified to add choices of "all" and "megares".
  - A little optimization of module VAR.
  - Replaced conda with mamba to update PGCGAP more quickly.
---

<center><strong>
<script async src="//busuanzi.ibruce.info/busuanzi/2.3/busuanzi.pure.mini.js"></script>
<span id="busuanzi_container_site_pv">Total visits: <span id="busuanzi_value_site_pv"></span> times</span>
<span class="post-meta-divider">|</span>
<span id="busuanzi_container_site_uv">Visitors: <span id="busuanzi_value_site_uv"></span> people</span>
</strong></center>