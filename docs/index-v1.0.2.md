---
sort: 998
---

# DOC for V1.0.2
---

![Platform](https://img.shields.io/badge/Platform-WSL%2FLinux%2FmacOS-green) [![License](https://img.shields.io/github/license/liaochenlanruo/pgcgap)](https://github.com/liaochenlanruo/pgcgap/blob/master/LICENSE) [![GitHubversion](https://anaconda.org/bioconda/pgcgap/badges/version.svg)](https://anaconda.org/bioconda/pgcgap) ![Downloads conda](https://img.shields.io/conda/dn/bioconda/pgcgap.svg?style=flat) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pgcgap/README.html) [![Bioinformatics](https://pub.idqqimg.com/wpa/images/group.png "945751012")](//shang.qq.com/wpa/qunwpa?idkey=fd4637eecd73bf0a5a8caa274843a07afdf1fbbc40a86630df5d4b029749cc7b)


<p><center>
<script type="text/javascript" src="//rf.revolvermaps.com/0/0/4.js?i=0ett3z77b0x&amp;m=0&amp;h=128&amp;c=ff0000&amp;r=0" async="async"></script>
&nbsp;&nbsp;&nbsp;&nbsp;
<script type="text/javascript" src="//rf.revolvermaps.com/0/0/0.js?i=0ett3z77b0x&amp;d=3&amp;p=1&amp;b=1&amp;w=293&amp;g=2&amp;f=arial&amp;fs=13&amp;r=0&amp;c0=ff8a00&amp;c1=0006ff&amp;c2=000000&amp;ic0=0&amp;ic1=0" async="async"></script>
</center></p>
-------------

[English Readme](https://liaochenlanruo.github.io/pgcgap) \| [Chinese](https://liaochenlanruo.github.io/post/848f.html)


          ____       ____      ____     ____       _        ____    
        U|  _"\ u U /"___|u U /"___| U /"___|u U  /"\  u  U|  _"\ u 
        \| |_) |/ \| |  _ / \| | u   \| |  _ /  \/ _ \/   \| |_) |/ 
         |  __/    | |_| |   | |/__   | |_| |   / ___ \    |  __/   
         |_|        \____|    \____|   \____|  /_/   \_\   |_|      
         ||>>_      _)(|_    _// \\    _)(|_    \\    >>   ||>>_    
        (__)__)    (__)__)  (__)(__)  (__)__)  (__)  (__) (__)__)   

-------------

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Required dependencies](#required-dependencies)
- [Usage](#usage)
- [Generating Input files](#generating-input-files)
  * [Working directory](#working-directory)
  * [Assemble](#assemble)
  * [ANI](#ani)
  * [CoreTree](#coretree)
  * [OrthoF](#orthof)
  * [Pan](#pan)
  * [COG](#cog)
  * [VAR](#var)
- [Output Files](#output-files)
  * [Assemble](#assemble-1)
  * [ANI](#ani-1)
  * [CoreTree](#coretree-1)
  * [OrthoF](#orthof-1)
  * [Pan](#pan-1)
  * [COG](#cog)
  * [VAR](#var-1)
- [License](#license)
- [Feedback and Issues](#feedback-and-issues)
- [Citation](#citation)
- [FAQ](#faq)
  * [Q1 VAR FOUNCTION ran failed to get annotated VCFs and Core results](#q1-var-founction-ran-failed-to-get-annotated-vcfs-and-core-results)




## Introduction

PGCGAP is a pipeline for prokaryotic comparative genomics analysis. It can take the pair-end reads as input. In addition to genome assembly, gene prediction and annotation, it can also get common comparative genomics analysis results such as phylogenetic trees of single-core proteins and core SNPs, pan-genome, whole-genome Average Nucleotide Identity (ANI), orthogroups and orthologs, COG annotations, substitutions (snps) and insertions/deletions (indels) with only one line of commands.

## Installation

The software was only tested on Linux64 platform, The MacOS may also be good in theory. However, Windows could not be supported. Because this software relies on a large number of other softwares, so it is recommended to install with __[Bioconda](https://bioconda.github.io/index.html)__. The main program and most of other dependencies can be installed with one command as shown below, __but the "Gubbins" should be installed separately__ because of it relys on python3, while PGCGAP relys on python2.


__Step1: Install Gubbins__

```
$conda install gubbins
```

__Step2: Install PGCGAP__

```
$conda create -n pgcgap python=2.7 anaconda
$conda activate pgcgap
$conda install pgcgap # By command "whereis pgcgap", users can find the installation path of pgcgap, and then you should add it into your environment variable.
$conda deactivate
```

__Step3: Setup COG database__ (Users should execute this after first installation of pgcgap)

```
$pgcgap --setup-COGdb
```

## Required dependencies


- [ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss/)
- [CD-HIT](http://weizhongli-lab.org/cd-hit/)
- [Diamond](https://github.com/bbuchfink/diamond)
- [FastANI](https://github.com/ParBLiSS/FastANI)
- [Fastme](http://www.atgc-montpellier.fr/fastme/binaries.php)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [Gubbins](https://github.com/sanger-pathogens/gubbins) >=2.3.4
- [Htslib](https://github.com/samtools/htslib)
- [Mafft](https://mafft.cbrc.jp/alignment/software/)
- [Mmseqs2](https://github.com/soedinglab/mmseqs2)
- [NCBI-blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [OrthoFinder](https://github.com/davidemms/OrthoFinder)
- [OpenJDK8](https://openjdk.java.net/)
- [PAL2NAL v14](http://www.bork.embl.de/pal2nal/)
- [Perl](http://www.perl.org/get.html) & the modules
  - [perl-bioperl](http://metacpan.org/pod/BioPerl)
  - [perl-data-dumper](http://metacpan.org/pod/Data::Dumper)
  - [perl-file-tee](http://metacpan.org/pod/File::Tee)
  - [perl-getopt-long](http://metacpan.org/pod/Getopt::Long)
  - [perl-pod-usage](http://search.cpan.org/~marekr/Pod-Usage-1.69/)
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
  - [plotrix](https://cran.r-project.org/web/packages/plotrix/)
- [Roary](https://sanger-pathogens.github.io/Roary/)
- [Sickle-trim](https://github.com/najoshi/sickle)
- [Snippy](https://github.com/tseemann/snippy)
- [Snp-sites](https://github.com/sanger-pathogens/snp-sites)

## Usage


- __Print the help messages:__
    ```
    $pgcgap --help
    ```
<br/>

- __General usage:__
    ```
    $pgcgap [Fuctions] [options]
    ```
<br/>

- __Setup COG database:__ (Users should execute this after first installation of pgcgap)
    ```
    $pgcgap --setup-COGdb
    ```
<br/>

- __FUCTIONS:__

  - __[--All]__                          Perform all functions with one command

  - __[--Assemble]__                     Assemble reads ( paired-end only ) into contigs,
                                         predict genes and annotate them

  - __[--CoreTree]__                     Construct single-core proteins tree and SNPs tree of single core genes

  - __[--Pan]__                          Run "roary" pan genome pipeline with gff3 files

  - __[--OrthoF]__                       Identify orthologous protein sequence
                                         families with "OrthoFinder"

  - __[--ANI]__                          Compute whole-genome Average Nucleotide Identity ( ANI )

  - __[--COG]__                          Run COG annotation for each strain (*.faa),
				                                 and generate a table containing the relative
				                                 abundance of each flag for all strains

  - **[--VAR]**                          Rapid haploid variant calling and core
                                         genome alignment with "Snippy"
<br/>

- __Global Options:__

  - __[--strain_num (INT)]__             [Required by "--All", "--CoreTree" and "--VAR"]
                                         The total number of strains used for analysis,
                                         including reference genomes

  - __[--ReadsPath (PATH)]__             [Required by "--All", "--Assemble" and "--VAR"]
                                         Reads of all strains as file paths ( Default ./Reads )

  - __[--AAsPath (PATH)]__               [Required by "--CoreTree" and "--COG"]
                                         Amino acids of all strains as fasta file paths,
                                         ( Default "./Results/Annotations/AAs" )

  - __[--reads1 (STRING)]__              [Required by "--All", "--Assemble" and
                                         "--VAR"] The suffix name of reads 1 ( for
                                         example: if the name of reads 1 is
                                         "YBT-1520_L1_I050.R1.clean.fastq.gz",
                                         "YBT-1520" is the strain same, so the
                                         suffix name should be ".R1.clean.fastq.gz")

  - __[--reads2 (STRING)]__              [Required by "--All", "--Assemble" and
                                         "--VAR"] The suffix name of reads 2( for
                                         example: if the name of reads 2 is
                                         "YBT-1520_2.fq", the suffix name should be "_2.fq" )

  - __[--codon (INT)]__                  [Required by "--All", "--Assemble", "--CoreTree" and "--Pan"] Translation table ( Default 11 )

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

  - __[--suffix_len (INT)]__             [Required by "--All", "--Assemble" and
                                         "--VAR"] __(Strongly recommended)__ The suffix length of the reads,
                                         that is the length of your reads name
                                         minus the length of your strain name. For
                                         example the --suffix_len of
                                         "YBT-1520_L1_I050.R1.clean.fastq.gz" is 26
                                         ( "YBT-1520" is the strain name ) ( Default 0 )

  - __[--logs (STRING)]__                Name of the log file ( Default Logs.txt )

  - __[--threads (INT)]__                Number of threads to be used ( Default 4 )
<br/>

- __Local Options:__

  - __--Assemble__
      - __[--kmmer (INT)]__              [Required] k-mer size for genome assembly ( Default 81 )

      - __[--genus (STRING)]__           Genus name of your strain ( Default "NA" )

      - __[--species (STRING)]__         Species name of your strain ( Default "NA")
  <br/>

  - **--CoreTree**

      - __[--CDsPath (PATH)]__           [Required] CDs of all strains as fasta file
                                         paths, ( Default "./Results/Annotations/CDs" )

      - __[-c (FLOAT)]__                 Sequence identity threshold, ( Default 0.5)

      - __[-n (INT)]__                   Word_length, see user's guide of CD-HIT
                                         for choosing it ( Default 2 )

      - __[-G (INT)]__                   Use global (set to 1) or local (set to 0)
                                         sequence identity, ( Default 0 )

      - __[-t (INT)]__                   Tolerance for redundance ( Default 0 )

      - __[-aL (FLOAT)]__                Alignment coverage for the longer
                                         sequence. If set to 0.9, the alignment
                                         must covers 90% of the sequence ( Default 0.9 )

      - __[-aS (FLOAT)]__                Alignment coverage for the shorter sequence.
                                         If set to 0.9, the alignment must covers
                                         90% of the sequence ( Default 0.9 )

      - __[-g (INT)]__                   If set to 0, a sequence is clustered to
                                         the first cluster that meet the threshold
                                         (fast cluster). If set to 1, the program
                                         will cluster it into the most similar
                                         cluster that meet the threshold (accurate
                                         but slow mode, Default 1)

      - __[-d (INT)]__                   length of description in .clstr file. if
                                         set to 0, it takes the fasta defline and
                                         stops at first space ( Default 0 )
  <br/>

  - __--Pan__

      - __[--GffPath (PATH)]__           [Required] Gff files of all strains as paths
                                         ( Default "./Results/Annotations/GFF" )
  <br/>

  - __--orthoF__

      - __[--Sprogram (STRING)]__        Sequence search program, Options: blast,
                                         mmseqs, blast_gz, diamond ( Default blast)
  <br/>

  - __--ANI__

      - __[--queryL (FILE)]__            [Required] The file containing paths to query
                                         genomes, one per line ( Default scaf.list )

      - __[--refL (FILE)]__              [Required] The file containing paths to reference genomes,
                                         one per line. ( Default scaf.list )

      - __[--ANIO (FILE)]__              The name of output file ( Default "Results/ANI/ANIs" )
  <br/>

  - **--VAR**

      - __[--refgbk (FILE)]__            [Required] The full path and name of
                                         reference genome in GENBANK format (
                                         recommended ), fasta format is also OK.
                                         For example: "/mnt/g/test/ref.gbk"

      - __[--qualtype (STRING)]__        [Required] Type of quality values (solexa
                                         (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7),
				         sanger (which is CASAVA >= 1.8)). ( Default sanger )

      - __[--qual (INT)]__               Threshold for trimming based on average
                                         quality in a window. ( Default 20 )

      - __[--length (INT)]__             Threshold to keep a read based on length
                                         after trimming. ( Default 20 )

      - __[--mincov (INT)]__             The minimum number of reads covering a
                                         site to be considered ( Default 10 )

      - __[--minfrac (FLOAT)]__          The minimum proportion of those reads
                                         which must differ from the reference ( Default 0.9 )

      - __[--minqual (INT)]__            The minimum VCF variant call "quality" ( Default 100 )

      - __[--ram (INT)]__                Try and keep RAM under this many GB ( Default 8 )

      - __[--tree_builder (STRING)]__    Application to use for tree building
                                         [raxml|fasttree|hybrid] ( Default fasttree)

      - __[--iterations (INT)]__         Maximum No. of iterations for gubbins ( Default 5 )
<br/>

- __Paths of external programs__

  Not needed if they were in the environment variables path. Users can check with the "--check-external-programs" option for the essential programs.
  <br/>

  - __[--abyss-bin (PATH)]__             Path to abyss binary file. Default tries if abyss is in PATH;

  - __[--prodigal-bin (PATH)]__          Path to prodigal binary file. Default tries if prodigal is in PATH;

  - __[--prokka-bin (PATH)]__            Path to prokka binary file. Default tries if prokka is in PATH;

  - __[--cd-hit-bin (PATH)]__            Path to cd-hit binary file. Default tries if cd-hit is in PATH;

  - __[--mafft-bin (PATH)]__             Path to mafft binary file. Default tries if mafft is in PATH;

  - __[--fasttree-bin (PATH)]__        Path to the fasttree binary file. Default tries if fasttree is in PATH;

  - __[--pal2nal-bin (PATH)]__           Path to the pal2nal.pl binary file. Default tries if pal2nal.pl is in PATH;

  - __[--snp-sites-bin (PATH)]__         Path to the snp-sites binary file. Default tries if snp-sites is in PATH;

  - __[--roary-bin (PATH)]__             Path to the roary binary file. Default tries if roary is in PATH;

  - __[--orthofinder-bin (PATH)]__       Path to the orthofinder binary file. Default tries if orthofinder is in PATH;

  - __[--fastANI-bin (PATH)]__           Path to the fastANI binary file. Default tries if fastANI is in PATH;

  - __[--gubbins-bin (PATH)]__           Path to the run_gubbins.py binary file. Default tries if run_gubbins.py is in PATH;

  - __[--snippy-bin (PATH)]__            Path to the snippy binary file. Default tries if snippy is in PATH;

  - __[--sickle-bin (PATH)]__            Path to the sickle-trim binary file. Default tries if sickle is in PATH.
<br/>

- __Setup COG database__

  - __[--setup-COGdb]__                  Users should execute this after first installation of pgcgap.
<br/>

- Check the required external programs (__It is strongly recommended that this step be performed after the installation of pcgp__):
    ```
    $pgcgap --check-external-programs
    ```
<br/>

- __EXAMPLES:__
  - __Example 1:__ Perform all functions, take the *bacillus thuringiensis* as an example, total 4 strains for analysis.<br/>

    __Notice__: For the sake of flexibility, The "VAR" function needs to be added additionally.<br/>

    ```
    $pgcgap --All --ReadsPath <PATH> --reads1 .R1.clean.fastq.gz --reads2 .R2.clean.fastq.gz --suffix_len <INT> --kmmer 81 --genus bacillus --species thuringiensis --codon 11 --strain_num 4 --threads 4 --VAR --refgbk <FILE> --qualtype <STRING>
    ```

  - __Example 2:__ Conduct reads assembly, gene predition and annotation.

    ```
    $pgcgap --Assemble --ReadsPath <PATH> --reads1 .R1.clean.fastq.gz --reads2 .R2.clean.fastq.gz --kmmer 81 --genus bacillus --species thuringiensis --codon 11 --threads 4
    ```

  - __Example 3:__ Constructing single core protein tree and core SNPs tree.

    ```
    $pgcgap --CoreTree --CDsPath <PATH> --AAsPath <PATH> --codon 11 --strain_num 3 --threads 4
    ```

  - __Example 4:__ Conduct pan-genome analysis.

    ```
    $pgcgap --Pan --codon 11 --threads 4 --GffPath <PATH>
    ```

  - __Example 5:__ Inference of orthologous gene groups.

    ```
    $pgcgap --orthoF --threads 4 --AAsPath <PATH>
    ```

  - __Example 6:__ Compute whole-genome Average Nucleotide Identity (ANI).

    ```
    $pgcgap --ANI --threads 4 --queryL <FILE> --refL <FILE> --ANIO <FILE>
    ```

  - __Example 7:__ Run COG annotation for each strain.

    ```
    $pgcgap --COG --threads 4 --AAsPath <PATH>
    ```

  - __Example 8:__ Varients calling and phylogenetic tree construction based on reference genome.

    ```
    $pgcgap --VAR --threads 4 --refgbk <FILE> --ReadsPath <PATH> --reads1 .R1.clean.fastq.gz --reads2 .R2.clean.fastq.gz --suffix_len <INT> --strain_num <INT> --qualtype <STRING>
    ```

## Generating Input files

### Working directory
The directory where the PGCGAP software runs.

### Assemble
Pair-end reads of all strains in a directory (default: ./Reads/ under the working directory).

### ANI

QUERY_LIST and REFERENCE_LIST files containing full paths to genomes, one per line (default: scaf.list under the working directory). If the "--Assemble" fuction was run first, the list file will be generated automatically.

### CoreTree
Amino acids file (With ".faa" as the suffix) and nucleotide (With ".ffn" as the suffix) file of each strain placed into two directorys (default: ./Results/Annotations/AAs/ and ./Results/Annotations/CDs/). If the "--Assemble" fuction was run first, the files will be generated automatically.

### OrthoF
A set of protein sequence files (one per species) in FASTA format under a directory (default: ./Results/Annotations/AAs/). If the "--Assemble" fuction was run first, the files will be generated automatically.

### Pan
GFF3 files (With ".gff" as the suffix) of each strain placed into a directory. They must contain the nucleotide sequence at the end of the file. All GFF3 files created by Prokka are valid (default: ./Results/Annotations/GFF/). If the "--Assemble" fuction was run first, the files will be generated automatically.

### COG
Amino acids file (With ".faa" as the suffix) of each strain placed into a directory (default: ./Results/Annotations/AAs/). If the "--Assemble" fuction was run first, the files will be generated automatically.

### VAR
- Pair-end reads of all strains in a directory (default: ./Reads/Over/ under the working directory).
<br/>

- The full path of reference genome in fasta format or genebank format (__must be provided__).

## Output Files

### Assemble

- **Results/Assembles/*_assembly**<br/>
Directorys contain assembly files and information of each strain.
<br/>

- __Results/Assembles/Scaf__<br/>
Directory contain contig/scaffold of all strains.
<br/>

- **Results/Annotations/*_annotation**<br/>
Directorys contain [annotation files](https://github.com/tseemann/prokka?_blank) of each strain.
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
query genome, reference genome, ANI value, count of bidirectional fragment mappings, total query fragments
<br/>

- __Results/ANI/ANIs.matrix__<br/>
file with identity values arranged in a [phylip-formatted lower triangular matrix](https://www.mothur.org/wiki/Phylip-formatted_distance_matrix?_blank)
<br/>

- __Results/ANI/ANIs.heatmap__<br/>
An ANI matrix of all strains
<br/>

- __Results/ANI/ANI_matrix.pdf__<br/>
The heatmap plot of "ANIs.heatmap"


### CoreTree

- __Results/CoreTrees/faa/ALL.core.protein.fast__<br/>
Concatenated and aligned sequences file of single-core proteins.
<br/>

- __Results/CoreTrees/faa2ffn/ALL.core.nucl.fasta__<br/>
Concatenated and aligned sequences file of single-core genes.
<br/>

- __Results/CoreTrees/faa2ffn/ALL.core.snp.fasta__<br/>
Core SNPs of single-core genes in fasta format.
<br/>

- __Results/CoreTrees/ALL.core.protein.nwk__<br/>
The phylogenetic tree file of single-copy proteins for all strains.
<br/>

- __Results/CoreTrees/faa2ffn/ALL.core.snp.nwk__ (or gubbins.core.snp.final_tree.tre)<br/>
The phylogenetic tree file of SNPs of single-copy genes for all strains.
<br/>

- __Results/CoreTrees/faa2ffn/gubbins.*__<br/>
See [gubbins](https://github.com/sanger-pathogens/gubbins#output-files?_blank) outputs
<br/>

- __Results/CoreTrees/"Other_files"__<br/>
Intermediate directorys and files.
<br/>

### OrthoF

- __Results/OrthoF/Results_orthoF__<br/>
Same as [OrthoFinder](https://github.com/davidemms/OrthoFinder?_blank) outputs
<br/>

### Pan

- __Results/PanGenome/Pangenome_Pie.pdf__<br/>
An 3D pie chart of the breakdown of genes and the number of isolate they are present in
<br/>

- __Results/PanGenome/pangenome_frequency.pdf__<br/>
A graph with the frequency of genes versus the number of genomes
<br/>

- __Results/PanGenome/Pangenome_matrix.pdf__<br/>
A figure showing the tree compared to a matrix with the presence and absence of core and accessory genes
<br/>

- __Results/PanGenome/Other_files__<br/>
see [roary](https://sanger-pathogens.github.io/Roary/?_blank) outputs
<br/>

### COG

- __*.COG.xml, *.2gi.table, *.2id.table, *.2Sid.table__<br/>
Intermediate files
<br/>

- __*.2Scog.table__<br/>
The super COG table of each strain
<br/>

- __*.2Scog.table.pdf__<br/>
A plot of super COG table in pdf format
<br/>

- __All_flags_relative_abundances.table__
A table containing the relative abundance of each flag for all strains

### VAR
- __Results/Varients/directory-named-in-strains__<br/>
Directorys containing substitutions (snps) and insertions/deletions (indels) of each strain. See [Snippy](https://github.com/tseemann/snippy?_blank) outputs for detail.

- __Results/Varients/Core__<br/>
  The directory containing Core SNP phylogeny files

  - __.aln__ : A core SNP alignment in the --aformat format (default FASTA)
  - __.full.aln__ : A whole genome SNP alignment (includes invariant sites)
  - __.nwk__ : Phylogenetic tree constructed with __FastTree__ (ignoring possible recombination)
  - **_tree.tre** : Phylogenetic tree constructed with __gubbins__ (get rid of recombination)

## License

PGCGAP is free software, licensed under GPLv3.

## Feedback and Issues

Please report any issues to the [issues page](https://github.com/liaochenlanruo/pgcgap/issues?_blank) or email us at [liaochenlanruo@webmail.hzau.edu.cn](mailto:liaochenlanruo@webmail.hzau.edu.cn).

## Citation

If you use this software please cite:

## FAQ

### Q1 VAR FOUNCTION ran failed to get annotated VCFs and Core results

Check the log file named in "strain_name.log" under Results/Varients/<strain_name>/ directory. If you find a sentence like "WARNING: All frames are zero! This seems rather odd, please check that 'frame' information in your 'genes' file is accurate." This is an snpEFF error. Users can install JDK8 to solve this problem. 

```
$conda install java-jdk=8.0.112
```

Click [here](https://github.com/tseemann/snippy/issues/259?_blank) for more solutions.

------------------------------------------------------------------------

<p><center><strong>
<script async src="//busuanzi.ibruce.info/busuanzi/2.3/busuanzi.pure.mini.js"></script>
<span id="busuanzi_container_site_pv">Total visits: <span id="busuanzi_value_site_pv"></span> times</span>
<span class="post-meta-divider">|</span>
<span id="busuanzi_container_site_uv">Visitors: <span id="busuanzi_value_site_uv"></span> people</span>
</strong></center></p>
