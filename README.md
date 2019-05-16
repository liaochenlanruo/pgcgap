# PGCGAP - the Prokaryotic Genomics and Comparative Genomics Analysis Pipeline

[English Readme](https://liaochenlanruo.github.io/pgcgap?_blank) | [Chinese Readme](https://liaochenlanruo.github.io/2019/04/28/PGCGAP%E4%B8%AD%E6%96%87%E8%AF%B4%E6%98%8E/?_blank)

## Contents

- [PGCGAP - The prokaryotic genomics and comparative genomics analysis pipeline](#PGCGAP---the-Prokaryotic-Genomics-and-Comparative-Genomics-Analysis-Pipeline)
    - [Introduction](#Introduction)
    - [Installation](#Installation)
    - [Required dependencies](#Required-dependencies)
    - [License](#License)
    - [Feedback/Issues](#Feedback/Issues)
    - [Citation](#Citation)
    - [Usages](#Usages)

## Introduction

PGCGAP is a pipeline for prokaryotic comparative genomics analysis. It can take the pair-end reads as input. In addition to genome assembly, gene prediction and annotation, it can also get common comparative genomics analysis results such as phylogenetic trees of single-core proteins and core SNPs, pan-genome, whole-genome Average Nucleotide Identity (ANI), orthogroups and orthologs, COG annotations, substitutions (snps) and insertions/deletions (indels) with only one line of commands.

## Installation

The software was only tested on Linux64 platform, The MacOS may also be good in theory. However, Windows could not be supported. Because this software relies on a large number of other softwares, so it is recommended to install with __[Bioconda](https://bioconda.github.io/index.html)__. The main program and most of other dependencies can be installed with one command as shown below, __but the "pan2nal" and "OrthoFinder" should be installed separately!__

- Bioconda - OSX/Linux
    ```
    $conda install pgcgap
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

## License

PCGAP is free software, licensed under GPLv3.

## Feedback/Issues

Please report any issues to the [issues page](https://github.com/liaochenlanruo/pcgap/issues) or email us at [liaochenlanruo@webmail.hzau.edu.cn](mailto:liaochenlanruo@webmail.hzau.edu.cn).

## Citation

If you use this software please cite:

## Usages
For more detial informations, please visit the webpage of [PGCGAP](https://liaochenlanruo.github.io/pgcgap?_blank).
