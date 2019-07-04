# PGCGAP - the Prokaryotic Genomics and Comparative Genomics Analysis Pipeline
![Platform](https://badgen.net/badge/platform/Linux,macOS?list=|)
![License](https://badgen.net/github/license/liaochenlanruo/pgcgap)
[![GitHubversion](https://badge.fury.io/gh/liaochenlanruo%2Fpgcgap.svg)](https://badge.fury.io/gh/liaochenlanruo%2Fpgcgap)
![Downloads conda](https://img.shields.io/conda/dn/bioconda/pgcgap.svg?style=flat)
![Watchers](https://badgen.net/github/watchers/liaochenlanruo/pgcgap)
![Contributors](https://badgen.net/github/contributors/liaochenlanruo/pgcgap)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pgcgap/README.html)


[English Readme](https://liaochenlanruo.github.io/pgcgap?_blank) | [Chinese Readme](https://liaochenlanruo.github.io/2019/04/28/PGCGAP%E4%B8%AD%E6%96%87%E8%AF%B4%E6%98%8E/?_blank)

![earth](<script type="text/javascript" src="//rf.revolvermaps.com/0/0/8.js?i=5r6yr3vsmt5&amp;m=6&amp;c=ff0000&amp;cr1=ffffff&amp;f=arial&amp;l=33" async="async"></script>)

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Required dependencies](#required-dependencies)
- [License](#license)
- [Feedback and Issues](#feedback-and-issues)
- [Citation](#citation)

- [Usages](#usages)

## Introduction

PGCGAP is a pipeline for prokaryotic comparative genomics analysis. It can take the pair-end reads as input. In addition to genome assembly, gene prediction and annotation, it can also get common comparative genomics analysis results such as phylogenetic trees of single-core proteins and core SNPs, pan-genome, whole-genome Average Nucleotide Identity (ANI), orthogroups and orthologs, COG annotations, substitutions (snps) and insertions/deletions (indels) with only one line of commands.

## Installation

The software was tested successfully on Linux x64 platform and macOS. However, Windows could not be supported. Because this software relies on a large number of other softwares, so it is recommended to install with __[Bioconda](https://bioconda.github.io/index.html)__. The main program and most of other dependencies can be installed with one command as shown below, __but the "Gubbins" should be installed separately__ because of it relys on python3, while PGCGAP relys on python2.


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
- [Coreutils](https://www.gnu.org/software/coreutils/)
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
- [wget](https://www.gnu.org/software/wget/)


## License

PCGAP is free software, licensed under GPLv3.

## Feedback and Issues

Please report any issues to the [issues page](https://github.com/liaochenlanruo/pcgap/issues) or email us at [liaochenlanruo@webmail.hzau.edu.cn](mailto:liaochenlanruo@webmail.hzau.edu.cn).

## Citation

If you use this software please cite: (__Please keep an eye on it as it will be noted soon!__)


## Usages
For more detial informations, please visit the webpage of [PGCGAP](https://liaochenlanruo.github.io/pgcgap?_blank).

