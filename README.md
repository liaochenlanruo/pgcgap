# PGCGAP - the Prokaryotic Genomics and Comparative Genomics Analysis Pipeline
![](https://img.shields.io/badge/Platform-WSL%2FLinux%2FmacOS-green)
[![GitHub license](https://img.shields.io/github/license/liaochenlanruo/pgcgap)](https://github.com/liaochenlanruo/pgcgap/blob/master/LICENSE)
[![Version](https://anaconda.org/bioconda/pgcgap/badges/version.svg)](https://anaconda.org/bioconda/pgcgap)
![Downloads conda](https://img.shields.io/conda/dn/bioconda/pgcgap.svg?style=flat)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pgcgap/README.html)
<a target="_blank" href="https://qm.qq.com/cgi-bin/qm/qr?k=f6SeZ8F9dDTaZwLZZBymkczRXUHLW6Q1&jump_from=webapi"><img border="0" src="https://pub.idqqimg.com/wpa/images/group.png" alt="生信之巅" title="生信之巅"></a>


[English Readme](https://liaochenlanruo.github.io/pgcgap?_blank) | [Chinese Readme](https://liaochenlanruo.github.io/2019/04/28/PGCGAP%E4%B8%AD%E6%96%87%E8%AF%B4%E6%98%8E/?_blank)


		  ____       ____      ____     ____       _        ____    
		U|  _"\ u U /"___|u U /"___| U /"___|u U  /"\  u  U|  _"\ u 
		\| |_) |/ \| |  _ / \| | u   \| |  _ /  \/ _ \/   \| |_) |/ 
		 |  __/    | |_| |   | |/__   | |_| |   / ___ \    |  __/   
		 |_|        \____|    \____|   \____|  /_/   \_\   |_|      
		 ||>>_      _)(|_    _// \\    _)(|_    \\    >>   ||>>_    
		(__)__)    (__)__)  (__)(__)  (__)__)  (__)  (__) (__)__)   


## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Required dependencies](#required-dependencies)
- [License](#license)
- [Feedback and Issues](#feedback-and-issues)
- [Citation](#citation)

- [Usages](#usages)

## Introduction

PGCGAP is a pipeline for prokaryotic comparative genomics analysis. It can take the pair-end reads as input. In addition to genome assembly, gene prediction and annotation, it can also get common comparative genomics analysis results such as phylogenetic trees of single-core proteins and core SNPs, pan-genome, whole-genome Average Nucleotide Identity (ANI), orthogroups and orthologs, COG annotations, substitutions (snps) and insertions/deletions (indels) and antimicrobial and virulence genes mining with only one line of commands.

## Installation

The software was tested successfully on Windows WSL, Linux x64 platform and macOS. Because this software relies on a large number of other softwares, so it is recommended to install with __[Bioconda](https://bioconda.github.io/index.html)__.


__Step1: Install PGCGAP__

```
$conda create -n pgcgap python=3
$conda activate pgcgap
$conda install pgcgap
```
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

__Step2: Setup COG database__ (Users should execute this after first installation of pgcgap)

```
$conda activate pgcgap
$pgcgap --setup-COGdb
$conda deactivate
```


Users with [docker container](https://hub.docker.com/) installed have another choice to install PGCGAP.

```
$docker pull quay.io/biocontainers/pgcgap:<tag>
```

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
- [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks_server.html)
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
- [trimAL](http://trimal.cgenomics.org/)
- [unicycler](https://github.com/rrwick/Unicycler)
- [wget](https://www.gnu.org/software/wget/)


## License

PGCGAP is free software, licensed under GPLv3.

## Feedback and Issues

Please report any issues to the [issues page](https://github.com/liaochenlanruo/pcgap/issues) or email us at [liaochenlanruo@webmail.hzau.edu.cn](mailto:liaochenlanruo@webmail.hzau.edu.cn).

## Citation

If you use this software please cite: Liu H, Xin B, Zheng J, Zhong H, Yu Y, Peng D, Sun M. Build a bioinformatics analysis platform and apply it to routine analysis of microbial genomics and comparative genomics. *Protocol exchange*, 2021. DOI: [10.21203/rs.2.21224/v5](https://doi.org/10.21203/rs.2.21224/v5)


## Usages


For more detial informations, please visit the webpage of [PGCGAP](https://liaochenlanruo.github.io/pgcgap?_blank) and [WIKI](https://github.com/liaochenlanruo/pgcgap/wiki).

