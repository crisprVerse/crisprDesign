# crisprDesign

## Table of content
- [1. Overview](#id-section1)
- [2. Software requirements](#id-section1)
- [3. Installation](#id-section2)
- [4. Demo](#id-section3)

<div id='id-section1'/>


## Overview


`crisprDesign` is a one-stop shop R package for CRISPR gRNA design. It provides a comprehensive suite of functions to design and annotate CRISPR guide RNA (gRNAs) sequences. This includes on- and off-target search, on-target efficiency scoring, off-target scoring, full gene and TSS
contextual annotations, and SNP annotation (human only). It currently support five types of CRISPR modalities (modes of perturbations): CRISPR knockout (CRISPRko), CRISPR activation (CRISPRa), CRISPR inhibition (CRISPRi), CRISPR base editing (CRISPRbe), and CRISPR knockdown (CRISPRkd). All types of CRISPR nucleases are supported, including DNA- and RNA-target nucleases such as Cas9, Cas12a, and Cas13d. All types of base editors are also supported. gRNA design can be performed on reference genomes, transcriptomes, and custom DNA and RNA sequences. 

Our work is described in a recent bioRxiv preprint: ["A comprehensive Bioconductor ecosystem for the design of CRISPR guide RNAs across nucleases and technologies"](https://www.biorxiv.org/content/10.1101/2022.04.21.488824v2)


<div id='id-section2'/>

## Software requirements

### OS Requirements

This package is supported for macOS, Linux and Windows machines.
Some functionalities are not supported for Windows machines.

Packages were developed and tested on R version 4.2.

### R Dependencies 

- crisprBase: https://github.com/Jfortin1/crisprBase

- crisprBowtie: https://github.com/Jfortin1/crisprBowtie

- crisprScoreData: https://github.com/Jfortin1/crisprScoreData

- crisprScore: https://github.com/Jfortin1/crisprScore


### Optional dependencies 

The following dependencies are needed to run BWA. Note that they are not available for Windows machines.

- Rbwa: https://github.com/Jfortin1/Rbwa

- crisprBwa: https://github.com/Jfortin1/crisprBwa
 
<div id='id-section3'/>

## Installation

`crisprDesign` and its dependencies can be installed by typing the following commands inside of an R session:

```r
install.packages("devtools")
library(devtools)
install_github("Jfortin1/crisprBase")
install_github("Jfortin1/crisprBowtie")
install_github("Jfortin1/crisprScore")
install_github("Jfortin1/crisprScoreData")
install_github("jfortin1/crisprDesign")
install_github("jfortin1/crisprDesignData")
```


<div id='id-section4'/>

## Demo 

A reproducible and comprehensive workflow can be found here:
https://github.com/Jfortin1/crisprDesign/blob/master/vignettes/intro.Rmd


## Creating bowtie and BWA indexes

Downloading fasta for human genome (GRCh38):

```{bash}
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
gunzip hg38.fa.gz
```

Creating bowtie index, assuming the bowtie binary is on the PATH, 
and naming the bowtie index hg38:

```{bash}
bowtie-build hg38.fa hg38
```

Creating a BWA index, assuming the BWA binary is on the PATH, 
and naming the BWA index hg38:

```{}
bwa index -p ./hg38 hg38.fa
```

## License 

This project is covered under the MIT License.

