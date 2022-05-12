# crisprDesign

## Table of content
- [1. Overview](#id-section1)
- [2. Software requirements](#id-section1)
- [3. Installation](#id-section2)
- [4. Demo](#id-section3)

<div id='id-section1'/>


## Overview


Functions for design and annotation of CRISPR single-guide RNAs.


All RNA- and DNA-targeting nucleases are supported, including SpCas9, 
AsCas12a, enAsCas12, RfxCas13d, etc. Base editors are also supported,
such as BE4max. 

The package contains advanced functionalities for CRISPR knockout (CRISPRko),
CRISPR activation (CRISPRa), CRISPR inhibition (CRISPRi) and
CRISPR base editing (CRISPRbe). 

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

