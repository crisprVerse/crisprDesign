# crisprDesign

Functions for design and annotation of CRISPR single-guide RNAs.


All RNA- and DNA-targeting nucleases are supported, including SpCas9, 
AsCas12a, enAsCas12, RfxCas13d, etc. Base editors are also supported,
such as BE4max. 

The package contains advanced functionalities for CRISPR knockout (CRISPRko),
CRISPR activation (CRISPRa), CRISPR inhibition (CRISPRi) and
CRISPR base editing (CRISPRbe). 

### Dependencies 

This is depending on the following R packages:

- crisprBase: https://github.com/Jfortin1/crisprBase

- crisprBowtie: https://github.com/Jfortin1/crisprBowtie

- crisprBwa: https://github.com/Jfortin1/crisprBwa

- crisprScore: https://github.com/Jfortin1/crisprScore

- Rbwa: https://github.com/Jfortin1/Rbwa

- Cas13design (optional): https://github.com/Jfortin1/Cas13design

### Example data

- crisprDesignData: https://github.com/Jfortin1/crisprDesignData

### Creating bowtie and BWA indexes

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




