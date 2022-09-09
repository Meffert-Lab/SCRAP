# SCRAP: a bioinformatic pipeline for the analysis of small chimeric RNA-seq data

| File Name &nbsp;                    | Description |
| -------------- | ---------- |
| SCRAP_environment.yml   | File for creating a Conda environment with the requisite tools for running SCRAP        |
| Reference_Installation.sh | Script for configuring reference files required for SCRAP |
| SCRAP.sh   | Script for processing raw FASTQ files to identify sncRNA and genomic alignment of reads        |
| Peak_Calling.sh   | Script for calling peaks using output files from SCRAP.sh        |
| Peak_Annotation.sh   | Script for annotating bed file produced by Peak_Calling.sh with gene names and features        |
| miRBase.fasta      | FASTA file containing miRNAs downloaded from [miRBase](https://www.mirbase.org/) (accessed July 15, 2022)        |
| miRBase.hairpin.fasta      | FASTA file containing miRNA hairpin sequences obtained from [miRBase](https://www.mirbase.org/) (accessed July 15, 2022) |
| GtRNAdb.fasta      | FASTA file containing tRNA sequences obtained from [GtRNAdb](http://gtrnadb.ucsc.edu) (accessed July 15, 2022)        |
| miRNA family file      | Add description        |


## Installation

Once in the directory where you would like the PACeR files (XX MB) to be installed, run:

    sudo apt install git
    git clone https://github.com/Meffert-Lab/SCRAP.git

Create the Conda environment by running (requires [Miniconda](https://docs.conda.io/en/latest/miniconda.html)):

    conda install -n base conda-forge::mamba
    mamba env create -f SCRAP/SCRAP_environment.yml -n SCRAP

Execute the `reference_installation.sh` script with the following command line parameters:

| Command Line Parameter | Description |
| ---------- | ---------- |
| `-r` | Path to reference directory (e.g. `SCRAP`) |
| `-m` | Three-letter miRBase species abbreviation |
| `-g` | Reference genome abbreviation |

Example code for configuring human references:

    bash SCRAP/Reference_Installation.sh \
        -r SCRAP/ \
        -m hsa \
        -g hg38

## Running PACeR

Ensure data files are in the following configuration

        │───SCRAP
        │         SCRAP_environment.yml
        │         Reference_Installation.sh
        │         SCRAP.sh
        │         Peak_Calling.sh
        │         Peak_Annotation.sh
        │         miRBase.fasta
        │         miRBase.hairpin.fasta
        │         GtRNAdb.fasta
        │  
        └───files 
             │       
             │
             └───sample1
             │       sample1_R1.fastq.gz
             │       sample1_R2.fastq.gz
             │
             └───sample2
             │       sample2_R1.fastq.gz
             │       sample2_R2.fastq.gz
             │
             └───sample3
                     sample3_R1.fastq.gz
                     sample3_R2.fastq.gz

Execute the `SCRAP.sh` script with the following command line parameters:

| Command Line Parameter | Description |
| ---------- | ---------- |
| `-d` | Path to directory containing sample directories |
| `-a` | Path to adapter file |
| `-p` | Denote wether samples are paired-end (`yes` or `no`) |
| `-f` | Indicate whether or not to filter out pre-miRNAs and tRNAs (`yes` or `no`) |
| `-r` | Path to reference directory (e.g. `SCRAP`) |
| `-m` | Three-letter miRBase species abbreviation |
| `-g` | Reference genome abbreviation |

Example code for analyzing CLASH data:

    bash PACeR.sh \
        -d CLASH_Human/ \
        -a CLASH_Human/CLASH_Human_Adapters.txt \
        -p no \
        -f yes \
        -r SCRAP/ \
        -m hsa \
        -g hg38
