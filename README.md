# SCRAP: a bioinformatic pipeline for the analysis of small chimeric RNA-seq data

| File Name &nbsp;                    | Description |
| -------------- | ---------- |
| SCRAP_environment.yml   | File for creating a Conda environment with the requisite tools for running SCRAP        |
| Reference_Installation.sh | Script for configuring reference files required for SCRAP |
| SCRAP.sh   | Script for processing compressed FASTQ files to generate a BAM file containing uniquely aligned reads annotated with the corresponding sncRNA within the read        |
| Peak_Calling.sh   | Script for calling peaks using output files from SCRAP.sh        |
| Peak_Annotation.sh   | Script for annotating bed file produced by Peak_Calling.sh with gene names and features        |
| miRBase.fasta      | FASTA file containing miRNAs downloaded from [miRBase](https://www.mirbase.org/) (accessed July 15, 2022)        |
| miRBase.hairpin.fasta      | FASTA file containing miRNA hairpin sequences obtained from [miRBase](https://www.mirbase.org/) (accessed July 15, 2022) |
| GtRNAdb.fasta      | FASTA file containing tRNA sequences obtained from [GtRNAdb](http://gtrnadb.ucsc.edu) (accessed July 15, 2022)        |
| miRNA family file      | Add description        |


## Installation

Once in the directory where you would like the PACeR files (XX MB) to be installed, run:

    sudo apt install git
    git clone https://github.com/Meffert-Lab/PACeR_New.git

Create the Conda environment by running:

    conda install -n base conda-forge::mamba
    mamba env create -f SCRAP/SCRAP_environment.yml -n SCRAP

Execute the `reference_installation.sh` script with the following command line parameters (requires [Miniconda](https://docs.conda.io/en/latest/miniconda.html)):

| Command Line Parameter | Description |
| ---------- | ---------- |
| `-r` | Path to reference directory (e.g. `SCRAP`) |
| `-m` | Three-letter miRBase species abbreviation |
| `-g` | Reference genome abbreviation |

Example code for configuring human references:

    bash SCRAP/reference_installation.sh \
        -r SCRAP/ \
        -m hsa \
        -g hg38

## Running PACeR

Ensure data files are in the following configuration

        │───SCRAP
        │         SCRAP_environment.yml
        │         hisat2_index_downloads.txt
        │         PACeR.sh
        │         PeakCalling.sh
        │         reference_installation.sh
        │         sncRNA.fasta
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

Execute the `PACeR.sh` script with the following command line parameters:

| Command Line Parameter | Description |
| ---------- | ---------- |
| `-d` | Path to directory containing sample directories |
| `-a` | Path to adapter file |
| `-p` | Denote wether samples are paired-end (`yes` or `no`) |
| `-r` | Path to reference directory (e.g. `PACeR_New`) |
| `-m` | Three-letter miRBase species abbreviation |
| `-g` | Reference genome abbreviation |

Example code for analyzing CLASH data:

    bash PACeR.sh \
        -d /media/sf_Ubuntu_Sharing_2022/CLASH_Human \
        -a /media/sf_Ubuntu_Sharing_2022/CLASH_Human/CLASH_Human_Adapters.txt \
        -p no \
        -r /media/sf_Ubuntu_Sharing_2022/PACeR_New \
        -m hsa \
        -g hg38
