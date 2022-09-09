# SCRAP: a bioinformatic pipeline for the analysis of small chimeric RNA-seq data

| File Name &nbsp;                    | Description |
| -------------- | ---------- |
| environment.yml   | YML file for creating Conda environment with packages required for running PACeR.        |
| hisat2_index_downloads.txt | Text file containing links to HISAT2-indexed genomes of several common model organisms. |
| miRBase.fasta      | FASTA file containing miRNAs downloaded from [miRBase](https://www.mirbase.org/) (accessed July 15, 2022).        |
| reference_installation.sh      | Shell script for downloading and configuring reference files.        |
| PACeR.sh      | Shell script for running PACeR.        |
| PACeR_CLEAR-CLIP.sh      | Modified shell script for running PACeR with data from [Moore et al. 2015](https://www.nature.com/articles/ncomms9864).        |
| PeakCalling_Total.sh      | Shell script for calling peaks for all sncRNAs from output of `PACeR.sh`.        |
| PeakCalling_Subset.sh      | Shell script for calling peaks for specific sncRNAs or sncRNA familes from output of `PACeR.sh`.        |

## Installation

Once in the directory where you would like the PACeR files (XX MB) to be installed, run:

    sudo apt install git
    git clone https://github.com/Meffert-Lab/PACeR_New.git

Create the Conda environment by running:

    conda config --append channels conda-forge
    conda env create -f PACeR_New/environment.yml -n PACeR

Execute the `reference_installation.sh` script with the following command line parameters (requires [Miniconda](https://docs.conda.io/en/latest/miniconda.html)):

| Command Line Parameter | Description |
| ---------- | ---------- |
| `-r` | Path to reference directory (e.g. `PACeR_New`) |
| `-m` | Three-letter miRBase species abbreviation |
| `-g` | Reference genome abbreviation |

Example code for configuring human references:

    bash PACeR_New/reference_installation.sh \
        -r /media/sf_Ubuntu_Sharing_2022/PACeR_New \
        -m hsa \
        -g hg38

## Running PACeR

Ensure data files are in the following configuration

        │───PACeR_New
        │         environment.yml
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
