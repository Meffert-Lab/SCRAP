# SCRAP: a bioinformatic pipeline for the analysis of small chimeric RNA-seq data

## Contents
| File or Directory Name &nbsp;                    | Description |
| -------------- | ---------- |
| adatpers/	| Contains adapter sequences for CLASH and CLEAR-CLIP	|
| annotation/	| Contains annotation files for human, mouse, and c elegans, as well as miRNA family file.	|
| bin/	| Contains SCRAP scripts	|
| fasta/	| contains miRBase and tRNA FASTAs	|
| PLATFORM-SETUP.md     | Preliminary setup instructions for each compatible platform   |
| README.md	| Contains instructions to install and run SCRAP	|
| SCRAP_environment.yml   | File for creating a Conda environment with the requisite tools for running SCRAP        |

### adapters
| File Name &nbsp;                    | Description |
| -------------- | ---------- |
| CLASH_Human_Adapters.txt	| Adapters used in the CLASH libraries	|
| CLEAR-CLIP_Mouse_Adapters.txt	| Adapters used in the CLEAR-CLIP libraries	|

### annotation
| File Name &nbsp;                    | Description |
| -------------- | ---------- |
| miR_Family.txt	| Tab-delimited list detailing which miRNA families contain which miRBase miRNAs        |
| mouse.annotation.bed.gz	| Mouse genome annotation file used to annotate peaks after calling	|
| human.annotation.bed.gz	| Human genome annotation file used to annotate peals after calling	|
| worm.annotation.bed.gz	| C. elegans genome annotation file used to annotate peaks after calling	|

### bin
| File Name &nbsp;                    | Description |
| -------------- | ---------- |
| Reference_Installation.sh | Script for configuring reference files required for SCRAP |
| SCRAP.sh   | Script for processing raw FASTQ files to identify sncRNA and genomic alignment of reads        |
| Peak_Calling.sh   | Script for calling peaks using output files from SCRAP.sh        |
| Peak_Annotation.sh   | Script for annotating bed file produced by Peak_Calling.sh with gene names and features        |

### fasta
| File Name &nbsp;                    | Description |
| -------------- | ---------- |
| miRBase.fasta      | FASTA file containing miRNAs downloaded from [miRBase](https://www.mirbase.org/) (accessed July 15, 2022)        |
| miRBase.hairpin.fasta      | FASTA file containing miRNA hairpin sequences obtained from [miRBase](https://www.mirbase.org/) (accessed July 15, 2022) |
| GtRNAdb.fasta      | FASTA file containing tRNA sequences obtained from [GtRNAdb](http://gtrnadb.ucsc.edu) (accessed July 15, 2022)        |
| tRFdb.fasta	| FASTA file containing tRNA fragment sequences obtained from [miRBase](https://www.mirbase.org) (accessed July 15, 2022)	|

## Installation

SCRAP is a cross-platform pipeline that can be used in Windows Subsystem for Linux, MacOS, and Ubuntu.
In order to install SCRAP, you need one of these platforms as well as Git and Miniconda.
See PLATFORM-SETUP.md for specific instructions.

Once in the directory where you would like the SCRAP source to be cloned, run:

    git clone https://github.com/Meffert-Lab/SCRAP.git

Create the Conda environment by running (requires [Miniconda](https://docs.conda.io/en/latest/miniconda.html), see PLATFORM-SETUP.md):

    conda install -n base conda-forge::mamba
    mamba env create -f SCRAP/SCRAP_environment.yml -n SCRAP

Execute the `Reference_Installation.sh` script with the following command line parameters:

| Flag | Description |
| ---------- | ---------- |
| `-r` | Path to reference directory (e.g. `SCRAP`) |
| `-m` | Three-letter miRBase species abbreviation |
| `-g` | Reference genome abbreviation |
| `-s` | Indicate species used for annotation (`human` (*H. sapiens*), `mouse` (*M. musculus*), or `worm` (*C. elegans*) |


### Three-letter miRBase Species Abbreviations
| Abbreviation	| Species	|
| hsa	| Homo sapiens	|
| mmu	| Mus musculus	|
| rno	| Rattus norvegicus	|
| dme	| Drosophila melanogaster	|
| cel	| Caenorhabditis elegans	|
| ath	| Arabidopsis thaliana	|

Example code for configuring human references:

    bash SCRAP/bin/Reference_Installation.sh \
        -r SCRAP/ \
        -m hsa \
        -g hg38 \
        -s human

## Running SCRAP

Ensure data files are in the following configuration

	│───SCRAP
	│	│
	│	│
	│	└───bin        
	│	│	SCRAP.sh
	│	│	Peak_Calling.sh
	│	│	Peak_Annotation.sh
	│	│	Reference_Installation.sh
	│	│
	│	└───fasta
	│	│	miRBase.fasta
	│	│	miRBase.hairpin.fasta
	│	│	GtRNAdb.fasta
	│	│	tRFdb.fasta
	│	└───annotation
	│		human.annotation.bed
	│		miR_Family.txt
	│		mouse.annotation.bed
	│		worm.annotation.bed
	│
	│
	│
	└───files 
		│
		│
		└───sample1
		│	sample1_R1.fastq.gz
		│	sample1_R2.fastq.gz
		│
		└───sample2
		│	sample2_R1.fastq.gz
		│	sample2_R2.fastq.gz
		│
		└───sample3
			sample3_R1.fastq.gz
			sample3_R2.fastq.gz

Execute the `SCRAP.sh` script with the following command line parameters:

| Flag | Description |
| ---------- | ---------- |
| `-d` | Path to directory containing sample directories |
| `-a` | Path to adapter file |
| `-p` | Denote wether samples are paired-end (`yes` or `no`) |
| `-f` | Indicate whether or not to filter out pre-miRNAs and tRNAs (`yes` or `no`) |
| `-r` | Path to reference directory (e.g. `SCRAP`) |
| `-m` | Three-letter miRBase species abbreviation |
| `-g` | Reference genome abbreviation |

The adapter file is a tab-delimited `.txt` file containing the sample name, 5' adapter, 3' adapter, 5' barcode, and 3' barcode. This file can be generated with a text editor or in Excel and [saved as a tab-delimited text file](https://support.microsoft.com/en-us/office/save-a-workbook-to-text-format-txt-or-csv-3e9a9d6c-70da-4255-aa28-fcacf1f081e6).

Example code for analyzing CLASH data:

    bash SCRAP.sh \
        -d CLASH_Human/ \
        -a CLASH_Human/CLASH_Human_Adapters.txt \
        -p no \
        -f yes \
        -r SCRAP/ \
        -m hsa \
        -g hg38

After data have been analyzed with `SCRAP.sh`, sample folders will contain a file ending in `.aligned.unique.sam`

## Peak Calling

The `.aligned.unique.sam` file produced by `SCRAP.sh` can be used to identify peaks where multiple sncRNAs or sncRNA familiy members bind to the same region of the genome.

Execute the `Peak_Calling.sh` script with the following command line parameters:

| Flag | Description |
| ---------- | ---------- |
| `-d` | Path to directory containing sample directories |
| `-a` | Path to adapter file |
| `-c` | Indicate the minimum number of reads required to identify a peak |
| `-l` | Indicate the minimum number of libraries that a peak must be supported by |
| `-f` | Indicate whether or not peaks should be called by grouping sncRNAs into families (`yes` or `no`) |
| `-r` | Path to reference directory (e.g. `SCRAP`) |
| `-m` | Three-letter miRBase species abbreviation |
| `-g` | Reference genome abbreviation |

The adapter file can be the same as the adapter file used when running SCRAP or simply a `.txt` file with one sample name per row.

Example code for calling peaks with CLASH data previously analyzed with `SCRAP.sh`:

    bash Peak_Calling.sh \
        -d CLASH_Human/ \
        -a CLASH_Human/CLASH_Human_Adapters.txt \
        -c 3 \
        -l 2 \
        -f no \
        -r SCRAP/ \
        -m hsa \
        -g hg38

Peak calling will generate a `peaks.bed` (or `peaks.family.bed`) and `peakcalling.summary.txt` (or `peakcalling.family.summary.txt`) file in the directory denoted with the `-d` flag.

## Peak Annotation

The `peaks.bed` (or `peaks.family.bed`) file produced by `Peak_Calling.sh` can be annotated with gene names and features. Currently, this function is only available for annotating peaks called for data from human (*H. sapiens*), mouse (*M. musculus*), or worm (*C. elegans*).

Execute the `Peak_Annotation.sh` script with the following command line parameters:

| Flag | Description  |
| ---- | --------------- |
| `-p` | Path to peaks.bed or peaks.family.bed file     |
| `-r` | Path to reference directory (e.g. `SCRAP`)     |
| `-s` | Indicate species used for annotation (`human` (*H. sapiens*), `mouse` (*M. musculus*), or `worm` (*C. elegans*)  |

Example code for annotating peaks with CLASH data idnetified using `Peak_Calling.sh`:

    bash Peak_Annotation.sh \
        -p CLASH_Human/peaks.bed \
        -r SCRAP/ \
        -s human
