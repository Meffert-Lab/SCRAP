#!/bin/sh

while getopts r:m:g:s: flag
do
    case "${flag}" in
        r) reference_directory=${OPTARG};;
        m) miRBase_species_abbreviation=${OPTARG};;
        g) genome_species_abbreviation=${OPTARG};;
	s) species=${OPTARG};;
    esac
done

if [ -z "$reference_directory" ]; then
    echo "Error: Path to reference directory not provided [-r]"
    exit 1
fi

if [ -z "$miRBase_species_abbreviation" ]; then
    echo "Error: miRBase species abbreviation not provided [-m]"
    exit 1
fi

if [ -z "$genome_species_abbreviation" ]; then
    echo "Error: Species genome abbreviation not provided [-g]"
    exit 1
fi

if [ -z "$species" ]; then
    echo "Error: Species not provided [-s]"
    exit 1
fi

	location=$(conda info | awk '/base environment/' | awk '{print $4}')
	source ${location}/etc/profile.d/conda.sh

	conda activate SCRAP

	cd ${reference_directory}
	cd fasta

	mkdir ${miRBase_species_abbreviation}
	mkdir ${genome_species_abbreviation}

#Subset miRNA sequences for species of interest
#miRBase.fasta obtained using:		wget https://www.mirbase.org/ftp/CURRENT/mature.fa.gz -O miRBase.fasta
#Final miRNA name format (example): miRNA-mmu-let-7c-5p

	awk '/^>/ {printf "%s%s ", pfx, $0; pfx="\n"; next} {printf "%s", $0} END {print ""}' miRBase.fasta | \
	grep "${miRBase_species_abbreviation}-" | \
	awk '{print $1"\n"$6}' | \
	sed 's/>/>miRNA-/g' \
	> ${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta

	grep -A1 "${miRBase_species_abbreviation}" tRFdb.fasta \
	> ${miRBase_species_abbreviation}/tRF_${miRBase_species_abbreviation}.fasta

	cat \
	${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta \
	${miRBase_species_abbreviation}/tRF_${miRBase_species_abbreviation}.fasta	\
	> ${miRBase_species_abbreviation}/sncRNA_${miRBase_species_abbreviation}.fasta

#Make BLAST database for sncRNA list

	makeblastdb \
	-in ${miRBase_species_abbreviation}/sncRNA_${miRBase_species_abbreviation}.fasta \
	-dbtype nucl

#Download genome

	wget http://hgdownload.cse.ucsc.edu/goldenpath/${genome_species_abbreviation}/bigZips/${genome_species_abbreviation}.fa.gz -O ${genome_species_abbreviation}/${genome_species_abbreviation}.fa.gz
	gunzip -c ${genome_species_abbreviation}/${genome_species_abbreviation}.fa.gz > ${genome_species_abbreviation}/${genome_species_abbreviation}.fa
	rm -r ${genome_species_abbreviation}/${genome_species_abbreviation}.fa.gz

#Index mouse genome for HISAT2 alignment
#Takes ~2 hours

	cd ${genome_species_abbreviation}

	hisat2-build ${genome_species_abbreviation}.fa ${genome_species_abbreviation}

	cd ..

#Download reference genome chromosome sizes

	wget http://hgdownload.cse.ucsc.edu/goldenpath/${genome_species_abbreviation}/bigZips/${genome_species_abbreviation}.chrom.sizes -O ${genome_species_abbreviation}/${genome_species_abbreviation}.chrom.sizes

#Download refFlat

	wget https://hgdownload.cse.ucsc.edu/goldenPath/${genome_species_abbreviation}/database/refFlat.txt.gz -O ${genome_species_abbreviation}/${genome_species_abbreviation}refFlat.txt

#Index mouse genome for GATK and Samtools

	gatk CreateSequenceDictionary -R ${genome_species_abbreviation}/${genome_species_abbreviation}.fa
	samtools faidx ${genome_species_abbreviation}/${genome_species_abbreviation}.fa

#Subset miRNA hairpin sequences for species of interest
#miRBase.hairpin.fasta obtained using:	wget https://www.mirbase.org/ftp/CURRENT/hairpin.fa.gz -O miRBase.hairpin.fasta

	awk '/^>/ {printf "%s%s ", pfx, $0; pfx="\n"; next} {printf "%s", $0} END {print ""}' miRBase.hairpin.fasta | \
	grep "${miRBase_species_abbreviation}-" | \
	awk '{print $1"\n"$7}' | \
	sed 's/>/>miRNA-hairpin-/g' \
	> ${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta

#Make BLAST database for miRNA hairpin sequences

	makeblastdb \
	-in ${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
	-dbtype nucl

#Subset tRNA sequences for species of interest

	grep -A1 ">${miRBase_species_abbreviation}_" GtRNAdb.fasta \
	> ${miRBase_species_abbreviation}/tRNA_${miRBase_species_abbreviation}.fasta

#Make BLAST database for miRNA hairpin sequences

	makeblastdb \
	-in ${miRBase_species_abbreviation}/tRNA_${miRBase_species_abbreviation}.fasta \
	-dbtype nucl

	cd ..
	cd annotation

	gunzip ${species}.annotation.bed.gz

	conda deactivate
