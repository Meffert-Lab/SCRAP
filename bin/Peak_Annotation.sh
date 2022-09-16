#!/bin/sh

while getopts p:r:s: flag
do
    case "${flag}" in
        p) peak_file=${OPTARG};;
        r) reference_directory=${OPTARG};;
        s) species=${OPTARG};;
    esac
done

if [ -z "$peak_file" ]; then
    echo "Error: Path to peak file not provided [-p]"
    exit 1
fi

if [ -z "$reference_directory" ]; then
    echo "Error: Path to reference directory not provided [-r]"
    exit 1
fi

if [ -z "$species" ]; then
    echo "Error: Species not provided [-s]"
    exit 1
fi

#Add / to the end of path to reference directory if user did not include one

reference_directory=$(echo "${reference_directory}" | sed 's![^/]$!&/!')

#For the most part, code modifications are not required below this line
##############################################################################

#Activate conda environemnt 

	location=$(conda info | awk '/base environment/' | awk '{print $4}')
	source ${location}/etc/profile.d/conda.sh

	conda activate SCRAP

	bedtools \
	intersect \
	-a ${peak_file} \
	-b ${reference_directory}annotation/${species}.annotation.bed \
	-s \
	-wa \
	-wb | \
	sed 's/\t/@/' | \
	sed 's/\t/@/' | \
	sed 's/\t/@/' | \
	awk '!x[$1]++' | \
	sed 's/@/\t/' | \
	sed 's/@/\t/' | \
	sed 's/@/\t/' \
	> ${peak_file}.tmp.annotation

	bedtools \
	intersect \
	-a ${peak_file} \
	-b ${reference_directory}annotation/${species}.annotation.bed \
	-s \
	-wa \
	-v | \
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,"--","--","--","--","intergenic","--"}' | \
	cat ${peak_file}.tmp.annotation - | \
	sort -k1,1 -k2,2n \
	> ${peak_file}.${species}.annotation.txt

	rm ${peak_file}.tmp.annotation

	conda deactivate
