#!/bin/sh

while getopts d:a:r:c:l:f:m:g: flag
do
    case "${flag}" in
        d) directory=${OPTARG};;
        a) adapter_file=${OPTARG};;
        c) minimum_reads=${OPTARG};;
        l) minimum_libraries=${OPTARG};;
        f) family=${OPTARG};;
        r) reference_directory=${OPTARG};;
        m) miRBase_species_abbreviation=${OPTARG};;
        g) genome_species_abbreviation=${OPTARG};;
    esac
done

if [ -z "$directory" ]; then
    echo "Error: Path to sample directories not provided [-d]"
    exit 1
fi

if [ -z "$adapter_file" ]; then
    echo "Error: Path to adapter file not provided [-a]"
    exit 1
fi

if [ -z "$minimum_reads" ]; then
    echo "Error: Minimum read count not set [-c]"
    exit 1
fi

if [ -z "$minimum_libraries" ]; then
    echo "Error: Minimum library count not set [-l]"
    exit 1
fi

if [ -z "$family" ]; then
    echo "Error: Call peaks by miRNA family? [-f]"
    exit 1
fi

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

#Add / to the end of path to directory if user did not include one

directory=$(echo "${directory}" | sed 's![^/]$!&/!')

#Add / to the end of path to reference directory if user did not include one

reference_directory=$(echo "${reference_directory}" | sed 's![^/]$!&/!')

#For the most part, code modifications are not required below this line
##############################################################################

#Activate conda environemnt 

	location=$(conda info | awk '/base environment/' | awk '{print $4}')
	source ${location}/etc/profile.d/conda.sh

	conda activate SCRAP

samples=$(awk '!/^#/' ${adapter_file} | awk '{print $1}')

#Delete and temporary files

	[ -e ${directory}PeakCalling/combined.rearranged.bed ] && rm ${directory}PeakCalling/combined.rearranged.bed
	[ -e ${directory}PeakCalling/potential.peaks.bed ] && rm ${directory}PeakCalling/potential.peaks.bed
	[ -e ${directory}PeakCalling/cycles.txt ] && rm ${directory}PeakCalling/cycles.txt
	[ -e ${directory}PeakCalling/filtered.peaks.bed ] && rm ${directory}PeakCalling/filtered.peaks.bed
	[ -e ${directory}PeakCalling/peakcalling.summary.txt ] && rm ${directory}PeakCalling/peakcalling.summary.txt

#Make a directory in which the intermediate files generated during peak calling will be stored

	mkdir ${directory}PeakCalling

	echo "Start: $(date)" >> ${directory}PeakCalling/peakcalling.summary.txt

for sample in ${samples}
do

#Convert bam file to sam file

	samtools view \
	-h \
	${directory}${sample}/${sample}.aligned.unique.bam \
	> ${directory}PeakCalling/${sample}.aligned.unique.sam

if [ $family == "yes" ]
then

	awk '!/miRNA-/' ${directory}${sample}/${sample}.aligned.unique.sam \
	> ${directory}${sample}/${sample}.sam.nonmiRNA.txt

	join \
	-a 1 \
	-1 2 -2 1 \
	<(awk '/miRNA-/' ${directory}${sample}/${sample}.aligned.unique.sam | sed 's/\.miRNA-/\.miRNA-\t/' | sort -k2,2) \
	<(grep "${miRBase_species}-" ${reference_directory}annotation/miR_Family.txt | sort -k1,1) | \
	awk '{$1=""}1' | \
	awk '{$1=$1}1' | \
	awk '{$1 = $1$NF; print}' | \
	awk 'NF{NF-=1};1' | \
	sed 's/ /\t/g' \
	> ${directory}${sample}/${sample}.aligned.unique.tmp.sam

	cat ${directory}${sample}/${sample}.sam.nonmiRNA.txt ${directory}${sample}/${sample}.aligned.unique.tmp.sam \
	> ${directory}PeakCalling/${sample}.aligned.unique.tmp.sam

	rm ${directory}${sample}/${sample}.sam.nonmiRNA.txt

else

#If not calling peaks by family, move and rename sam file to peak calling folder

	cp ${directory}${sample}/${sample}.aligned.unique.sam ${directory}PeakCalling/${sample}.aligned.unique.tmp.sam

fi

#Convert sam file to bam file

	samtools view \
	-S \
	-h \
	-b ${directory}PeakCalling/${sample}.aligned.unique.tmp.sam \
	> ${directory}PeakCalling/${sample}.aligned.unique.bam

#Split reads across unmapped regions (e.g. introns)

	gatk SplitNCigarReads \
	-R ${reference_directory}fasta/${genome_species_abbreviation}/${genome_species_abbreviation}.fa \
	-I ${directory}PeakCalling/${sample}.aligned.unique.bam \
	-O ${directory}PeakCalling/${sample}.aligned.unique.split.bam

#Convert bam file to bed file

	bedtools \
	bamtobed \
	-i ${directory}PeakCalling/${sample}.aligned.unique.split.bam \
	> ${directory}PeakCalling/${sample}.aligned.unique.split.bed

#Append sncRNA (or family) name to chromosome number to allow for all peaks to be called at once

	sed 's/\./\t/' ${directory}PeakCalling/${sample}.aligned.unique.split.bed | \
	awk '{print $1"."$5"\t"$2"\t"$3"\t""X""\t""X""\t"$7}' \
	> ${directory}PeakCalling/${sample}.aligned.unique.split.rearranged.bed

#Combine bed files for each sample into a single file

	cat "${directory}PeakCalling/${sample}.aligned.unique.split.rearranged.bed" >> ${directory}PeakCalling/combined.rearranged.bed

done

#Sort combined bed file for easier processing

	sort -k1,1 -k2,2n ${directory}PeakCalling/combined.rearranged.bed > ${directory}PeakCalling/combined.rearranged.sorted.bed

#Create a temporary copy of the combined bed file that will be modified during successive rounds of peak calling

	cp ${directory}PeakCalling/combined.rearranged.sorted.bed ${directory}PeakCalling/combined.rearranged.sorted.tmp.bed

#Perform cycles of peak calling

x=1
while [[ $x > 0 ]]
do

#Merge overlapping reads

	bedtools \
	merge \
	-i ${directory}PeakCalling/combined.rearranged.sorted.tmp.bed \
	-s \
	-c 6 \
	-o distinct \
	> ${directory}PeakCalling/combined.rearranged.sorted.merged

#Determine coverage across overlapping regions
#Remove regions with fewer than minimum read threshold (assigned at beginning of script)
#Remaining regions termed "windows"

	bedtools \
	coverage \
	-a <(awk '{print $1"\t"$2"\t"$3"\t""X""\t""X""\t"$4}' ${directory}PeakCalling/combined.rearranged.sorted.merged) \
	-b ${directory}PeakCalling/combined.rearranged.sorted.tmp.bed \
	-s \
	-d | \
	awk -v var=${minimum_reads} '$8 >= var' \
	> ${directory}PeakCalling/combined.rearranged.sorted.merged.coverage

#Count the number of windows with minimal read coverage to see if another round of peak calling should be completed

	x=$(wc -l ${directory}PeakCalling/combined.rearranged.sorted.merged.coverage | awk '{print $1}')

#Find start and stop position of maximal coverage within window

	paste \
	<(sort -k1,1 -k2,2n -k8,8nr -k7,7n ${directory}PeakCalling/combined.rearranged.sorted.merged.coverage | awk '!window[$1, $2, $3, $6]++' | cut -f 1-3,6,7) \
	<(sort -k1,1 -k2,2n -k8,8nr -k7,7nr ${directory}PeakCalling/combined.rearranged.sorted.merged.coverage | awk '!window[$1, $2, $3, $6]++' | cut -f 7) | \
	awk '{print $1"\t"($2+$5-1)"\t"($2+$6)"\t""X""\t""X""\t"$4}' \
	> ${directory}PeakCalling/combined.rearranged.sorted.merged.coverage.filtered

	cat ${directory}PeakCalling/combined.rearranged.sorted.merged.coverage.filtered >> ${directory}PeakCalling/potential.peaks.bed

#Remove reads that overlap with peak

	bedtools \
	subtract \
	-a ${directory}PeakCalling/combined.rearranged.sorted.tmp.bed \
	-b ${directory}PeakCalling/combined.rearranged.sorted.merged.coverage.filtered \
	-s \
	-A \
	> ${directory}PeakCalling/combined.rearranged.sorted.merged.coverage.filtered.outside

#Replace temporary bed file with bed file that reads were removed from

	mv ${directory}PeakCalling/combined.rearranged.sorted.merged.coverage.filtered.outside ${directory}PeakCalling/combined.rearranged.sorted.tmp.bed

	echo "Cycle complete" >> ${directory}PeakCalling/cycles.txt

done

#Sort bed file for easier processing

	sort -k1,1 -k2,2n ${directory}PeakCalling/potential.peaks.bed > ${directory}PeakCalling/potential.peaks.sorted.bed

#Determine coverage for each individual sample

for sample in ${samples}
do

	bedtools \
	coverage \
	-a ${directory}PeakCalling/potential.peaks.sorted.bed \
	-b <(sort -k1,1 -k2,2n ${directory}PeakCalling/${sample}.aligned.unique.split.rearranged.bed) \
	-s \
	> ${directory}PeakCalling/${sample}.aligned.unique.split.rearranged.coverage.bed

done

#Combine coverages for each sample into a single file
#Each column represents a sample

	first=$(echo "${samples}" | sed -r '/^\s*$/d' | awk 'NR==1')

	[ -e ${directory}PeakCalling/tmp.tmp ] && rm ${directory}PeakCalling/tmp.tmp

	cut -f 1-6 ${directory}PeakCalling/${first}.aligned.unique.split.rearranged.coverage.bed > ${directory}PeakCalling/tmp.tmp

for sample in ${samples}
do

	cut \
	-f 7 \
	${directory}PeakCalling/${sample}.aligned.unique.split.rearranged.coverage.bed | \
	paste ${directory}PeakCalling/tmp.tmp - \
	> ${directory}PeakCalling/tmp2.tmp

	mv ${directory}PeakCalling/tmp2.tmp ${directory}PeakCalling/tmp.tmp

done

	mv ${directory}PeakCalling/tmp.tmp ${directory}PeakCalling/potential.peaks.sorted.coverage.bed

#Remove peaks supported by fewer than minimum number of libraries (assigned at beginning of script)

	awk \
	-v var=${minimum_libraries} \
	-v var2=$((6 + $(echo "$samples" | sed -r '/^\s*$/d' | wc -l) )) \
	'{nz=0; for(i=7;i<=var2;i++) nz+=($i!=0)} nz>=var' \
	${directory}PeakCalling/potential.peaks.sorted.coverage.bed \
	>> ${directory}PeakCalling/filtered.peaks.bed

	bedtools \
	coverage \
	-a <(awk '{print $1"\t"$2"\t"$3"\t""X""\t""X""\t"$6}' ${directory}PeakCalling/filtered.peaks.bed | sort -k1,1 -k2,2n) \
	-b <(awk '{print $1"\t"$2"\t"$3"\t""X""\t""X""\t"$6}' ${directory}PeakCalling/combined.rearranged.bed | sort -k1,1 -k2,2n) \
	-s | \
	sed 's/\./\t/' | \
	awk '{print $1"\t"$3"\t"$4"\t"$2"\t"$8"\t"$7}' \
	> ${directory}PeakCalling/peaks.bed

if [ $family == "yes" ]
then

	cp ${directory}PeakCalling/peaks.bed ${directory}peaks.family.bed

else

	cp ${directory}PeakCalling/peaks.bed ${directory}peaks.bed

fi

	echo "Finish: $(date)" >> ${directory}PeakCalling/peakcalling.summary.txt
	echo "Samples:" >> ${directory}PeakCalling/peakcalling.summary.txt
	echo "${samples}" >> ${directory}PeakCalling/peakcalling.summary.txt
	echo "miRBase species abbreviation: ${miRBase_species_abbreviation}" >> ${directory}PeakCalling/peakcalling.summary.txt
	echo "Genome species abbreviation: ${genome_species_abbreviation}" >> ${directory}PeakCalling/peakcalling.summary.txt
	echo "Peaks called by family? ${family}" >> ${directory}PeakCalling/peakcalling.summary.txt
	echo "Minimum reads: ${minimum_reads}" >> ${directory}PeakCalling/peakcalling.summary.txt
	echo "Minimum libraries: ${minimum_libraries}" >> ${directory}PeakCalling/peakcalling.summary.txt
	echo "Cycles complete: $(wc -l ${directory}PeakCalling/cycles.txt | awk '{print $1}')" >> ${directory}PeakCalling/peakcalling.summary.txt
	echo "Total peaks: $(wc -l ${directory}PeakCalling/peaks.bed | awk '{print $1}')" >> ${directory}PeakCalling/peakcalling.summary.txt

if [ $family == "yes" ]
then

	cp ${directory}PeakCalling/peakcalling.summary.txt ${directory}peakcalling.family.summary.txt

else

	cp ${directory}PeakCalling/peakcalling.summary.txt ${directory}peakcalling.summary.txt

fi

	rm -r ${directory}PeakCalling

conda deactivate
