#!/bin/sh

#Add links to tRNA FASTA files from GtRNAdb (http://gtrnadb.ucsc.edu)

hsa_link="http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa" #Homo Spapiens (hg38)
mmu_link="http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.fa" #Mus musculus (mm10)
rno_link="http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Rnorv7/rn7-tRNAs.fa" #Rattus norvegicus (rn7)
dre_link="http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Dreri11/danRer11-tRNAs.fa" #Danio rerio (danRer11)
xtr_link="http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Xtrop9/xenTro9-tRNAs.fa" #Xenopus tropicalis (xenTro9)
dme_link="http://gtrnadb.ucsc.edu/genomes/eukaryota/Dmela6/dm6-tRNAs.fa" #Drosophila melanogaster (dm6)
cel_link="http://gtrnadb.ucsc.edu/genomes/eukaryota/Celeg11/ce11-tRNAs.fa" #Caenorhabditis elegans (ce11)
ath_link="http://gtrnadb.ucsc.edu/genomes/eukaryota/Athal10/araTha1-tRNAs.fa" #Arabidopsis thaliana (TAIR10)

list_of_species="
hsa
mmu
rno
dre
xtr
cel
dme
ath"

#File contains line breaks within FASTA seqeunce that can cause downstream issues. Remove line breaks as described here: https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file. Add "hsa_" to the beginning of the read name so that it can be retrieved later. Print only the read name and sequence (remove SE ID, length, chromosome position, etc.).

for species in $list_of_species
do

url="${species}_link"

wget ${!url} -O ${species}.fa

awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${species}.fa | \
    sed 's/>/>'${species}'_/' | \
    awk '{print $1}' \
    > ${species}_tRNA.fasta

rm ${species}.fa

done

#Combine individual fasta files to create GtRNAdb reference FASTA file

touch tmp.combined.fasta

for species in $list_of_species
do

cat tmp.combined.fasta ${species}_tRNA.fasta > tmp2.combined.fasta

mv tmp2.combined.fasta tmp.combined.fasta

done

mv tmp.combined.fasta GtRNAdb.fasta
