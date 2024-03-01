#!/bin/sh

#Download all tRFs (http://genome.bioch.virginia.edu/trfdb/index.php)

#Convert to FASTA

	awk -F ', ' 'NR!=1 {print ">"$1"\n"$9}' tRF-1.csv > tRF-1.fasta
	awk -F ', ' 'NR!=1 {print ">"$1"\n"$9}' tRF-3.csv > tRF-3.fasta
	awk -F ', ' 'NR!=1 {print ">"$1"\n"$9}' tRF-5.csv > tRF-5.fasta

#Dedupe 

	awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tRF-1.fasta | \
	awk '!x[$1]++' | \
	sed 's/ /\n/' | \
	sed 's/>/>tRF-/g' \
	> tRF-1.deduped.fasta

	awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tRF-3.fasta | \
	awk '!x[$1]++' | \
	sed 's/ /\n/' | \
	sed 's/>/>tRF-/g' \
	> tRF-3.deduped.fasta

	awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tRF-5.fasta | \
	awk '!x[$1]++' | \
	sed 's/ /\n/' | \
	sed 's/>/>tRF-/g' \
	> tRF-5.deduped.fasta

	cat \
	tRF-1.deduped.fasta \
	tRF-3.deduped.fasta \
	tRF-5.deduped.fasta \
	> tRF.fasta
