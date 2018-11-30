#!/bin/bash  
### Trimm primers and adapters with bbmap
###use bbmap from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/

# run like: bash ./trimmadpat.sh  > log.txt 2>&1 &

# force trimm V1V2 primers
#SSU-F_04 GCTTGTCTCAAAGATTAAGCC
#SSU-R_22 GCCTGCTGCCTTCCTTGGA
#this seqeunces were added to $bbpath/resources/adapters.fa 

# "for f in *_R1_*.fastq; do ..."
# check if xou have gzipped files
# then "for f in *_R1_*.fastq.gz; do ..."


#####################################
### EDIT THESE OPTIONS            ###
#####################################

	# bbmap path
	bbpath='/home/pmartinez/Downloads/bbmap'

	#Primer V1V2
	primer='GCTTGTCTCAAAGATTAAGCC,GCCTGCTGCCTTCCTTGGA'

	# kmer size for searching primer
	# should be at most size of primer
	k=15
	
	# Hamming distance allowed
	hd=1

#####################################
### Do NOT EDIT BEXOND THIS LINE  ###
#####################################



	for f in *_R1_*.fastq.gz; do
	    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
 	    s=$(cut -d_ -f1 <<< "$f")
	   	   

	echo ====================================	
	echo "Step 1 : Trimm adapter overhang using paired-end reads and primer sequence"
	echo ====================================
	echo
	#$bbpath/bbduk.sh in1=$f in2=$r out1=$s\_R1_tr.fastq out2=$s\_R2_tr.fastq ref=$bbpath/resources/adapters.fa \
	#ktrim=l k=$k hdist=$hd rcomp=t tbo overwrite=true

	$bbpath/bbduk.sh in1=$f in2=$r out1=$s\_R1_tr.fastq out2=$s\_R2_tr.fastq literal=$primer \
	ktrim=l k=$k hdist=$hd rcomp=t tbo overwrite=true
	
	echo 
	echo

	done


