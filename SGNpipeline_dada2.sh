#!/bin/bash                                                                       

# SGN pipeline for mapping dada2 results to blast hits
# create 5 subfolders
# /scripts (holding the scripts, execute from this directory)
# /input (save the dada2 files here)
# /output for hoding the output files (note that they will be overwritten if present)
# /output/blast inside the output folder for holding the output of blast and the final Taxon table
# /blastdb for holding your blast databases and own library
#
#	for cluster with local installation of vsearch uncomment and edit following line
#	export vsearch='/home/yourloginname/vsearch/bin/vsearch'
#	edit variable VSEARCH=$vsearch
#
#	edit variable THREADS = number of threads or cores for parallel computing
     

# run in ubuntu with bash ./SGNpipeplinemarch_dada2.sh

# to log the output and errors to a file do:
#  bash ./SGNpipeline.sh  >> ../output/logfile.txt 2>&1 &


#####################################
### EDIT THESE OPTIONS            ###
#####################################
# running options
	THREADS=8
	PERL=$(which perl)
	VSEARCH=$(which vsearch)

# DADA2 files
	DADA2OUT=nonchim.dada2.txt # file with non chimeric sequences
	DADA2Tab=full.nonchim.dada2.txt # Community table from dada2 without col and row names

# Amplicon option
	GENE=V1V2 # Fragment name

# options for blast
	NBALIGN=10 # how many hits to keep from genebank blast
	OWNLIB=Ref_Library_Metabarcoding_18s.fasta # the name of your own library concatenated with genebank blast hits
	BLASTP='/usr/bin'
	BLASTDB='/home/pmartinez/metabarcoding/vsearch/blastdb'
#####################################
### DO NOT EDIT BEYOND THIS LINE  ###
#####################################


	
# cleaning up
	echo
	echo ====================================
	echo Removing old files from output directory
	echo ====================================

	date
	rm -r ../output
	mkdir ../output
	mkdir ../output/blast
	mkdir ../output/taxa
	
	echo
	echo ====================================
	echo Processing gene $GENE from DADA2
	echo ====================================

# Enter subdirectory with files                                         

	echo
	echo Checking FASTQ format version 
	cd ../input

# edit file to create fasta format
	sed -i 's/"x"//g' $DADA2OUT #exclude initial "x"
	sed -i '/^\s*$/d' $DADA2OUT #exclude empty line
	sed -i 's/"//g' $DADA2OUT # exclude all "
	sed -i 's/^/>/g' $DADA2OUT # create a > at begin of line
	sed -E -i  's/[0-9]{1,10}/&\n/' $DADA2OUT # all new line after last numeral
	sed -i 's/^ //g' $DADA2OUT # remove space at begin of line

# edit community table file to remove colum names
	sed -E -i 's/[AGTC"]{1,1000}/&\n/' $DADA2Tab
	sed  -i '/^"/d' $DADA2Tab


# Blast against Genebank, combine with own library and blast again
# Change the variant name with Taxon name
 
	echo ###################################	
	echo query agains blastdb
	echo ###################################	
	
##### Adapt here to sahar $BLASTDB
#SAHAR	
	cd ../blastdb
	blastn -db nt -query ../input/$DADA2OUT \
	-out ../output/blast/all.otus.$GENE.dada2.blast.txt \
	-num_alignments $NBALIGN  -num_threads $THREADS \
	-outfmt "6  qseqid pident length  sblastnames sscinames sacc evalue sseq"


	
	echo ###################################	
	echo deduplicating fasta from blast results
	echo ###################################	

	awk -F'\t' '{print ">"$6"|"$4"|"$5}'  \
	../output/blast/all.otus.$GENE.dada2.blast.txt \
	| sort -uk1,1  > ../output/blast/all.otus.$GENE.dada2.blast.dedup.txt

	# exclude hit with low taxonomic resolution
	grep -v  'eukaryotes|uncultured\|eukaryotes|eukaryote\|animals|uncultured\|fungi|uncultured' \
	../output/blast/all.otus.$GENE.dada2.blast.dedup.txt > \
	../output/blast/all.otus.$GENE.dada2.blast.dedup.clean.txt



	# create list of accession numbers
	sed 's/ /_/g' ../output/blast/all.otus.$GENE.dada2.blast.dedup.clean.txt \
	| awk -F'[>\|]' '{print $2}' > ../output/blast/all.otus.$GENE.dada2.blast.sacc.txt	 

	
##!! SAHAR
    # retrieve full sequences with accession numbers
#	$BLASTP/blastdbcmd -db $BLASTDB/nt \
#	-entry_batch ../output/blast/all.otus.$GENE.dada2.blast.sacc.txt \
#	-outfmt "%s" > ../output/blast/all.otus.$GENE.dada2.blast.fullseq.txt	
#PEDRO
    # retrieve full sequences with accession numbers
	#-target_only # with this option produces no result
	#blastdbcmd -db nt -entry MG253167 -outfmt "%a|%s" -target_only 
	blastdbcmd -db nt \
	-entry_batch ../output/blast/all.otus.$GENE.dada2.blast.sacc.txt \
	-outfmt "%a|%s"  > ../output/blast/all.otus.$GENE.dada2.blast.fullseq.txt	

	# remove version number
	sed -i -e 's/\..\?//g' ../output/blast/all.otus.$GENE.dada2.blast.fullseq.txt


    	echo ###################################	
	echo deduplicate and create fasta_file from blast results
	echo ###################################	
	
	awk -F'[>\|]' 'NR==FNR{a[$1]=$2;next} ($2 in a) {print ">"$2"|"$3"|"$4"\n"a[$2]}'  \
	../output/blast/all.otus.$GENE.dada2.blast.fullseq.txt \
	../output/blast/all.otus.$GENE.dada2.blast.dedup.clean.txt \
	> ../output/blast/all.otus.$GENE.dada2.blast.dedup.fasta	
	#change space to _	
	sed -i 's/ /_/g' ../output/blast/all.otus.$GENE.dada2.blast.dedup.fasta  	
	
	echo ###################################
	echo concatenate blast hits with own library
	echo ###################################

###SAHAR
#	cat ../output/blast/all.otus.$GENE.dada2.blast.dedup.fasta $BLASTDB/$OWNLIB \
#	    > $BLASTDB/finaldb.fasta

	cat ../output/blast/all.otus.$GENE.dada2.blast.dedup.fasta $OWNLIB \
	    > finaldb.fasta

	echo ###################################
	echo making blast database ...
	echo ###################################

###SAHAR
#	$BLASTP/makeblastdb -in $BLASTDB/finaldb.fasta -parse_seqids -dbtype nucl
	
	makeblastdb -in finaldb.fasta -parse_seqids -dbtype nucl

	echo ###################################
	echo quey against own blast db all.otus.$GENE.dada2.fasta ...
	echo ###################################
##SAHAR
#	$BLASTP/blastn -db $BLASTDB/finaldb.fasta -query \
#	    ../input/$DADA2OUT \
#	    -out ../output/blast/all.otus.$GENE.dada2.final.blast.txt -num_alignments 1  \
#	    -num_threads $THREADS \
#	    -outfmt "10  qseqid pident qcovs evalue length sacc "

	
	blastn -db finaldb.fasta -query \
	    ../input/$DADA2OUT \
	    -out ../output/blast/all.otus.$GENE.dada2.final.blast.txt -num_alignments 1  \
	    -num_threads $THREADS \
	    -outfmt "10  qseqid pident qcovs evalue length sacc"

	#replace "|" ,  "," and ";" by space
	sed -i 's/|/ /g;s/,/ /g;s/;/ /g' ../output/blast/all.otus.$GENE.dada2.final.blast.txt 

### added 12.5.2018
	echo ###################################	
	echo deduplicating blast results extracting only first species name provided by blast
	echo ###################################	
	
	awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8}' ../output/blast/all.otus.$GENE.dada2.final.blast.txt | \
	sort -uk1,1 -V > ../output/blast/all.otus.$GENE.dada2.final.dedup.blast.txt

	#format as table
	column -t ../output/blast/all.otus.$GENE.dada2.final.dedup.blast.txt \
	    > ../output/blast/all.otus.$GENE.dada2.final.blast.tab.txt

	column -t ../input/$DADA2Tab \
	    > ../output/blast/all.otutab.$GENE.dada2.tab.txt

####################################################################################TRy since here
	#rename the otu table

	paste ../output/blast/all.otus.$GENE.dada2.final.blast.tab.txt \
	 ../output/blast/all.otutab.$GENE.dada2.tab.txt \
	 > ../output/blast/Taxon_table.$GENE.dada2.txt
  
#	sed -i 's/^ ID/accn Group OTU Species pident qcovs length/g' ../output/blast/Taxon_table.$GENE.dada2.txt

	#create fasta file with original OTU sequence and description from blast
	grep '^[AGTC]' ../input/$DADA2OUT > ../output/blast/otuseqs.$GENE.dada2.txt

	paste -d '\t' ../output/blast/all.otus.$GENE.dada2.final.dedup.blast.txt ../output/blast/otuseqs.$GENE.dada2.txt \
	| awk  '{print ">"$6"|"$7"|"$8"|"$1"|"$2"|"$3"|"$4"|"$5"\n"$9}' \
	>  ../output/taxa/All.OTUS.blast.$GENE.dada2.fasta



	#extract by taxon
	declare -a taxa=($(awk -F '|' '{print $2}' ../output/taxa/All.OTUS.blast.$GENE.dada2.fasta | sort -uk1,1))
	
	for taxon in "${taxa[@]}"; do
	grep -A1 --no-group-separator $taxon ../output/taxa/All.OTUS.blast.$GENE.dada2.fasta \
	 > ../output/taxa/$taxon.$GENE.dada2.fasta

	done #done extract by taxon

	# clean up temporary files
	#rm ../output/blast/otuseqs.$GENE.dada2.txt
	#rm ../output/blast/all.otus.$GENE.dada2.final.blast.tab.txt 
	#rm ../output/blast/all.otus.$GENE.dada2.final.dedup.blast.txt
	#rm ../output/blast/all.otus.$GENE.dada2.blast.dedup.fasta
	#rm ../output/blast/all.otutab.$GENE.dada2.tab.txt 
	#rm ../output/blast/all.otus.$GENE.dada2.blast.sacc.txt 
	#rm ../output/blast/all.otus.$GENE.dada2.blast.fullseq.txt	


#done # done OTUSIM loop


# Game over
	echo
	echo Done
	date


