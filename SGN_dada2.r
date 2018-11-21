# Pedro Martinez Arbizu
# R script to process data using dada2
# dada2 commands based on https://benjjneb.github.io/dada2/tutorial.html

##########################################
# Edit this variables
##########################################
#working directory
setwd('/home/pmartinez/projects/metabarcoding/vsearch/mangan16meio/dadaR/')
#path to fastq files
path <- "/home/pmartinez/projects/metabarcoding/vsearch/mangan16meio/sequences/trim"
#file pattern
patternF <- "_R1_tr.fastq"
patternR <- "_R2_tr.fastq"
# length of reads
readL <- 300
# bases to truncate due to bad quality atend of read
truncN <- 40
#maximum expected errors allowed on forward
maxEF <- 2
#maximum expected errors allowed on reverse
maxER <- 3
###########################################
# Do not edit beyond this line
###########################################

#load packages
library(dada2)
library(vegan)

#check for access to files
list.files(path)

# here files names contain '_tr.' change to '_001.' if needed
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern=patternF, full.names = TRUE))
fnRs <- sort(list.files(path, pattern=patternR, full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#save the workspace
save.image('dada2.RData')

#if needed
# check quality of reads
pdf(file='quality_profileF.pdf')
plotQualityProfile(fnFs[1:2])
dev.off()

pdf(file='quality_profileR.pdf')
plotQualityProfile(fnRs[1:2])
dev.off()


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(readL-truncN,readL-truncN), maxN=0, maxEE=c(maxEF,maxER), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
#check result
head(out)

#learning errors on forward
errF <- learnErrors(filtFs, multithread=TRUE)
#learning errors on reverse
errR <- learnErrors(filtRs, multithread=TRUE)

#plot the trained error model
pdf(file='errorModelF.pdf')
plotErrors(errF, nominalQ=TRUE)
dev.off()
#plot the trained error model
pdf(file='errorModelR.pdf')
plotErrors(errR, nominalQ=TRUE)
dev.off()

# dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#inspect objects
dadaFs[[1]]
dadaRs[[1]]

#dadaFs[[1]]
#dada-class: object describing DADA2 denoising results
#1350 sequence variants were inferred from 46038 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

#dadaRs[[1]]
#dada-class: object describing DADA2 denoising results
#727 sequence variants were inferred from 139266 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
 
#merge paired reads 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#build seqeunce table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#[1]   19 3621

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#Identified 23391 bimeras out of 27012 input sequences.
#> dim(seqtab.nochim)
#[1]   19 3621


# track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

sum(seqtab.nochim)/sum(seqtab)
#[1] 0.8710768

#export the track data
write.csv(track,'trackdada2.cvs')

#export the dada2 sequences
write.table(colnames(seqtab.nochim),'nonchim.dada2.txt')

#export the community table
write.table(t(seqtab.nochim),'full.nonchim.dada2.txt',col.names=FALSE)

#save the workspace
save.image('dada2.RData')
###STOP here and assign taxonomy with dada2script on server


