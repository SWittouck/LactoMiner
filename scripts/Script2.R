#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#Load package
library(dada2); packageVersion("dada2")
#library(ggplot2)
library(tidyverse)

#-----------------------
#PARAMETERS;
loc1 =paste(args[1],"/FastqFiles", sep='')
loc2=paste(args[1],"/DB/BigData.fa.gz",sep='')
loc3=paste(args[2], "/taxa_SRA.csv",sep='')
loc4=paste(args[2], "/abundances_SRA.csv",sep='')
loc5=paste(args[2], "/abundances2_SRA.csv",sep='')
pathFastq <- toString(loc1)
pathDataB1 <- toString(loc2)
pathCSV <- toString(loc3)
patternF <- "_1.fastq" #Forward filenames
patternR <- "_2.fastq" #Reverse filenames
a <- c(150,150) #truncLen
b <- 0 #MaxN
c <- c(2,2) #MaxEE
d <- 2  #truncQ
e <- TRUE #rm.phix
#f <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) 
#g <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim)) # If processing a single sample

#-----------------------

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(pathFastq, patternF, full.names = TRUE))
fnRs <- sort(list.files(pathFastq, patternR, full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# a) basename removes all of the path up to and including the last path separator (if any).
# b) Split the elements of a character vector x into substrings according to the matches to substring split within them
#    met strsplit(x, split)
# c) sapply returns a list of the same length as X, each element of which is the result of applying FUN to the corresponding element of X.
#    met sapply(X, FUN). Meer info:https://stackoverflow.com/questions/19260951/using-square-bracket-as-a-function-for-lapply-in-r
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) x[1])

#Kwaliteit
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
#Place filtered files in filtered/ subdirectory
#Eerst FilterAndTrim, dan zippen (gz)
filtFs <- file.path(pathFastq, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(pathFastq, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=a,
                     maxN=b, maxEE=c, truncQ=d, rm.phix=e,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

#TIJDBESPARING
#filtFs=filtFs[1:2]
#filtRs=filtRs[1:2]

#Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#plot
plotErrors(errF, nominalQ=TRUE) + ggtitle("Errors Forward") + scale_y_continuous(breaks=c(0.001,0.100),trans='log10')
plotErrors(errR, nominalQ=TRUE) + ggtitle("Errors Reverse") + scale_y_continuous(breaks=c(0.001,0.100),trans='log10')

#dereplication
sample.namesF <- sapply(strsplit(basename(filtFs), "_"), function(x) x[1])
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), function(x) x[1])
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.namesF
names(derepRs) <- sample.namesR

#Sample inference aka Denoising
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
out1 <- out[out[, 2]!=0,]
track <- cbind(out1, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.namesF
head(track)

#Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, pathDataB1, multithread=TRUE)

#print taxonomy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Save taxa as .csv
abundances <- matrix(nrow = ncol(seqtab.nochim)*nrow(seqtab.nochim),ncol = 3)
sample<-rownames(seqtab.nochim)
sequence <- colnames(seqtab.nochim)
a=1
for (i in seq(1,nrow(seqtab.nochim))){
	  for (j in seq(1,ncol(seqtab.nochim))){
		      abundances[a,1] <- sample[i]
    abundances[a,2] <- seqtab.nochim[i,j]
        abundances[a,3] <- sequence[j]
        a=a+1
	  }
}

taxa <- cbind(rownames(taxa),taxa.print)
colnames(taxa)<- c("taxon_id","kingdom","phylum","class","order","family","genus")
colnames(abundances) <- c("sample_id","abundance","taxon_id")
abundances_df <- as.data.frame(abundances)
abundances2 <- spread(abundances_df,taxon_id,abundance)

write.csv2(abundances, loc4,row.names=FALSE)
write.csv2(taxa,pathCSV)
write.csv2(abundances2, loc5)
