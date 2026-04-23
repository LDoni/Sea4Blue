# get ready bastards https://benjjneb.github.io/dada2/tutorial.html
#https://compbiocore.github.io/metagenomics-workshop/assets/DADA2_tutorial.html

#!!!!!!!!!!!!!!!!!!!!!!!    ---> leggere questo https://blogs.oregonstate.edu/earthmotes/2021/09/28/dada2-pipeline-for-16s-datasets-in-r/



library(DECIPHER)
library(phangorn)
library(dada2)
library(ShortRead)
library(Biostrings)
library(ShortRead)
path <- "RAWREADS/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


#
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#plot quality
plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnRs[1:10])






#check for primers
#515F FWD: GTGYCAGCMGCCGCGGTAA   #19 bp
#806R REV: GGACTACNVGGGTWTCTAAT   #20bp
#amplicon size: 806-515=291

FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACNVGGGTWTCTAAT"  ##





allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#prefiltering for remove Ns

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)


primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#                Forward Complement Reverse RevComp

#!!! non ci sono i primers

#da terminale:
#python3 ~/figaro/figaro/figaro.py -i . -o figaro -a 252 -f 1 -r 1

#consiglia: {"trimPosition": [148, 163], "maxExpectedError": [1, 2], "readRetentionPercent": 94.57, "score": 93.56757419671327}


#python3 ~/figaro/figaro/figaro.py -i . -o figaro1 -a 252 -f 1 -r 1
#{"trimPosition": [124, 150], "maxExpectedError": [1, 2], "readRetentionPercent": 96.02, "score": 95.02452783909737}

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered_figaro", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered_figaro", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names



out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(124,150),
                     maxN=0, maxEE=c(1,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
out<-data.frame(out)
out$perc<-(out$reads.out/out$reads.in*100)
head(out)



#learn errors
errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)
#visualize
plotErrors(errF, nominalQ=TRUE)

#sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#merge pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


#construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
#percentage of chimeras
dim(seqtab.nochim)[2]/dim(seqtab)[2]

#chimera abundances 
1-sum(seqtab.nochim)/sum(seqtab)



#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))


# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track<-data.frame(track)
options(digits=3)
track$percInpOut= (track$nonchim/ track$input*100 )
head(track)
write.csv(track, "dada2_denoising_stats.csv")


#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

#inspect
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

grep("Vibrio",taxa.print)

#Evaluating DADA2’s accuracy on the mock community:
rownames(seqtab.nochim)
unqs.mock <- seqtab.nochim["Lapospositive",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")




##comparare le seq con i fasta mock
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")




save.image("DADA2.rdata")
# Close R, Re-open R
#load("DADA2.rdata")




#Construct Phylogenetic Tree
sequences<-getSequences(seqtab.nochim)
names(sequences)<-sequences
#Run Sequence Alignment (MSA) using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
#Change sequence alignment output into a phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
#Create distance matrix
dm <- dist.ml(phang.align)
#Perform Neighbor joining
treeNJ <- NJ(dm) # Note, tip order != sequence order
#Internal maximum likelihood
fit = pml(treeNJ, data=phang.align)
#negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", 
                    control = pml.control(trace = 0))

head(fitGTR$tree)


#https://github.com/dermilke/Pacific_Bacterioplankton

library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(read.csv("metadata.csv",row.names=1, check.names=FALSE,sep = ",")), 
               tax_table(taxa),phy_tree(fitGTR$tree)) #

#save refseq nella slot seqs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

sample_data(ps)
#save phyloseq object
saveRDS(ps, "ps.rds")
#ps <- readRDS("ps.rds")




phyloseq_obj_css


# Extract abundance matrix from the phyloseq object




OTU<-(as.data.frame(otu_table(phyloseq_obj_css)))
OTU$otu<-row.names(OTU)


TAX<-(as.data.frame(tax_table(phyloseq_obj_css)))
TAX$otu<-row.names(TAX)
head(TAX)

write_csv(merge(OTU,TAX, by="otu"),"otu_tax.csv")
