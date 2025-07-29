# Microbiome16S
Recommend making a environment just for Dada2.

# Load required packages
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(ggrepel)

BiocManager::install("decontam") 
library(decontam)  # Decontamination package

# Set path to your data
Make sure - data is trimmed 
```path <- "/Users/yourpath"  ## CHANGE to your fastq directory
list.files(path)```

# Identify forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_trimmed.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)  # Extracts "CP_S104" from "CP_S104_R1.fastq"

# Define control sample name **without** R1/R2
control_sample_name <- "CP_S104"

# Quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filter and trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(275,275),
                     maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy
*Implement with your own taxa assignments
taxa <- assignTaxonomy(seqtab.nochim, "/Users/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1")
rownames(taxa) <- getSequences(seqtab.nochim)
head(taxa)
