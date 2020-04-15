# Process ITS2 reads #
starttime <- Sys.time()

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(tidyverse)
library(phyloseq)
library(decontam)


# find files ####
path <- "../RawSeqsProcessing/DM-FAYE-LOG"

fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

fnFs <- fnFs[grep("ITS2",fnFs)]
fnRs <- fnRs[grep("ITS2",fnRs)]

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname){strsplit(basename(fname), "_")[[1]][1]}
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)


# inspect read qualities 
# plotQualityProfile(fnFs[1:4])
# plotQualityProfile(cutRs[1:4])


# filter and trim ####
filtFs <- file.path(path,"ITS2_filts",basename(fnFs))
filtRs <- file.path(path,"ITS2_filts",basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(3, 3), 
                     truncQ = 10, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)

# learn errors ####
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)


# de-replicate
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# sample inferrence ####
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# save intermediate output
saveRDS(dadaFs,"./output/dadaFs.RDS")
saveRDS(dadaRs,"./output/dadaRs.RDS")

# merge reads ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


# make sequence table ####
seqtab.mergers <- makeSequenceTable(mergers)
seqtab <- makeSequenceTable(dadaFs) # FWD seqs only
dim(seqtab)
dim(seqtab.mergers)

# remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.mergers.nochim <- removeBimeraDenovo(seqtab.mergers, method="consensus", multithread=TRUE, verbose=TRUE)

# inspect dist. of seq lengths ####
plot(table(nchar(getSequences(seqtab.nochim))))
plot(table(nchar(getSequences(seqtab.mergers.nochim))))

# track reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

saveRDS(seqtab.nochim,"./output/seqtab.nochim.RDS")
seqtab.nochim <- readRDS("./output/seqtab.nochim.RDS")

saveRDS(seqtab.mergers.nochim,"./output/seqtab.mergers.nochim.RDS")
seqtab.nochim <- readRDS("./output/seqtab.mergers.nochim.RDS")

# Assign taxonomy ####


# UNITE general FASTA release for eukaryotes 2 - Created: 2020-02-04 - DOI: 10.15156/BIO/786371
unite.ref <- "./taxonomy/UNITE_Euk_2020-02-04_non-dev.fasta.gz"  # CHANGE ME to location on your machine

# test
# testseqs <- getSequences(seqtab.nochim)[1:3]
# testtaxa <- assignTaxonomy(testseqs, unite.ref, multithread = TRUE, tryRC = FALSE, verbose = TRUE)

taxa.merged <- assignTaxonomy(seqtab.mergers.nochim, unite.ref, multithread = TRUE, tryRC = FALSE, verbose = TRUE)
saveRDS(taxa.merged,"./output/taxa.merged.RDS")


taxa.FWD <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = FALSE, verbose = TRUE)
saveRDS(taxa.FWD,"./output/taxa.FWD.RDS")




taxa.print.FWD <- taxa.FWD  # Removing sequence rownames for display only
rownames(taxa.print.FWD) <- NULL
head(taxa.print.FWD)


endtime <- Sys.time()

-difftime(starttime,endtime)

# Create phyloseq objects ####

# Load ALL metadata
meta <- read.csv("./metadata/full_metadata.csv",stringsAsFactors = FALSE)

meta <- meta[which(meta$Sample.Name. %in% row.names(seqtab)),]  # subset to samples present

# one for fwd-only table

# one for merged reads table

colnames(seqtab.mergers.nochim)
row.names(meta) <- meta$Sample.Name.

ps.merged <- phyloseq(sample_data(meta),
         otu_table(seqtab.mergers.nochim,taxa_are_rows = FALSE),
         tax_table(taxa.merged))


# ID and remove contaminant seqs ####
contams <- isContaminant(ps.merged,method = "prevalence",neg = "NegativeControl")
noncontams <- row.names(contams)[which(contams$contaminant == FALSE)]

ps.merged
ps.merged.noncontam <- subset_taxa(ps.merged,taxa_names(ps.merged) %in% noncontams)


# save them as RDS ####
saveRDS(ps.merged.noncontam,"./output/ITS2_phyloseq_object_no-contams.RDS")


# subset to LOG samples only
ps.juniper <- subset_samples(ps.merged.noncontam, Project == "Juniper" & NegativeControl == FALSE)


# remove non-fungi
ps.juniper <- subset_taxa(ps.juniper,Kingdom == "k__Fungi")

# save
saveRDS(ps.juniper,"./output/juniper_phyloseq_object_clean_ITS2.RDS")


