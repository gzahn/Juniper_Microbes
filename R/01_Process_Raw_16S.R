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

fnFs <- fnFs[grep("16S",fnFs)]
fnRs <- fnRs[grep("16S",fnRs)]

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname){strsplit(basename(fname), "_")[[1]][1]}
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)


# inspect read qualities 
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])


# filter and trim ####
filtFs <- file.path(path,"16S_filts",basename(fnFs))
filtRs <- file.path(path,"16S_filts",basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)

# reassign filts for lost samples, if any
filtFs <- list.files(file.path(path,"16S_filts"), pattern = "R1_001.fastq.gz", full.names = TRUE)
filtRs <- list.files(file.path(path,"16S_filts"), pattern = "R2_001.fastq.gz", full.names = TRUE)


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
saveRDS(dadaFs,"./output/dadaFs_16S.RDS")
saveRDS(dadaRs,"./output/dadaRs_16S.RDS")

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
write.csv(track,"./output/tracked_reads_16S.csv",row.names = FALSE,quote = FALSE)

saveRDS(seqtab.nochim,"./output/seqtab.nochim.RDS")
seqtab.nochim <- readRDS("./output/seqtab.nochim.RDS")

saveRDS(seqtab.mergers.nochim,"./output/seqtab.mergers.nochim.RDS")
seqtab.nochim <- readRDS("./output/seqtab.mergers.nochim.RDS")

# Assign taxonomy ####

taxa <- assignTaxonomy(seqtab.mergers.nochim, "./taxonomy/rdp_train_set_16.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "./taxonomy/rdp_species_assignment_16.fa.gz")

saveRDS(taxa,"./output/taxa.merged.RDS")

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


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
                      tax_table(taxa))


# ID and remove contaminant seqs ####
contams <- isContaminant(ps.merged,method = "prevalence",neg = "NegativeControl")
noncontams <- row.names(contams)[which(contams$contaminant == FALSE)]

ps.merged
ps.merged.noncontam <- subset_taxa(ps.merged,taxa_names(ps.merged) %in% noncontams)


# save them as RDS ####
saveRDS(ps.merged.noncontam,"./output/16S_phyloseq_object_no-contams.RDS")


# subset to LOG samples only
ps.juniper <- subset_samples(ps.merged.noncontam, Project == "Juniper" & NegativeControl == FALSE)


# remove non-fungi
ps.juniper <- subset_taxa(ps.juniper,Kingdom == "Bacteria")

# save
saveRDS(ps.juniper,"./output/juniper_phyloseq_object_clean_16S.RDS")


