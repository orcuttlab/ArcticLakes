#Adopted from the DADA2 Pipeline: https://benjjneb.github.io/dada2/tutorial.html
#These reads were sequenced at Integrated Microbiome Resource (IMR) using the 515FB / 926R primer pair (Parada et al., 2015, Walters et al., 2015)
#Raw read files are uploaded to the European Nucleotide Archive (ENA) under accession number PRJEB25188
#Processed data saved in RDS file format for upload and analysis in R are included in the repo
path <- "~/INSERT DIRECTORY" # Directory where files are located
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names
plotQualityProfile(fnFs[31:32]) # assess quality, decide where to trim
plotQualityProfile(fnRs[31:32]) # assess quality, decide where to trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200), 
              maxN=0, maxEE=c(2,2), trimLeft = 30, truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out) # Note trimming parameters 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(340,365)] # Selection of Amplicon length to keep
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/silva_nr_v138_train_set.fa.gz", multithread=TRUE) # Silva v138 Taxonomy Database
saveRDS(seqtab.nochim, file = "arctic_all_asv") # Files uploaded to Github
saveRDS(taxa, file = "arctic_all_taxa")# Files uploaded to Github

