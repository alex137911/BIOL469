# Preliminary analysis of genes impacting cancers

# Remove objects in workspace
rm(list = ls())

# Required packages
suppressMessages(library(rtracklayer))

# -------------------------------------------------------------------
# Import MANE data
# Downloaded May 16, 2023 (Release 1.1, hg38 genome assembly)
# https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/

inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/4A/BIOL 469/Final Project/Data/Input"
inDir  <- sprintf("%s/Input", dirname(inpath))
setwd(inDir)

MANE   <- readGFF("MANE.GRCh38.v1.1.ensembl_genomic.gff.gz")
MANEdf <- as.data.frame(MANE)
message(sprintf("Loaded %i transcripts", nrow(MANEdf)))                                 # 521 393 transcripts

# -------------------------------------------------------------------
# PRE-PROCESSING

# Filter MANE data for start codons
startCodons <- subset(MANEdf, MANEdf$type == "start_codon")
message(sprintf("%i number of start codons", nrow(startCodons)))                        # 19 825 transcripts

# Drop values after the decimal
startCodons$transcript_id <- sub("\\..*", "", startCodons$transcript_id)

# Remove sequences with "MANE_Plus_Clinical" tag
# Want only canonical transcripts (so that there is one transcript per gene)
# "MANE_Plus_Clinical" tag includes transcripts with 
# clinical significance (e.g. disease-causing mutations)
startCodons_filtered <- startCodons[ !grepl("Clinical", startCodons$tag), ]
message(sprintf("%i transcripts", nrow(startCodons_filtered)))                          # 19 222 transcripts







