# Preliminary analysis of genes impacting cancers
# Determine canonical transcripts of implicated cancer genes
# Calculate promoter positions of these transcripts

# Remove objects in workspace
rm(list = ls())

# Required packages
suppressMessages(library(rtracklayer))
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))

# -------------------------------------------------------------------
# Import MANE data
# Downloaded May 16, 2023 (Release 1.1, hg38 genome assembly)
# https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/

inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/4A/BIOL 469/Final Project/BIOL469/Data/Input"
inDir  <- sprintf("%s/Input", dirname(inpath))
setwd(inDir)

MANE   <- readGFF("MANE.GRCh38.v1.1.ensembl_genomic.gff.gz")
MANEdf <- as.data.frame(MANE)
message(sprintf("Loaded %i transcripts", nrow(MANEdf)))                                 # 521 393 transcripts

cancerGenes <- read_excel("Cancer Genes.xlsx")
cancerGenes <- as.data.frame(cancerGenes)
cancerGenes <- cancerGenes%>% rename("transcript_id" = "Ensembl Transcript ID")

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

# Merge cancer genes with MANE data
cancerGenes_MANE <- left_join(cancerGenes, startCodons_filtered, by = "transcript_id")

# -------------------------------------------------------------------
# PROMOTER POSITIONS
chr <- cancerGenes_MANE$seqid

genomicRanges <- GRanges(chr,
                         IRanges(cancerGenes_MANE$start, cancerGenes_MANE$end),
                         strand = cancerGenes_MANE$strand,
                         gene_name = cancerGenes_MANE$`HGNC Symbol`,
                         transcript_id = cancerGenes_MANE$transcript_id)

cancerGenes_promoters <- promoters(genomicRanges, upstream = 1000, downstream = 100)
cancerGenes_promoters <- as.data.frame(cancerGenes_promoters)




