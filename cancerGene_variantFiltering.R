# Preliminary analysis of genes impacting cancers
# Determine canonical transcripts of implicated cancer genes
# Calculate promoter positions of these transcripts
# Identify ClinVar variants overlapping promoter regions

# Remove objects in workspace
rm(list = ls())

# Required packages
suppressMessages(library(rtracklayer))
suppressMessages(library(readxl))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(vcfR))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(caret))
suppressMessages(library(glmnet)) # For logistic regression model

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

# Import genes implicated in cancer (breast, bladder, colorectal, endometrial, kidney)
cancerGenes <- read_excel("Cancer Genes.xlsx")
cancerGenes <- as.data.frame(cancerGenes)
cancerGenes <- cancerGenes %>% rename("transcript_id" = "Ensembl Transcript ID")

# Read in green genes from NHS Genomics England PanelApp
# Pulled from: https://nhsgms-panelapp.genomicsengland.co.uk/entities
# Current as of November 10, 2023 (3640 genes)
greenGenes_NHS <- read_table(file = sprintf("%s/greenGenes_NHS.txt", inDir), col_names = FALSE)
greenGenes_NHS <- greenGenes_NHS %>% gather(value = "gene_name") %>% select(-key)
excludedNHS_genes <- data.frame(gene_name = 
                                  cancerGenes$`HGNC Symbol`[!cancerGenes$`HGNC Symbol` %in% greenGenes_NHS$gene_name])

# Import clinVar VCF
# Downloaded November 10, 2023 (hg38 genome assembly)
# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
clinVar_VCF <- read.vcfR("clinvar.vcf.gz")                                              # 2 300 621 variants

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
                         transcript_id = cancerGenes_MANE$transcript_id,
                         cancer = cancerGenes_MANE$Cancer)
genomicRanges <- genomicRanges[ !duplicated(genomicRanges$transcript_id) ]

# GRanges object for promoter regions
cancerGenes_promotersGR <- promoters(genomicRanges, upstream = 1000, downstream = 100)
cancerGenes_promoters   <- as.data.frame(cancerGenes_promotersGR)

# -------------------------------------------------------------------
# OVERLAPPING CLINVAR VARIANTS

# Extract necessary columns from clinVar_VCF
chrom <- clinVar_VCF@fix[, "CHROM"]
pos   <- as.integer(clinVar_VCF@fix[, "POS"])
ref   <- clinVar_VCF@fix[, "REF"]

# Convert clinVar data to GRanges
clinVar_GR <- GRanges(seqnames = as.factor(chrom), 
                      ranges = IRanges(start = pos, end = pos + nchar(as.character(ref)) - 1))

# Prepend 'chr' to chromosome names in clinVar_GR
seqlevels(clinVar_GR) <- paste0("chr", seqlevels(clinVar_GR))

# Extract overlapping variants (query is promoters, subject is clinVar variants)
cancerGenes_promoterOlap   <- findOverlaps(cancerGenes_promotersGR, clinVar_GR)

# Extract indices of overlapping variants and promoters
queryIndices   <- queryHits(cancerGenes_promoterOlap)
subjectIndices <- subjectHits(cancerGenes_promoterOlap)

# Extract relevant rows from clinVar_GR using subjectIndices
clinVarRelevant <- data.frame(
  seqnames     = seqnames(clinVar_GR)[subjectIndices],
  variantStart = start(clinVar_GR)[subjectIndices],
  variantEnd   = end(clinVar_GR)[subjectIndices],
  width        = width(clinVar_GR)[subjectIndices],
  #strand       = strand(clinVar_GR)[subjectIndices],
  INFO         = clinVar_VCF@fix[subjectIndices, "INFO"]
)

# Extract relevant rows from cancerGenes_promoters using queryIndices
promoterRelevant <- cancerGenes_promoters[queryIndices, ]
promoterRelevant <- promoterRelevant %>% rename("promoterStart" = "start")
promoterRelevant <- promoterRelevant %>% rename("promoterEnd" = "end")
promoterRelevant <- promoterRelevant %>% rename("promoterWidth" = "width")
promoterRelevant <- subset(promoterRelevant, select = -seqnames)

# Merge the two dataframes (5903 variants)
cancerGenes_olapVariants <- cbind(clinVarRelevant, promoterRelevant)

# Convert INFO column to character type
cancerGenes_olapVariants$INFO <- as.character(cancerGenes_olapVariants$INFO)

# Parse INFO column and extract clinical significance
extractCLNSIG <- function(infoString) {
  clnsigPattern <- "CLNSIG=[^;]+"
  matches <- gregexpr(clnsigPattern, infoString)
  clnsigValues <- regmatches(infoString, matches)
  
  if(length(clnsigValues[[1]]) > 0){
    # Join multiple CLNSIG values with a delimiter
    return(paste(sapply(clnsigValues[[1]], function(x) sub("CLNSIG=", "", x)), collapse = ";"))
  } 
  else{
    # Return "NA" or "None" if no CLNSIG is found
    return("NA")
  }
}

# Parse INFO column and extract molecular consequence
extractInfoField <- function(infoString, fieldName){
  # Split the INFO string into key-value pairs
  keyValuePairs <- strsplit(infoString, ";")[[1]]
  
  # Find the key-value pair that starts with the desired field name
  fieldPair <- keyValuePairs[grep(paste0("^", fieldName), keyValuePairs)]
  
  # Extract the value part of the key-value pair
  if(length(fieldPair) > 0){
    fieldValue <- sub(paste0(fieldName, "="), "", fieldPair)
  } 
  else{
    fieldValue <- NA  # Return NA if the field is not found
  }
  return(fieldValue)
}

# Extract clinical significance and molecular consequence
cancerGenes_olapVariants$clinicalSignificance <- as.factor(sapply(cancerGenes_olapVariants$INFO, extractCLNSIG))
cancerGenes_olapVariants$molecularConsequence <- as.factor(sapply(cancerGenes_olapVariants$INFO, extractInfoField, "MC"))

# -------------------------------------------------------------------
# PROCESS CLINICAL VARIANTS
# Clinical Significance
levels(cancerGenes_olapVariants$clinicalSignificance)

# Drop rows with ambiguous/unprovided significance (2369 variants)
exclude_clnicalSignificance <- c("Conflicting_interpretations_of_pathogenicity", 
                                 "not_provided","Uncertain_significance", "other")

cancerGenes_olapVariants <- cancerGenes_olapVariants[!cancerGenes_olapVariants$clinicalSignificance 
                                                     %in% exclude_clnicalSignificance, ]

# Check
cancerGenes_olapVariants$clinicalSignificance <- droplevels(cancerGenes_olapVariants$clinicalSignificance)
levels(cancerGenes_olapVariants$clinicalSignificance)

# Reclassify into either Benign/Likely Benign or Pathogenic/Likely Pathogenic
cancerGenes_olapVariants$reclassifiedSignificance <- ifelse(
  cancerGenes_olapVariants$clinicalSignificance %in% c("Benign", "Benign/Likely_benign", "Likely_benign"),
  "Benign/Likely_benign",
  "Pathogenic/Likely_pathogenic"
)

# Check new levels
table(cancerGenes_olapVariants$reclassifiedSignificance)
table(cancerGenes_olapVariants$clinicalSignificance)

# Molecular Consequence
levels(cancerGenes_olapVariants$molecularConsequence)

# Define a mapping function to reclassify based on sequence ontology
map_toCategory <- function(so_term){
  if(grepl("splice", so_term)){
    return("Splice Variant")
  } 
  else if(grepl("missense|nonsense", so_term)){
    return("Missense/Nonsense Variant")
  }
  else if(grepl("synonymous", so_term)){
    return("Synonymous Variant")
  }
  else if(grepl("frameshift", so_term)){
    return("Frameshift Variant")
  }
  else if(grepl("inframe", so_term)){
    return("Inframe Variant")
  }
  else if(grepl("intron", so_term)){
    return("Intron Variant")
  }
  else if(grepl("UTR", so_term)){
    return("UTR Variant")
  } 
  else if(grepl("upstream", so_term)){
    return("Upstream Variant")
  }
  else{
    return("Other")
  }
}

cancerGenes_olapVariants$generalOntology <- sapply(cancerGenes_olapVariants$molecularConsequence, map_toCategory)

# Check the new levels
table(cancerGenes_olapVariants$generalOntology)

# Extract unique levels
uniqueLevels <- unique(cancerGenes_olapVariants$molecularConsequence)

# Apply the mapping function to each unique level
mappedCategories <- sapply(uniqueLevels, map_toCategory)

# Identify levels that are mapped to "Other"
otherLevels <- uniqueLevels[mappedCategories == "Other"]

# Print the levels that fall under "Other"
print(otherLevels)

# -------------------------------------------------------------------
# CLEAN DATASET
cancerGenes_olapFinal <- subset(cancerGenes_olapVariants, select = -INFO)

# Output directory
outpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/4A/BIOL 469/Final Project/BIOL469/Data/"
outDir  <- sprintf("%s/Data/", dirname(outpath))
if(!file.exists(outDir)) dir.create(outDir)

write.table(cancerGenes_olapFinal, file = sprintf("%s/cancerGenes_olapFinal.tsv", outDir),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)