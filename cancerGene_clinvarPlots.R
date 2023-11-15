# Graphs and figures

# Remove objects in workspace
rm(list = ls())

# Required packages
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# -------------------------------------------------------------------
inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/4A/BIOL 469/Final Project/BIOL469/Data/Input"
inDir  <- sprintf("%s/Input", dirname(inpath))
setwd(inDir)

# Read in dataset
cancerGenes_olapFinal <- read_delim("cancerGenes_olapFinal.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# -------------------------------------------------------------------
# Renaming values in the clinicalSignificance column
cancerGenes_olapFinal <- cancerGenes_olapFinal %>%
  mutate(clinicalSignificance = case_when(
    clinicalSignificance == "Likely_benign" ~ "Likely Benign",
    clinicalSignificance == "Benign/Likely_benign" ~ "Benign/Likely Benign",
    clinicalSignificance == "Likely_pathogenic" ~ "Likely Pathogenic",
    clinicalSignificance == "Pathogenic/Likely_pathogenic" ~ "Pathogenic/Likely Pathogenic",
    TRUE ~ clinicalSignificance  # Keeps all other values as they are
  ))

# Calculate frequencies
cancerGenes_olapFinal$Count <- 1
aggregateData <- aggregate(Count ~ generalOntology + clinicalSignificance + cancer, 
                           data = cancerGenes_olapFinal, FUN = sum)
# Convert factors
aggregateData$generalOntology      <- factor(aggregateData$generalOntology)
aggregateData$clinicalSignificance <- factor(aggregateData$clinicalSignificance)
aggregateData$cancer               <- factor(aggregateData$cancer)

# Define fixed positions for each clinicalSignificance level
significance_levels <- c("Benign", "Benign/Likely Benign", "Likely Benign", 
                         "Likely Pathogenic", "Pathogenic/Likely Pathogenic", "Pathogenic")
# Increment for spacing
increment = 0.5
significance_positions <- setNames(seq(from = 1, by = increment, 
                                       length.out = length(significance_levels)), significance_levels)

# Add fixed positions to the data
adjustedData <- aggregateData %>%
  mutate(clinicalSignificance = factor(clinicalSignificance, levels = significance_levels)) %>%  # Set factor levels
  mutate(Position = significance_positions[clinicalSignificance]) %>%
  group_by(generalOntology, clinicalSignificance) %>%
  mutate(Position = ifelse(is.na(Position), max(Position, na.rm = TRUE) + increment, Position)) %>%
  ungroup()

# Reorder the levels of generalOntology
adjustedData$generalOntology <- factor(adjustedData$generalOntology, 
                                       levels = c("Other", "UTR Variant", "Upstream Variant",
                                                  "Intron Variant", "Splice Variant",
                                                  "Synonymous Variant", "Inframe Variant",
                                                  "Frameshift Variant", "Missense/Nonsense Variant"))

# Create the plot
ggplot(adjustedData, aes(x = Position, y = generalOntology, size = Count, color = clinicalSignificance)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ cancer, scales = "free_x", ncol = 5, strip.position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1),
        # Remove vertical gridlines
        panel.grid.major.x = element_blank(),  
        panel.grid.minor.x = element_blank(),
        # Reduce vertical and horizontal spacing
        panel.spacing = unit(0.5, "lines"),
        # Light border around each panel
        panel.border = element_rect(colour = "grey", fill = NA, linewidth = 0.5)) +  
  labs(y = "Mutation Class", size = "Frequency", color = "Clinical Significance") +
  scale_color_manual(values = c("Benign" = "#77DA94", "Benign/Likely Benign" = "#80DFE1", 
                                "Likely Benign" = "#AFCCFF", "Likely Pathogenic" = "#FAB0F1",
                                "Pathogenic/Likely Pathogenic" = "#D8CC77", "Pathogenic" = "#FAAEA8"))
