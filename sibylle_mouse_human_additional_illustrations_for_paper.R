######################################################################
### SIBYLLE PAPER - ALL ADDITIONAL FIGURES TO ILLUSTRATE THE PAPER ###
######################################################################


### Libraries
library(tidyverse)
library(GEOquery)
library(openxlsx)
library(ggvenn)


### Parameters and paths are defined in the YAML file
params <- yaml::read_yaml(file = file.path(this.path::here(),
                                           "sibylle_mouse_human_data_preprocessing_parameters.yaml"))
suppressWarnings(params %>% attach) # To enable accessing the parameters by their name without using params$...
tissues <- list()
tissues[["mouse"]] <- tissuesOfInterest$mouse
tissues[["human"]] <- tissuesOfInterest$human
# Set appropriate working directory
setwd(projectPath)




#################
### LOAD DATA ###
#################

### Mouse data

# Metadata
mouseMetadata <- read.table(file.path(projectPath,
                                      dataPath,
                                      mousePath,
                                      "biovacsafe_mouse_metadata_for_analysis.txt"),
                            sep = "\t",
                            header = T,
                            stringsAsFactors = F)

# Count data before low count genes filtering
mouseCount <- lapply(tissues$mouse, function(tiss) {
  x <- read.table(file.path(projectPath,
                            dataPath,
                            mousePath,
                            paste0("biovacsafe_mouse_",
                                   tiss,
                                   "_count_data_for_analysis.txt")),
                  sep = "\t",
                  header = T,
                  stringsAsFactors = F)
})
names(mouseCount) <- tissues$mouse

# Count data after low count genes filtering
mouseCountAfterQC <- lapply(tissues$mouse, function(tiss) {
  x <- read.csv(file.path(projectPath,
                          outputPath,
                          mousePath,
                          "QC",
                          paste0("biovacsafe_mouse_count_data_filtering_low_count_genes_",
                                 tiss,
                                 ".csv")),
                header = T,
                stringsAsFactors = F) %>%
    column_to_rownames(var = "X")
})
names(mouseCountAfterQC) <- tissues$mouse

# Features before low count genes filtering
mouseFeatures <- read.xlsx(file.path(projectPath,
                                      dataPath,
                                      mousePath,
                                      "biovacsafe_mouse_features_protein_coding_without_duplicates_for_analysis.xlsx"),
                            sheet = 1)


### Human data

# Metadata
humanMetadata <- read.table(file.path(projectPath,
                                      dataPath,
                                      humanPath,
                                      "biovacsafe_human_metadata_for_analysis.txt"),
                            sep = "\t",
                            header = T,
                            stringsAsFactors = F)

# Count data before low count genes filtering
humanCount <- lapply(tissues$human, function(tiss) {
  x <- read.table(file.path(projectPath,
                            dataPath,
                            humanPath,
                            paste0("biovacsafe_human_",
                                   tiss,
                                   "_count_data_for_analysis.txt")),
                  sep = "\t",
                  header = T,
                  stringsAsFactors = F)
})
names(humanCount) <- tissues$human

# Count data after low count genes filtering
humanCountAfterQC <- lapply(tissues$human, function(tiss) {
  x <- read.csv(file.path(projectPath,
                          outputPath,
                          humanPath,
                          "QC",
                          paste0("biovacsafe_human_count_data_filtering_low_count_genes_",
                                 tiss,
                                 ".csv")),
                header = T,
                stringsAsFactors = F) %>%
    column_to_rownames(var = "X")
})
names(humanCountAfterQC) <- tissues$human

# Features before low count genes filtering 
humanFeatures <- read.xlsx(file.path(projectPath,
                                      dataPath,
                                      humanPath,
                                      "biovacsafe_human_features_protein_coding_without_duplicates_for_analysis.xlsx"),
                            sheet = 1)


### Mouse-human orthology

mouseHumanOrthologs <- read.csv(file.path(projectPath,
                                          outputPath,
                                          "orthology",
                                          "mouse_human_orthology_only_unique_pairs.csv"),
                                header = T,
                                stringsAsFactors = F)




########################################################################
### INTERSECTION BETWEEN GENES IN MOUSE MUSCLE/BLOOD AND HUMAN BLOOD ###
########################################################################

### Venn diagram of probes between the three compartiments, NOT including low count genes filtering

# List of probes in mouse muscle, with orthologs match and unique ID for future Venn
tmpMM <- mouseFeatures %>%
  dplyr::select(SPOT_ID) %>%
  rename(probeName.mouse = SPOT_ID) %>%
  left_join(mouseHumanOrthologs[,c("probeName.mouse",
                                   "probeName.human")],
            by = "probeName.mouse") %>%
  mutate(customProbeID = paste0("m_",
                                probeName.mouse,
                                "_h_",
                                probeName.human))

# List of probes in mouse blood, with orthologs match and unique ID for future Venn
tmpMB <- mouseFeatures %>%
  dplyr::select(SPOT_ID) %>%
  rename(probeName.mouse = SPOT_ID) %>%
  left_join(mouseHumanOrthologs[,c("probeName.mouse",
                                   "probeName.human")],
            by = "probeName.mouse") %>%
  mutate(customProbeID = paste0("m_",
                                probeName.mouse,
                                "_h_",
                                probeName.human))

# List of probes in human, with orthologs match and unique ID for future Venn
tmpHB <- humanFeatures %>%
  dplyr::select(SPOT_ID) %>%
  rename(probeName.human = SPOT_ID) %>%
  left_join(mouseHumanOrthologs[,c("probeName.mouse",
                                   "probeName.human")],
            by = "probeName.human") %>%
  mutate(customProbeID = paste0("m_",
                                probeName.mouse,
                                "_h_",
                                probeName.human))

# Quick checks on numbers
tmpMM %>%
  filter(!is.na(probeName.human)) %>%
  nrow() # 5371 mouse muscle probes with human blood probes orthologs
tmpMB %>%
  filter(!is.na(probeName.human)) %>%
  nrow()  # 3389 mouse blood probes with human blood probes orthologs
tmpHB %>%
  filter(!is.na(probeName.mouse)) %>%
  nrow() # 4256 human blood probes with mouse probes orthologs
# Intersection between mouse muscle and mouse blood ?
length(intersect(tmpMM$probeName.mouse,
                 tmpMB$probeName.mouse)) # 19266
# Intersection between mouse muscle, mouse blood and human blood
verif <- tmpMM %>%
  filter(probeName.mouse %in% tmpMB$probeName.mouse) %>%
  filter(probeName.human %in% tmpHB$probeName.human) %>%
  filter(!is.na(probeName.human))
nrow(verif) # 3067

# Final data for Venn diagram
dataForVenn <- list(mouseMuscle = tmpMM$customProbeID,
                    mouseBlood = tmpMB$customProbeID,
                    humanBlood = tmpHB$customProbeID)
names(dataForVenn) <- c("mouse muscle",
                        "mouse blood",
                        "human blood")

# Do Venn diagram
gp <- ggvenn(dataForVenn,
             fill_color = c("firebrick", "purple", "#EFC000FF"),
             text_size = 4.5)
ggsave(filename = file.path(projectPath,
                            outputPath,
                            "orthology",
                            "mouse_human_orthology_venn_diagram.png"),
       plot = gp,
       width = 6,
       height = 6)

# Cleaning
rm(tmpMM, tmpMB, tmpHB, dataForVenn, gp)


### Venn diagram of probes between the three compartiments, including low count genes filtering

# List of probes in mouse muscle, with orthologs match and unique ID for future Venn
tmpMM <- data.frame(probeName.mouse = rownames(mouseCountAfterQC$muscle),
                   stringsAsFactors = F) %>%
  left_join(mouseHumanOrthologs[,c("probeName.mouse",
                                   "probeName.human")],
            by = "probeName.mouse") %>%
  mutate(customProbeID = paste0("m_",
                                probeName.mouse,
                                "_h_",
                                probeName.human))

# List of probes in mouse blood, with orthologs match and unique ID for future Venn
tmpMB <- data.frame(probeName.mouse = rownames(mouseCountAfterQC$blood),
                   stringsAsFactors = F) %>%
  left_join(mouseHumanOrthologs[,c("probeName.mouse",
                                   "probeName.human")],
            by = "probeName.mouse") %>%
  mutate(customProbeID = paste0("m_",
                                probeName.mouse,
                                "_h_",
                                probeName.human))

# List of probes in human, with orthologs match and unique ID for future Venn
tmpHB <- data.frame(probeName.human = rownames(humanCount$blood),
                   stringsAsFactors = F) %>%
  left_join(mouseHumanOrthologs[,c("probeName.mouse",
                                   "probeName.human")],
            by = "probeName.human") %>%
  mutate(customProbeID = paste0("m_",
                                probeName.mouse,
                                "_h_",
                                probeName.human))

# Quick checks on numbers
tmpMM %>%
  filter(!is.na(probeName.human)) %>%
  nrow() # 3031 mouse muscle probes with human blood probes orthologs
tmpMB %>%
  filter(!is.na(probeName.human)) %>%
  nrow()  # 3298 mouse blood probes with human blood probes orthologs
tmpHB %>%
  filter(!is.na(probeName.mouse)) %>%
  nrow() # 3281 human blood probes with mouse probes orthologs
# Intersection between mouse muscle and mouse blood ?
length(intersect(tmpMM$probeName.mouse,
                 tmpMB$probeName.mouse)) # 16139
# Intersection between mouse muscle, mouse blood and human blood
verif <- tmpMM %>%
  filter(probeName.mouse %in% tmpMB$probeName.mouse) %>%
  filter(probeName.human %in% tmpHB$probeName.human) %>%
  filter(!is.na(probeName.human))
nrow(verif) # 2139

# Final data for Venn diagram
dataForVenn <- list(mouseMuscle = tmpMM$customProbeID,
                    mouseBlood = tmpMB$customProbeID,
                    humanBlood = tmpHB$customProbeID)
names(dataForVenn) <- c("mouse muscle",
                        "mouse blood",
                        "human blood")

# Do Venn diagram
gp <- ggvenn(dataForVenn,
             fill_color = c("firebrick", "purple", "#EFC000FF"),
             text_size = 4.5)
ggsave(filename = file.path(projectPath,
                            outputPath,
                            "orthology",
                            "mouse_human_orthology_venn_diagram_including_low_count_genes_filtering.png"),
       plot = gp,
       width = 6,
       height = 6)

# Cleaning
rm(tmpMM, tmpMB, tmpHB, dataForVenn, gp)