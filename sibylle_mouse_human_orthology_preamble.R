##########################################################################################
### PREPARING NECESSARY DATA TO RETRIEVE TRANSCRIPTS ORTHOLOGY BETWEEN MOUSE AND HUMAN ###
##########################################################################################




# Orthology will be based on chromosomal location, using UCSC LiftOver platform
# => see https://genome.ucsc.edu/cgi-bin/hgLiftOver





#### Library load
library(tidyverse)
library(GEOquery)
library(openxlsx)


### Parameters and paths are defined in the YAML file
params <- yaml::read_yaml(file = file.path(this.path::here(),
                                           "sibylle_mouse_human_data_preprocessing_parameters.yaml"))
suppressWarnings(params %>% attach) # To enable accessing the parameters by their name without using params$...
tissues <- tissuesOfInterest$human
# Set appropriate working directory
setwd(projectPath)




### Load features information for mouse and human study

# Mouse
load(file.path(projectPath,
               dataPath,
               mousePath,
               "biovacsafe_mouse_gse120661_processed_from_GEO.RData"))
mouseFeatures <- fData(gset)
rm(gset)

# Human
load(file.path(projectPath,
               dataPath,
               humanPath,
               "biovacsafe_human_gse124533_processed_from_GEO.RData"))
humanFeatures <- fData(gset)
rm(gset)




### Write list of chromosomal locations for mouse to be used in LiftOver tool

# Mouse
toWrite <- mouseFeatures %>%
  dplyr::select(SPOT_ID,
                CHROMOSOMAL_LOCATION,
                GENE_SYMBOL,
                ENSEMBL_ID,
                REFSEQ)
toWrite$chr <- str_split(toWrite$CHROMOSOMAL_LOCATION,
                         pattern = ":",
                         simplify = T)[,1]
toWrite$start <- str_split(str_split(toWrite$CHROMOSOMAL_LOCATION,
                                     pattern = ":",
                                     simplify = T)[,2],
                           pattern = "-",
                           simplify = T)[,1]
toWrite$end <- str_split(str_split(toWrite$CHROMOSOMAL_LOCATION,
                                   pattern = ":",
                                   simplify = T)[,2],
                         pattern = "-",
                         simplify = T)[,2]
toWrite <- toWrite %>%
  mutate(start = as.integer(start),
         end = as.integer(end)) %>%
  mutate(positionForLiftOver = case_when(
    start < end ~ paste(chr, start - 1, end),
    start >= end ~ paste(chr, end - 1, start)
  )) # For some reason in LiftOver all start positions in matches are shifted by 1 compared to initial start value
write.table(toWrite,
            file.path(projectPath,
                      dataPath,
                      mousePath,
                      "biovacsafe_mouse_features_chromosomal_location_with_probe_ids.txt"),
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)
toWrite <- toWrite %>%
  dplyr::select(positionForLiftOver) %>%
  filter(!is.na(positionForLiftOver))
write.table(toWrite,
            file.path(projectPath,
                      dataPath,
                      mousePath,
                      "biovacsafe_mouse_features_chromosomal_location_without_probe_ids.txt"),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)
rm(toWrite)


# Human
toWrite <- humanFeatures %>%
  dplyr::select(SPOT_ID,
                CHROMOSOMAL_LOCATION,
                GENE_SYMBOL,
                ENSEMBL_ID,
                REFSEQ) %>%
  distinct() # Exclude duplicated probes
toWrite$chr <- str_split(toWrite$CHROMOSOMAL_LOCATION,
                         pattern = ":",
                         simplify = T)[,1]
toWrite$start <- str_split(str_split(toWrite$CHROMOSOMAL_LOCATION,
                                     pattern = ":",
                                     simplify = T)[,2],
                           pattern = "-",
                           simplify = T)[,1]
toWrite$end <- str_split(str_split(toWrite$CHROMOSOMAL_LOCATION,
                                   pattern = ":",
                                   simplify = T)[,2],
                         pattern = "-",
                         simplify = T)[,2]
toWrite <- toWrite %>%
  mutate(start = as.integer(start),
         end = as.integer(end)) %>%
  mutate(positionForLiftOver = case_when(
    start < end ~ paste(chr, start - 1, end),
    start >= end ~ paste(chr, end - 1, start)
  )) # For some reason in LiftOver all start positions in matches are shifted by 1 compared to initial start value
write.table(toWrite,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_features_chromosomal_location_with_probe_ids.txt"),
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)
toWrite <- toWrite %>%
  dplyr::select(positionForLiftOver) %>%
  filter(!is.na(positionForLiftOver))
write.table(toWrite,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_features_chromosomal_location_without_probe_ids.txt"),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)
rm(toWrite)


### Cleaning
rm(list = ls())
gc()




### Genome assemply used for the BioVacSafe micro-arrays ?
# => based on Agilent website informations
# - mouse : RefSeq37 / NCBI37 (https://www.agilent.com/store/fr_FR/Prod-G4852A/G4852A)
# - human : hg18 or hg 19 ? (https://www.agilent.com/cs/library/brochures/5990-3368en_lo.pdf)
# => after verification of a few probes, it seems more to be hg19 than hg18 (contrary to what is written in the doc)
