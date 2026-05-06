###########################################################
### SIBYLLE PAPER - ANALYSIS OF REACTOGENICITY READOUTS ###
###########################################################


### Libraries
library(tidyverse)
library(GEOquery)
library(openxlsx)
library(ggplot2)


### Parameters and paths are defined in the YAML file
params <- yaml::read_yaml(file = file.path(this.path::here(),
                                           "sibylle_mouse_human_data_preprocessing_parameters.yaml"))
suppressWarnings(params %>% attach) # To enable accessing the parameters by their name without using params$...
tissues <- list(mouse = tissuesOfInterest$mouse,
                human = tissuesOfInterest$human)
# Set appropriate working directory
setwd(projectPath)


### List of AEs with classification specific / non-specific
x <- read.table(file.path(projectPath,
                          dataPath,
                          humanPath,
                          "biovacsafe_human_ae_data_including_subjective_with_categories_for_analysis.txt"),
                sep = "\t",
                header = T,
                stringsAsFactors = F) %>%
  dplyr::select(aeTerm, aeCategory, aeCategory2, include) %>%
  distinct() %>%
  rename(isSpecific = include)
write.xlsx(x,
           file.path(projectPath,
                     outputPath,
                     humanPath,
                     "Supplementary_Table_1.xlsx"),
           quote = F,
           rowNames = F)
rm(x)


### Ortholog probes with gene name etc
x <- read.csv(file.path(projectpath,
                        outputPath,
                        "orthology",
                        "mouse_human_orthology_unique_pairs.csv"),
              header = T,
              stringsAsFactors = F) %>%
  dplyr::select(probeName.mouse,
                geneSymbol.mouse,
                ensemblTranscriptID.mouse,
                refseq.mouse,
                probeName.human,
                geneSymbol.human,
                ensemblTranscriptID.human,
                refseq.human)
write.xlsx(x,
           file.path(projectPath,
                     outputPath,
                     "orthology",
                     "Supplementary_Table_4.xlsx"),
           quote = F,
           rowNames = F)
rm(x)


