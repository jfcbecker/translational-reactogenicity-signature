###########################################################
### BIOVACSAFE MOUSE DATA - DATA LOAD AND PREPROCESSING ###
###########################################################




### Libraries
library(tidyverse)
library(GEOquery)
library(openxlsx)


### Parameters and paths are defined in the YAML file
params <- yaml::read_yaml(file = file.path(this.path::here(),
                                           "sibylle_mouse_human_data_preprocessing_parameters.yaml"))
suppressWarnings(params %>% attach) # To enable accessing the parameters by their name without using params$...
tissues <- tissuesOfInterest$mouse
# Set appropriate working directory
setwd(projectPath)




#########################
### INITIAL DATA LOAD ###
#########################

### Count data - load from GEO (if necessary)
mouseGeoDataFile <- file.path(projectPath,
                              dataPath,
                              mousePath,
                              "biovacsafe_mouse_gse120661_processed_from_GEO.RData")
if (!(file.exists(mouseGeoDataFile)))
{
  Sys.setenv(VROOM_CONNECTION_SIZE = 10000072) # otherwise error 'The size of the connection buffer (131072) was not large enough'
  
  # Retrieving data from GEO (using R script from website tutorial)
  geoID <- "GSE120661"
  gset <- getGEO(geoID,
                 GSEMatrix = T, # FALSE to get the full data from SOFT files (not only the count table) => takes longer, but more complete (though Agilent flags are not there)
                 getGPL = T, # To get features / probes / genes annotations
                 parseCharacteristics = T) 
  if (length(gset) > 1)
  {
    idx <- grep("GPL10787", attr(gset, "names"))
  } else
  {
    idx <- 1
  }
  gset <- gset[[idx]]
  rm(idx)
  
  # Save full loaded object from GEO
  save(gset,
       file = file.path(projectPath,
                        dataPath,
                        mousePath,
                        "biovacsafe_mouse_gse120661_processed_from_GEO.RData"))
} else
{
  load(file.path(projectPath,
                 dataPath,
                 mousePath,
                 "biovacsafe_mouse_gse120661_processed_from_GEO.RData"))
}
rm(mouseGeoDataFile)

# Get count data
countData <- exprs(gset)
# => NB: This is background corrected and normalised data (see Weiner et al. 2019)

# Get metadata
metadata <- pData(gset)

# Get features / probes info
features <- fData(gset)




######################
### DATA WRANGLING ###
######################

### Formatting features data

# Combine probe name with gene name to create featureID
features <- features %>%
  mutate(featureID = if_else(GENE_SYMBOL == "",
                             SPOT_ID,
                             paste0(GENE_SYMBOL,
                                    "_",
                                    SPOT_ID))) %>%
  dplyr::select(featureID, ID, SPOT_ID, GENE_SYMBOL, everything())

  
### Formatting metadata

# Focusing on most important information
metadataFull <- metadata
clnms <- colnames(metadata)
charaClnms <- clnms[grepl("characteristics", clnms)]
metadata <- metadata %>%
  dplyr::select(all_of(c("title", "geo_accession", charaClnms))) %>%
  dplyr::rename(sampleID = geo_accession)

# Proper naming and formatting of characteristics
for (cl in charaClnms)
{
  vals <- metadata[,cl]
  if (sum(vals == "") != length(vals)) # Non empty columns
  {
    tmp <- str_split(vals, pattern = ": ", simplify = T)
    newvals <- tmp[,2]
    newcl <- unique(tmp[,1])
    metadata <- metadata %>%
      mutate(newcol = newvals) %>%
      rename(!!newcl := newcol) %>%
      dplyr::select(-all_of(cl))
  } else
  {
    metadata <- metadata %>%
      dplyr::select(-all_of(cl))
  }
  rm(vals, tmp, newvals, newcl)
}
rm(cl)
metadata <- metadata %>%
  dplyr::rename("animal" = "animal replicate",
                "strain" = "strain/background")

# Definition of additional / corrected variables
tps <- paste0(sort(unique(as.integer(metadata$timepoint))), "h")
metadata <- metadata %>%
  # Proper naming of tissue
  mutate(tissue = case_when(
    tissue == "muscle" ~ "muscle",
    tissue == "PBMC" ~ "blood",
    tissue == "draining medial iliac lymph nodes (MLN)" ~ "MLN")) %>%
  mutate(tissueAcronym = case_when(
    tissue == "blood" ~ "BL",
    tissue == "muscle" ~ "MU",
    tissue == "MLN" ~ "MLN")) %>%
  # Proper naming of treatments
  mutate(treatment = Hmisc::capitalize(treatment)) %>%
  mutate(treatmentAcronym = case_when(
    treatment %in% c("PolyIC", "IFA", "LPS") ~ treatment,
    treatment == "Engerix B" ~ "ENG",
    treatment == "Pentavac SD" ~ "PERT",
    treatment == "Saline" ~ "SAL",
    treatment == "Agrippal" ~ "TriFLU",
    treatment == "Fluad" ~ "TriFLU.MF59"
  )) %>%
  # Reactogenicity class
  mutate(reactoClass = case_when(
    treatment %in% lowReactoClass ~ "low",
    treatment %in% mediumReactoClass ~ "medium",
    treatment %in% highReactoClass ~ "high"
  )) %>%
  mutate(reactoClass = factor(reactoClass,
                              levels = c("low", "medium", "high"),
                              ordered = T)) %>%
  # Reshaping timepoint
  mutate(timepoint = factor(paste0(timepoint, "h"),
                            levels = tps)) %>%
  # Proper animal ID
  mutate(animalID = paste0(treatmentAcronym,
                            "_",
                            timepoint,
                            "_",
                            animal)) %>%
  # Define new sample ID based on tissue, treatment, timepoint and and animal ID
  mutate(originalSampleID = sampleID) %>%
  mutate(sampleID = paste0(treatmentAcronym,
                           "_",
                           tissueAcronym,
                           "_",
                           timepoint,
                           "_",
                           animal))


### Formatting count data

# Arrange samples by tissue type, treatment and timepoint in count data and metadata
metadata <- metadata %>%
  arrange(tissue, treatment, timepoint)
countData <- countData[, metadata$originalSampleID]

# Arrange features in same order in count data and feature list
features <- features %>%
  arrange(SPOT_ID)
countData <- countData[features$SPOT_ID,]

# Adding features ID to data frame
originalCountData <- countData
countData <- countData %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID")
# Correct colname = sample ID with metadata included
colnames(countData) <- c("ID", metadata$sampleID)




#####################
### DATA CLEANING ###
#####################

### Exclude useless probes for analysis

# Features to keep, based on information from GEO
featsToKeep <- features %>%
  # Replace "" by NA
  mutate(across(.cols = everything(),
                .fns = ~ if_else(.x == "", NA_character_, as.character(.x)))) %>%
  # Remove control probes
  filter(CONTROL_TYPE == F) %>%
  # Remove probes with no gene symbol or ENSEMBL ID
  filter(!is.na(GENE_SYMBOL)
         | !is.na(ENSEMBL_ID))

# Features to keep, without duplicates
featsToKeepNoDuplik <- featsToKeep %>%
  dplyr::select(SPOT_ID, GENE_SYMBOL, GENE_NAME, ENSEMBL_ID, REFSEQ) %>%
  distinct()

# Features to keep, based on biotype from Biomart
biomartMatch <- read.table(file.path(projectPath,
                                     dataPath,
                                     mousePath,
                                     "biomart_match_mouse_features_agilent_8x60K_v1_v2.txt"),
                           sep = ",",
                           header = T,
                           stringsAsFactors = F)
biomartMatch <- biomartMatch %>%
  rename(geneID = Gene.stable.ID,
         transcriptID = Transcript.stable.ID,
         geneType = Gene.type,
         transcriptType = Transcript.type,
         spotIDv1 = AGILENT.SurePrint.G3.GE.8x60k.probe,
         spotIDv2 = AGILENT.SurePrint.G3.GE.8x60k.v2.probe)

# Match on ENSEMBL
tomatch <- featsToKeepNoDuplik
tmp <- tomatch %>%
  filter(!is.na(ENSEMBL_ID)) %>%
  rename(transcriptID = ENSEMBL_ID)
tmp2 <- biomartMatch %>%
  filter(!is.na(transcriptID)) %>%
  filter(transcriptType == "protein_coding") %>%
  dplyr::select(transcriptID, transcriptType) %>%
  distinct()
matchEnsembl <- tmp %>%
  left_join(tmp2, by = c("transcriptID")) %>%
  rename(ENSEMBL_ID = transcriptID)
nrow(matchEnsembl) == nrow(tmp) # Check if not introducing duplicates
rm(tmp, tmp2)

# Exclude non proteing coding transcripts
proteinCodingFeats <- matchEnsembl %>%
  filter(!is.na(transcriptType)
         & transcriptType == "protein_coding") %>%
  dplyr::select(SPOT_ID, transcriptType)
featsToKeepNoDuplikProtCoding <- featsToKeepNoDuplik %>%
  inner_join(proteinCodingFeats, by = "SPOT_ID")

# Summary table of filtering steps
summaryNbFeats <- data.frame(filteringSteps = c("initial",
                                                "excluding controls and missing symbol or ensembl or entrez ID",
                                                "excluding duplicates",
                                                "excluding non protein coding features"),
                             n = c(nrow(features),
                                   nrow(featsToKeep),
                                   nrow(featsToKeepNoDuplik),
                                   nrow(featsToKeepNoDuplikProtCoding)),
                             stringsAsFactors = F)

# Filter in count data
countDataFiltered <- countData %>%
  filter(ID %in% featsToKeep$ID) %>%
  left_join(featsToKeep[, c("ID", "SPOT_ID")], by = "ID") %>%
  dplyr::select(ID, SPOT_ID, everything()) %>%
  filter(SPOT_ID %in% featsToKeepNoDuplikProtCoding$SPOT_ID)


### Average gene expression values for each probe on all replicates
dataAllTissues <- lapply(unique(metadata$tissue), function(tiss) {
  curMeta <- metadata %>%
    filter(tissue == tiss)
  curCount <- countDataFiltered[,c("ID", "SPOT_ID", curMeta$sampleID)]
  tmp <- curCount %>%
    dplyr::select(-all_of(c("ID", "SPOT_ID")))
  ids <- curCount$SPOT_ID
  res <- limma::avereps(tmp, ID = ids)
  # => output = matrix with rows = unique probes (rownames = SPOT_ID already set), cols = samples
  return(res)
})
names(dataAllTissues) <- unique(metadata$tissue)




##############################
### SAVING DATA BEFORE QCs ###
##############################

# All in one RData file
save(dataAllTissues,
     metadata,
     featsToKeep,
     featsToKeepNoDuplik,
     featsToKeepNoDuplikProtCoding,
     summaryNbFeats,
     file = file.path(projectPath,
                      dataPath,
                      mousePath,
                      "biovacsafe_mouse_data_for_analysis.RData"))

# Counts (filtered -> no control probes left)
for (tiss in unique(metadata$tissue))
{
  write.table(dataAllTissues[[tiss]],
              file.path(projectPath, 
                        dataPath,
                        mousePath,
                        paste0("biovacsafe_mouse_",
                               tiss,
                               "_count_data_for_analysis.txt")),
              sep = "\t",
              quote = F,
              row.names = T)
}
rm(tiss)


# Metadata
write.table(metadata,
            file.path(projectPath,
                      dataPath,
                      mousePath,
                      "biovacsafe_mouse_metadata_for_analysis.txt"),
            sep = "\t",
            quote = F,
            row.names = F)

# Features / mapping info (filtered -> no control probes left)
write.xlsx(featsToKeepNoDuplik,
            file.path(projectPath,
                      dataPath,
                      mousePath,
                      "biovacsafe_mouse_features_without_duplicates_for_analysis.xlsx"),
            quote = F)
write.xlsx(featsToKeepNoDuplikProtCoding,
            file.path(projectPath,
                      dataPath,
                      mousePath,
                      "biovacsafe_mouse_features_protein_coding_without_duplicates_for_analysis.xlsx"),
            quote = F)

# Summary of number of features across filtering
write.table(summaryNbFeats,
            file.path(projectPath,
                      dataPath,
                      mousePath,
                      "biovacsafe_mouse_features_filtering_summary_table.txt"),
            sep = "\t",
            quote = F,
            row.names = F)




#################################
### LOW COUNT GENES FILTERING ###
#################################

### Create appropriate directory for QCs and corrections
dir.create(file.path(projectPath,
                     outputPath,
                     mousePath),
           showWarnings = F)

dir.create(file.path(projectPath,
                     outputPath,
                     mousePath,
                     "QC"),
           showWarnings = F)


### Preliminary investigation to find the appropriate threshold for each tissue type
for (tiss in tissues)
{
  # Select data for considered tissues
  tmpData <- matrix(dataAllTissues[[tiss]])
  # Calculate gene-wise median
  tmpMedian <- matrixStats::rowMedians(tmpData)
  # Histogram of medians
  png(file.path(projectPath,
                outputPath,
                mousePath,
                "QC",
                paste0("gene_wise_median_",
                       tiss,
                       "_histogram.png")),
      width = 640,
      height = 480)
  hist(tmpMedian,
       breaks = 100,
       main = paste0("Gene-wise median of normalised counts - ", tiss),
       xlab = "Normalised count value",
       ylab = "Frequency")
  abline(v = lowCountThreshold[["mouse"]][[tiss]], col = "red")
  dev.off()
  # Cleaning
  rm(tmpData, tmpMedian)
}
rm(tiss)


### Perform filtering for each tissue type
# Criterion 1: keep only genes with > X counts in at least K samples, with K = nb samples in smallest tissue-TP-treatment group
# => see https://f1000research.com/articles/5-1384"
# Criterion 2: keep only genes with sum of counts across all samples > X*K with X the low count threshold and K the nb of samples in smallest tissue-TP-treatment group
# => HERE WE USE CRITERION 2
genesToKeepPerTissue <- lapply(tissues, function(tiss) {
  # Number of samples per treatment - timepoint group
  tmp <- metadata %>%
    filter(tissue == tiss) %>%
    group_by(treatment, timepoint) %>%
    dplyr::summarise(n = n())
  K <- floor(min(tmp$n))
  # List of genes passing the filtering criterion
  genesToKeep1 <- rownames(dataAllTissues[[tiss]])[rowSums(dataAllTissues[[tiss]] > lowCountThreshold[["mouse"]][[tiss]]) >= K]
  genesToKeep2 <- rownames(dataAllTissues[[tiss]])[rowSums(dataAllTissues[[tiss]]) > lowCountThreshold[["mouse"]][[tiss]]*K]
  res <- list(lowCountThreshold = lowCountThreshold[["mouse"]][[tiss]],
              minNbSamplesPassingThreshold = K,
              genesToKeep = genesToKeep2)
  return(res)
})
names(genesToKeepPerTissue) <- tissues


## Filter low count genes
# All data
dataFilteredPerTissue <- lapply(tissues, function(tiss) {
  # Genes to keep for considered tissue
  genesToKeep <- genesToKeepPerTissue[[tiss]]$genesToKeep
  # Filtered data
  dat <- dataAllTissues[[tiss]][genesToKeep,]
  return(dat)
})
names(dataFilteredPerTissue) <- tissues
# Save count data excluding low count genes before filtering outliers
for (tiss in names(dataFilteredPerTissue))
{
  write.csv(dataFilteredPerTissue[[tiss]],
            file.path(projectPath,
                      outputPath,
                      mousePath,
                      "QC",
                      paste0("biovacsafe_mouse_count_data_filtering_low_count_genes_",
                             tiss,
                             ".csv")),
            quote = F,
            row.names = T)
}
rm(tiss)




####################
### QUICK CHECKS ###
####################

runQuickChecks <- F

if (runQuickChecks)
{
  ### Investigate Engerix B samples at D0 (muscle)
  
  # Samples of interest
  sampsToCheck <- metadata %>%
    filter(tissue == "muscle"
           & treatment == "Engerix B"
           & timepoint == "0h") %>%
    dplyr::select(sampleID,
                  originalSampleID,
                  tissue,
                  treatment,
                  timepoint,
                  everything())
  
  # Data before avereps
  tmp1 <- countDataFiltered[,c("SPOT_ID", sampsToCheck$sampleID)] %>%
    pivot_longer(!SPOT_ID,
                 names_to = "sampleID",
                 values_to = "countStep1")
  # Data after avereps but before low count genes filtering
  tmp2 <- dataAllTissues[["muscle"]][,sampsToCheck$sampleID] %>%
    as.data.frame() %>%
    rownames_to_column(var = "SPOT_ID") %>%
    pivot_longer(!SPOT_ID,
                 names_to = "sampleID",
                 values_to = "countStep2")
  # Data after low count genes filtering
  tmp3 <- dataFilteredPerTissue[["muscle"]][,sampsToCheck$sampleID] %>%
    as.data.frame() %>%
    rownames_to_column(var = "SPOT_ID") %>%
    pivot_longer(!SPOT_ID,
                 names_to = "sampleID",
                 values_to = "countStep3")
  
  # Merge the three data types together
  dataToCheck <- tmp1 %>%
    full_join(tmp2,
              by = c("SPOT_ID",
                     "sampleID")) %>%
    full_join(tmp3,
              by = c("SPOT_ID",
                     "sampleID"))
  rm(tmp1, tmp2, tmp3)
  
  # Compare boxplot of gene expression per sample at each step
  figs <- lapply((1:3), function(stepNb) {
    gp <- ggplot(data = dataToCheck) +
      geom_boxplot(aes(x = sampleID,
                       y = get(paste0("countStep", stepNb))),
                   fill = "dodgerblue") +
      theme_bw(base_size = 14) +
      labs(x = "Sample",
           y = "Normalised count",
           title = paste0("Step ", stepNb),
           fill = "") +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5,
                                       size = 5,
                                       face = "bold"))
    return(gp)
  })
  ggpubr::ggarrange(plotlist = figs,
                    nrow = 1)
  
  # Cleaning
  rm(sampsToCheck, dataToCheck, figs)
  
  
  ### Verify gene expression boxplots per sample in muscle, step by step
  
  # On norm data from GEO
  toplot0 <- originalCountData %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    pivot_longer(!ID,
                 names_to = "originalSampleID",
                 values_to = "count") %>%
    left_join(metadata[,c("originalSampleID", "sampleID", "treatment", "timepoint", "tissue")],
              by = "originalSampleID") %>%
    filter(tissue == "muscle")
  tmp <- toplot0 %>%
    dplyr::select(sampleID, treatment, timepoint) %>%
    distinct() %>%
    arrange(treatment, timepoint) %>%
    dplyr::select(sampleID) %>%
    unlist() %>%
    unname()
  toplot0$sampleID <- factor(toplot0$sampleID,
                             levels = tmp)
  rm(tmp)
  gp0 <- ggplot(data = toplot0) +
    geom_boxplot(aes(x = sampleID,
                     y = count,
                     fill = treatment)) +
    theme_bw(base_size = 14) +
    labs(x = "Sample",
         y = "Normalised count",
         title = "Original count data",
         fill = "") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 5,
                                     face = "bold"))
  ggsave(file.path(projectPath,
                   outputPath,
                   mousePath,
                   "QC",
                   "verif_gene_expression_muscle_from_GEO_0.png"),
         plot = gp0,
         width = 20,
         height = 6)
  rm(toplot0, gp0)
  
  # On norm data from GEO limiting to protein coding probes
  toplot1 <- countDataFiltered %>%
    dplyr::select(-SPOT_ID) %>%
    pivot_longer(!ID,
                 names_to = "sampleID",
                 values_to = "count") %>%
    left_join(metadata[,c("originalSampleID", "sampleID", "treatment", "timepoint", "tissue")],
              by = "sampleID") %>%
    filter(tissue == "muscle")
  tmp <- toplot1 %>%
    dplyr::select(sampleID, treatment, timepoint) %>%
    distinct() %>%
    arrange(treatment, timepoint) %>%
    dplyr::select(sampleID) %>%
    unlist() %>%
    unname()
  toplot1$sampleID <- factor(toplot1$sampleID,
                             levels = tmp)
  rm(tmp)
  gp1 <- ggplot(data = toplot1) +
    geom_boxplot(aes(x = sampleID,
                     y = count,
                     fill = treatment)) +
    theme_bw(base_size = 14) +
    labs(x = "Sample",
         y = "Normalised count",
         title = "Original count data - protein coding only",
         fill = "") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 5,
                                     face = "bold"))
  ggsave(file.path(projectPath,
                   outputPath,
                   mousePath,
                   "QC",
                   "verif_gene_expression_muscle_from_GEO_1.png"),
         plot = gp1,
         width = 20,
         height = 6)
  rm(toplot1, gp1)
  
  # On norm data from GEO after avereps
  toplot2 <- dataAllTissues[["muscle"]] %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    pivot_longer(!ID,
                 names_to = "sampleID",
                 values_to = "count") %>%
    left_join(metadata[,c("originalSampleID", "sampleID", "treatment", "timepoint", "tissue")],
              by = "sampleID") %>%
    filter(tissue == "muscle")
  tmp <- toplot2 %>%
    dplyr::select(sampleID, treatment, timepoint) %>%
    distinct() %>%
    arrange(treatment, timepoint) %>%
    dplyr::select(sampleID) %>%
    unlist() %>%
    unname()
  toplot2$sampleID <- factor(toplot2$sampleID,
                             levels = tmp)
  rm(tmp)
  gp2 <- ggplot(data = toplot2) +
    geom_boxplot(aes(x = sampleID,
                     y = count,
                     fill = treatment)) +
    theme_bw(base_size = 14) +
    labs(x = "Sample",
         y = "Normalised count",
         title = "Count data after avereps - protein coding only",
         fill = "") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 5,
                                     face = "bold"))
  ggsave(file.path(projectPath,
                   outputPath,
                   mousePath,
                   "QC",
                   "verif_gene_expression_muscle_from_GEO_2.png"),
         plot = gp2,
         width = 20,
         height = 6)
  rm(toplot2, gp2)
  
  # On norm data from GEO after avereps + filtering of low count genes
  toplot3 <- dataFilteredPerTissue[["muscle"]] %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    pivot_longer(!ID,
                 names_to = "sampleID",
                 values_to = "count") %>%
    left_join(metadata[,c("originalSampleID", "sampleID", "treatment", "timepoint", "tissue")],
              by = "sampleID") %>%
    filter(tissue == "muscle")
  tmp <- toplot3 %>%
    dplyr::select(sampleID, treatment, timepoint) %>%
    distinct() %>%
    arrange(treatment, timepoint) %>%
    dplyr::select(sampleID) %>%
    unlist() %>%
    unname()
  toplot3$sampleID <- factor(toplot3$sampleID,
                             levels = tmp)
  rm(tmp)
  gp3 <- ggplot(data = toplot3) +
    geom_boxplot(aes(x = sampleID,
                     y = count,
                     fill = treatment)) +
    theme_bw(base_size = 14) +
    labs(x = "Sample",
         y = "Normalised count",
         title = "Count data after avereps - filtering low count genes",
         fill = "") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 5,
                                     face = "bold"))
  ggsave(file.path(projectPath,
                   outputPath,
                   mousePath,
                   "QC",
                   "verif_gene_expression_muscle_from_GEO_3.png"),
         plot = gp3,
         width = 20,
         height = 6)
  rm(toplot3, gp3)
}




######################
### QUALITY CHECKS ###
######################

### Retrieve color information for figures
treatmentColors <- colorsPerTreatment$mouse$colors
names(treatmentColors) <- colorsPerTreatment$mouse$treatments


### Boxplots of gene expression per sample

bpPerTissue <- lapply(tissues, function(tiss) {
  # The data to plot
  toplot <- dataFilteredPerTissue[[tiss]] %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    pivot_longer(!ID,
                 names_to = "sampleID",
                 values_to = "count") %>%
    left_join(metadata[,c("originalSampleID", "sampleID", "treatment", "timepoint", "tissue")],
              by = "sampleID")
  # Arrange to have all samples of same group side by side
  tmp <- toplot %>%
    dplyr::select(sampleID, treatment, timepoint) %>%
    distinct() %>%
    arrange(treatment, timepoint) %>%
    dplyr::select(sampleID) %>%
    unlist() %>%
    unname()
  toplot <- toplot %>%
    mutate(sampleID = factor(sampleID,
                             levels = tmp))
  rm(tmp)
  # Make the boxplot
  gp <- ggplot(data = toplot) +
    geom_boxplot(aes(x = sampleID,
                     y = count,
                     fill = treatment)) +
    theme_bw(base_size = 14) +
    scale_fill_manual(values = treatmentColors) +
    labs(x = "Sample",
         y = "Normalised count",
         title = "",
         fill = "") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 5,
                                     face = "bold"))
  # Save the boxplot
  ggsave(file.path(projectPath,
                   outputPath,
                   mousePath,
                   "QC",
                   paste0("gene_expression_boxplot_normalised_",
                          tiss,
                          ".png")),
         plot = gp,
         width = 20,
         height = 6)
  return(gp)
})
rm(bpPerTissue)


### PCA projection colored per variable of interest
pcaVars <- c("treatment", "timepoint")

# Loop on tissues of interest
pcaOutput <- lapply(tissues, function(tiss) {
  # The data to use for PCA
  data <- dataFilteredPerTissue[[tiss]]
  # To run PCA, expected to have data as matrix with row <-> sample and col <-> gene
  pcaRes <- prcomp(t(data), center = T, scale. = F)
  pcaPcVarExplained <- pcaRes$sdev^2/sum(pcaRes$sdev^2)
  tmp <- pcaRes$x %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleID")
  # Data to plot PCA
  toplot <- metadata %>%
    inner_join(tmp, by = "sampleID")
  rm(tmp)
  # Plot PCA per variable of interest
  for (pcaVar in pcaVars)
  {
    gp <- ggplot(data = toplot) +
      theme_bw(base_size = 14) +
      geom_point(aes(x = PC1,
                     y = PC2,
                     colour = get(pcaVar))) +
      stat_ellipse(aes(x = PC1,
                       y = PC2,
                       colour = get(pcaVar))) +
      labs(x = paste0("PC1 (",
                      round(100*pcaPcVarExplained[1],
                            digits = 2), 
                            "%)"),
           y = paste0("PC2 (",
                      round(100*pcaPcVarExplained[2],
                            digits = 2),
                      "%)"),
           color = "",
           title = paste0(Hmisc::capitalize(pcaVar), " - ", tiss),
           subtitle = paste0("(", nrow(data), " genes included)"))
    if (pcaVar == "treatment")
    {
      gp <- gp +
        scale_color_manual(values = treatmentColors)
    }
    ggsave(filename = file.path(projectPath,
                                outputPath,
                                mousePath,
                                "QC",
                                paste0("PCA_",
                                       tiss,
                                       "_with_all_timepoints_colored_by_",
                                       pcaVar,
                                       ".png")),
           plot = gp,
           width = 8, height = 6)
  }
  rm(pcaVar, gp)
  # Return PCA
  return(pcaRes)
})
names(pcaOutput) <- tissues


### Cleaning
rm(list = ls())
