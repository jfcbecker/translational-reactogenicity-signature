###########################################################
### BIOVACSAFE HUMAN DATA - DATA LOAD AND PREPROCESSING ###
###########################################################


### Libraries
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




#########################
### INITIAL DATA LOAD ###
#########################


### Count data - load from GEO (if necessary)
humanGeoDataFile <- file.path(projectPath,
                              dataPath,
                              humanPath,
                              "biovacsafe_human_gse124533_processed_from_GEO.RData")
if (!(file.exists(humanGeoDataFile)))
{
  Sys.setenv(VROOM_CONNECTION_SIZE = 800072) # otherwise error 'The size of the connection buffer (131072) was not large enough'
  
  # Retrieving data from GEO (using R script from website tutorial)
  geoID <- "GSE124533"
  gset <- getGEO(geoID,
                 GSEMatrix = T, # FALSE to get the full data from SOFT files (not only the count table) => takes longer, but more complete (though Agilent flags are not there)
                 getGPL = T, # To get features / probes / genes annotations
                 parseCharacteristics = T) 
  if (length(gset) > 1)
  {
    idx <- grep("GPL21272", attr(gset, "names"))
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
                        humanPath,
                        "biovacsafe_human_gse124533_processed_from_GEO.RData"))
} else
{
  load(file.path(projectPath,
                 dataPath,
                 humanPath,
                 "biovacsafe_human_gse124533_processed_from_GEO.RData"))
}

# Get count data
countData <- exprs(gset)
# => NB: This is background corrected and normalised data (see Weiner et al. 2019)

# Get metadata
metadata <- pData(gset)

# Get features / probes info
features <- fData(gset)

# NB: 1606 samples, 62975 transcripts
rm(gset)


### Laboratory data - read downloaded file from GEO
labData <- data.table::fread(file.path(projectPath,
                                       dataPath,
                                       humanPath,
                                       "GSE124533_CDISCSDTM_CRC305ABC_LB.csv.gz"))


### Adverse events - read downloaded file from GEO
aeData <- data.table::fread(file.path(projectPath,
                                      dataPath,
                                      humanPath,
                                      "GSE124533_CDISCSDTM_CRC305ABC_AE.csv.gz"))


### PTX3 supplementary data (from Giuseppe)
ptx3File <- file.path(projectPath,
                      dataPath,
                      humanPath,
                      "PTX3_Results_updated.xlsx")
ptx3Raw <- read.xlsx(ptx3File,
                     sheet = 1)
rm(ptx3File)


### Adverse events supplementary classification (from Laurent B.)
aeFile <- file.path(projectPath,
                    dataPath,humanPath,
                    "biovacsafe_human_blood_only_distinct_ae_with_counts_and_reviewed_categories.xlsx")
aeGroupsReviewed <- read.xlsx(aeFile,
                              sheet = 1)
rm(aeFile)




######################
### DATA WRANGLING ###
######################

### Formatting metadata

# Focusing on most important information
metadataFull <- metadata
clnms <- colnames(metadata)
charaClnms <- clnms[grepl("characteristics", clnms)]
metadata <- metadata %>%
  dplyr::select(all_of(c("title", "geo_accession", charaClnms))) %>%
  rename(sampleID = geo_accession)

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

# Proper handling of duration (to get it numeric)
metadata$age <- gsub("y", "", metadata$age)

# Additional handling of numeric fields
numClnms <- c("day", "age")
metadata <- metadata %>%
  mutate(across(.cols = all_of(numClnms),
                .fns = as.numeric))
rm(numClnms)

# Definition of additional / corrected variables
tps <- paste0(sort(unique(as.integer(metadata$day*24))), "h")
metadata <- metadata %>%
  # Proper naming of tissue
  mutate(tissue = "blood") %>%
  mutate(tissueAcronym = "BL") %>%
  # Proper naming of treatments
  mutate(treatment = case_when(
    treatment %in% c("AGRIPPAL", "VARILRIX", "STAMARIL") ~ Hmisc::capitalize(tolower(treatment)),
    treatment == "ENGERIXB1" ~ "Engerix B",
    treatment == "ENGERIXB3" ~ "Engerix B3",
    treatment == "FLUADC" ~ "Fluad",
    treatment == "PLACEBOAB1C" ~ "Saline",
    treatment == "PLACEBOB3" ~ "Saline B3"
  )) %>%
  mutate(treatmentAcronym = case_when(
    treatment == "Engerix B" ~ "ENG",
    treatment == "Engerix B3" ~ "ENGB3",
    treatment == "Varilrix" ~ "VZV",
    treatment == "Stamaril" ~ "YFV",
    treatment == "Saline" ~ "SAL",
    treatment == "Agrippal" ~ "TriFLU",
    treatment == "Fluad" ~ "TriFLU.MF59",
    treatment == "Saline B3" ~ "SALB3",
  )) %>%
  # Reshaping timepoint
  mutate(timepoint = factor(paste0(day*24, "h"),
                            levels = tps),
         timepointDay = paste0("D", day)) %>%
  # Define new sample ID based on tissue, treatment, timepoint and and animal ID
  mutate(originalSampleID = sampleID) %>%
  mutate(sampleID = paste0(treatmentAcronym,
                           "_",
                           tissueAcronym,
                           "_",
                           timepoint,
                           "_",
                           participant))
rm(tps)


### Formatting count data

# Arrange samples by tissue type, treatment and timepoint in count data and metadata
metadata <- metadata %>%
  arrange(tissue, treatment, timepoint)
countData <- countData[, metadata$originalSampleID]

# Arrange features in same order in count data and feature list
features <- features %>%
  mutate(ID = as.character(ID)) %>%
  arrange(ID)
countData <- countData[features$ID,]

# Adding features ID to data frame
originalCountData <- countData
countData <- countData %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID")
# Correct colname = sample ID with metadata included
colnames(countData) <- c("ID", metadata$sampleID)


### Formatting lab data

labData <- labData %>%
  # Rename some columns to simplify
  rename("labTest" = "LBTESTCD",
         "day" = "LBDY",
         "timeOfDay" = "LBTPT",
         "labTestValue" = "LBSTRESN",
         "labTestUnit" = "LBSTRESU")


### Formatting PTX3 data as other lab data

# Extract info on patient, study and visits for each patient
patientVisitDetails <- labData %>%
  # Extract details on visits
  dplyr::select(STUDYID,
                USUBJID,
                VISITNUM,
                VISIT,
                LBDTC,
                day) %>%
  distinct() %>%
  # Extract "simplified" patient ID (this is the one that is used in PTX3 data)
  mutate(simplePatientID = str_split(USUBJID, pattern = "-", simplify = T)[,2])
patientIDs <- patientVisitDetails %>%
  dplyr::select(STUDYID,
                USUBJID,
                simplePatientID) %>%
  distinct()

# Reformat PTX3 data as appropriate
clnms <- colnames(ptx3Raw)
clnmsToPivot <- clnms[!(clnms %in% c("Treatment", "Subject.ID"))]
ptx3CleanData <- ptx3Raw %>%
  # Set patient ID as character
  mutate(Subject.ID = as.character(Subject.ID)) %>%
  # Add information of lab testing
  mutate(labTest = "PTX3",
         labTestUnit = "ng/mL") %>%
  # Set time of sampling in one column, and all measurements in one single column (pivot_longer)
  pivot_longer(cols = all_of(clnmsToPivot),
               names_to = "timepoint",
               values_to = "labTestValue") %>%
  # Rename some columns
  rename("treatment" = "Treatment",
         "simplePatientID" = "Subject.ID") %>%
  # Add full patient ID and study ID (from lab data)
  left_join(patientIDs,
            by = "simplePatientID") %>%
  # Proper formatting of timepoint / visit / day information
  mutate(timepointChar = timepoint,
         timepoint = str_remove(timepointChar,
                                pattern = "D")) %>%
  mutate(day = if_else(str_detect(timepoint,
                                  pattern = ":"),
                       as.numeric(str_split(timepoint,
                                            pattern = "_",
                                            simplify = T)[,1]),
                       as.numeric(timepoint)),
         timeOfDay = if_else(str_detect(timepoint,
                                        pattern = ":"),
                             str_split(timepoint,
                                       pattern = "_",
                                       simplify = T)[,2],
                             NA)) %>%
  # Proper formatting of treatment column
  mutate(treatment = if_else(treatment == "PLACEBO",
                             "Saline",
                             Hmisc::capitalize(tolower(treatment)))) %>%
  # Define participant as defined in lab data
  mutate(participant = str_remove(USUBJID, "-")) %>%
  # Rearrange columns
  dplyr::select(-timepoint) %>%
  dplyr::select(STUDYID,
                USUBJID,
                simplePatientID,
                treatment,
                labTest,
                labTestValue,
                labTestUnit,
                timepointChar,
                day,
                timeOfDay)
rm(clnms, clnmsToPivot)


### Formatting adverse events data

# AE data
aeData <- aeData %>%
  # Rename some columns to simplify
  rename("severity" = "AESEV",
         "dayAfterVacc" = "AESTDY",
         "duration" = "AEDUR",
         "relatedToVacc" = "AEREL",
         "aeTerm" = "AETERM",
         "aeCategory" = "AEHLGT",
         "aeCategory2" = "AESOC") %>%
  # Convert severity (mild, moderate, severe) to numeric value (1, 2, 3)
  mutate(severityNum = case_when(
    severity == "MILD" ~ 1,
    severity == "MODERATE" ~ 2,
    severity == "SEVERE" ~ 3
  ))

# AE categories
aeGroupsReviewed <- aeGroupsReviewed %>%
  # Select relevant columns
  dplyr::select(all_of(c("aeTerm",
                         "aeCategory",
                         "aeCategory2",
                         "include",
                         "reviewedAeCategory"))) %>%
  # Convert 'include' to boolean
  mutate(include = case_when(
    include == "T" ~ TRUE,
    include == "F" ~ FALSE
  ))

# Combine both
aeData <- aeData %>%
  left_join(aeGroupsReviewed,
            by = c("aeTerm", "aeCategory", "aeCategory2"))
rm(aeGroupsReviewed)
  
  
  
  
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
  # Remove probes with no gene symbol or ENTREZ gene ID or ENSEMBL ID
  filter(!is.na(ENTREZ_GENE_ID)
         | !is.na(GENE_SYMBOL)
         | !is.na(ENSEMBL_ID))

# Features to keep, without duplicates
featsToKeepNoDuplik <- featsToKeep %>%
  dplyr::select(SPOT_ID, ENTREZ_GENE_ID, GENE_SYMBOL, GENE_NAME, ENSEMBL_ID, REFSEQ) %>%
  distinct()

# Features to keep, based on biotype from Biomart
biomartMatch <- read.table(file.path(projectPath,
                                     dataPath,
                                     humanPath,
                                     "biomart_match_human_features_agilent_8x60K_v1_v2.txt"),
                           sep = ",",
                           header = T,
                           stringsAsFactors = F)
biomartMatch <- biomartMatch %>%
  rename(geneID = Gene.stable.ID,
         transcriptID = Transcript.stable.ID,
         geneType = Gene.type,
         transcriptType = Transcript.type,
         geneSymbol = HGNC.symbol,
         entrezGeneID = NCBI.gene..formerly.Entrezgene..ID,
         spotIDv1 = AGILENT.SurePrint.G3.GE.8x60k.probe,
         spotIDv2 = AGILENT.SurePrint.G3.GE.8x60k.v2.probe)

# Try hybrid match, using gene symbol, then ENSEMBL, then Entrez gene ID
tomatch <- featsToKeepNoDuplik
# Match first on gene symbol
tmp <- tomatch %>%
  filter(!is.na(GENE_SYMBOL)) %>%
  rename(geneSymbol = GENE_SYMBOL)
tmp2 <- biomartMatch %>%
  filter(!is.na(geneSymbol)) %>%
  filter(transcriptType == "protein_coding") %>%
  dplyr::select(geneSymbol, transcriptType) %>%
  distinct()
matchSymbol <- tmp %>%
  left_join(tmp2, by = c("geneSymbol")) %>%
  rename(GENE_SYMBOL = geneSymbol)
nrow(matchSymbol) == nrow(tmp) # Check if not introducing duplicates
rm(tmp, tmp2)
# Then match remaining features on ensembl
tmp <- tomatch %>%
  filter(!(SPOT_ID %in% matchSymbol$SPOT_ID)) %>%
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
# Finally match remaining features on Entrez gene ID
# => actually not needed, all matches found with gene symbol and ENSEMBL
# Merge the two versions of match (no unmatched features)
fullMatch <- matchSymbol %>%
  bind_rows(matchEnsembl)
nrow(fullMatch) == nrow(tomatch) # Check if good number of transcripts
rm(matchSymbol, matchEnsembl, tomatch)
# Exclude non proteing coding transcripts
proteinCodingFeats <- fullMatch %>%
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
tmp <- countDataFiltered %>%
  dplyr::select(-ID) %>%
  dplyr::select(-SPOT_ID)
ids <- countDataFiltered$SPOT_ID
countDataFilteredAveraged <- limma::avereps(tmp,
                                            ID = ids)
# => output = matrix with rows = unique probes (rownames = SPOT_ID already set), cols = samples
rm(tmp, ids)




#########################################################
### FILTERING OUT INFORMATION NOT NEEDED FOR ANALYSIS ###
#########################################################

### Metadata
tps <- metadata %>%
  filter(day <= 7) %>%
  dplyr::select(day, timepoint) %>%
  distinct() %>%
  arrange(day) %>%
  dplyr::select(timepoint) %>%
  unlist() %>%
  unname() %>%
  as.character()
metadata <- metadata %>%
  # Exclude samples of HBV3 and corresponding placebo samples
  filter(!(treatment %in% c("Engerix B3", "Saline B3"))) %>%
  # Exclude samples after D7
  filter(day <= 7) %>%
  # Redefine properly timepoints levels
  mutate(timepoint = factor(timepoint,
                            levels = tps))
rm(tps)


### Count data

# Exclude samples of HBV3 and corresponding placebo samples, as well as samples before D0 and after D7
countDataFilteredAveraged <- countDataFilteredAveraged[,metadata$sampleID]
dataAllTissues <- list()
dataAllTissues[["blood"]] <- countDataFilteredAveraged


### Lab data

tmp <- metadata %>%
  dplyr::select(participant,
                treatment) %>%
  distinct()
labData <- labData %>%
  # Select columns of interest
  dplyr::select(all_of(c("STUDYID",
                         "participant",
                         "USUBJID",
                         "labTest",
                         "day",
                         "timeOfDay",
                         "labTestValue",
                         "labTestUnit"))) %>%
  # Limit to CRP
  filter(labTest == "CRP") %>%
  # Limit to valid timepoint before D7
  filter(!is.na(day)
         & day <= 7) %>%
  # Add treatment information
  left_join(tmp, by = "participant") %>%
  # Exclude HBV3 and corresponding placebo samples
  filter(!is.na(treatment))
rm(tmp)

# Do the same cleaning with PTX3
ptx3CleanData <- ptx3CleanData %>%
  # Limit to valid timepoint before and D7
  filter(!is.na(day)
         & day <= 7)


### AE data with AE categories
aeData <- aeData %>%
  # Select columns of interest
  dplyr::select(all_of(c("STUDYID",
                         "participant",
                         "USUBJID",
                         "relatedToVacc",
                         "severity",
                         "severityNum",
                         "dayAfterVacc",
                         "duration",
                         "aeTerm",
                         "aeCategory",
                         "aeCategory2",
                         "reviewedAeCategory",
                         "include"))) %>%
  # Limit to adverse events starting between D0 and D7 after immunisation
  filter(!is.na(dayAfterVacc)
         & dayAfterVacc >= 0
         & dayAfterVacc <= 7) %>%
  # Limit to adverse events related to immunisation
  filter(!str_detect(relatedToVacc, "NOT RELATED")) %>%
  # Exclude local complications from indwelling cannula
  filter(!str_detect(aeTerm, "CANNULA SITE")) %>%
  # Exclude laboratory AEs
  filter(!str_detect(aeCategory, "INVESTIGATIONS")
         | str_detect(aeCategory, "CARDIAC")) %>%
  filter(!str_detect(aeCategory, "PROTEIN AND CHEMISTRY ANALYSES NEC")) %>%
  # Calculate reactoscore per AE
  mutate(reactoscore = severityNum*duration)
# If need be, exclude other AEs based on review by Laurent B.
aeDataFull <- aeData
aeData <- aeData %>%
    filter(include == T)




###############################
### REACTOSCORE CALCULATION ###
###############################

reactoscoreCategory <- c("all", "systemic", "local")

### Calculate individual reactoscore (overall, systemic or local)
patients <- unique(metadata$participant)
allIndividualReactoscores <- lapply(list(aeDataFull, aeData), function(x) {
  curIndividualReactoscores <- lapply(reactoscoreCategory, function(aeCat) {
    # Loop on patients
    doer <- lapply(patients, function(pat) {
      # Metadata of interest
      curMeta <- metadata %>%
        filter(participant == pat) %>%
        dplyr::select(participant, treatment) %>%
        distinct()
      # Limit to AEs of interest
      if (aeCat == "all")
      {
        touse <- x
      } else
      {
        touse <- x %>%
          filter(reviewedAeCategory == aeCat)
      }
      touse <- touse %>%
        filter(participant == pat)
      # Calculate individual reactoscore
      if (nrow(touse) > 0)
      {
        patReactoscore = sum(touse$reactoscore, na.rm = T)
      }
      else
      {
        patReactoscore = 0
      }
      res <- curMeta %>%
        mutate(aeCategory = aeCat,
               reactoscore = patReactoscore)
      return(res)
    })
    res <- do.call("bind_rows", doer)
    return(res)
  })
  curIndividualReactoscores <- do.call("bind_rows", curIndividualReactoscores)
  return(curIndividualReactoscores)
})
names(allIndividualReactoscores) <- c("full", "restricted")
individualReactoscores <- do.call("bind_rows", allIndividualReactoscores[["restricted"]]) %>%
  arrange(treatment,
          participant,
          aeCategory)
individualReactoscoresFull <- do.call("bind_rows", allIndividualReactoscores[["full"]]) %>%
  arrange(treatment,
          participant,
          aeCategory)
rm(patients)


### Calculate aggregated reactoscore metrics (overall, systemic or local)

# Mean individual reactoscore per vaccine and reactosum per vaccine
allAggregReactoscorePerVaccine <- lapply(c("full", "restricted"), function(subj) {
  if (subj == "full")
  {
    x <- aeDataFull
    y <- individualReactoscoresFull
  } else if (subj == "restricted")
  {
    x <- aeData
    y <- individualReactoscores
  }
  curAggregReactoscorePerVaccine <- lapply(reactoscoreCategory, function(aeCat) {
    # Loop on treatments
    doer <- lapply(unique(metadata$treatment), function(vacc) {
      # Patients for considered treatment
      vaccPatients <- metadata %>%
        filter(treatment == vacc) %>%
        dplyr::select(participant) %>%
        distinct() %>%
        unlist() %>%
        unname()
      # AEs for considered patients
      if (aeCat == "all")
      {
        vaccAEs <- x
        vaccIndivAEs <- y
      } else
      {
        vaccAEs <- x %>%
          filter(reviewedAeCategory == aeCat)
        vaccIndivAEs <- y %>%
          filter(aeCategory == aeCat)
      }
      vaccAEs <- vaccAEs %>%
        filter(participant %in% vaccPatients)
      vaccIndivAEs <- vaccIndivAEs %>%
        filter(participant %in% vaccPatients)
      # Calculate mean individual reactoscore
      meanIndividualReactoscore <- mean(vaccIndivAEs$reactoscore, na.rm = T)
      # Calculate reactosum
      reactosum <- sum(vaccAEs$reactoscore, na.rm = T)
      # Result to keep
      res <- data.frame(treatment = vacc,
                        reactosum = reactosum,
                        meanPatientReactoscore = meanIndividualReactoscore,
                        numPatients = length(vaccPatients),
                        numPatientsWithAE = nrow(vaccIndivAEs),
                        numAE = nrow(vaccAEs),
                        stringsAsFactors = F)
      return(res)
    })
    doer <- do.call("bind_rows", doer) %>%
      mutate(aeCategory = aeCat)
    return(doer)
  })
  curAggregReactoscorePerVaccine <- do.call("bind_rows", curAggregReactoscorePerVaccine)
  return(curAggregReactoscorePerVaccine)
})
names(allAggregReactoscorePerVaccine) <- c("full", "restricted")
aggregReactoscorePerVaccine <- do.call("bind_rows", allAggregReactoscorePerVaccine[["restricted"]])
aggregReactoscorePerVaccineFull <- do.call("bind_rows", allAggregReactoscorePerVaccine[["full"]])




###################
### SAVING DATA ###
###################

# All in one RData file
save(countDataFilteredAveraged,
     metadata,
     labData,
     ptx3CleanData,
     aeData,
     featsToKeep,
     featsToKeepNoDuplik,
     featsToKeepNoDuplikProtCoding,
     summaryNbFeats,
     file = file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_data_for_analysis.RData"))

# Counts (filtered -> no control probes left)
write.table(countDataFilteredAveraged,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_blood_count_data_for_analysis.txt"),
            sep = "\t",
            quote = F,
            row.names = T)

# Metadata
write.table(metadata,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_metadata_for_analysis.txt"),
            sep = "\t",
            quote = F,
            row.names = F)

# Laboratory data
write.table(labData,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_crp_data_for_analysis.txt"),
            sep = "\t",
            quote = F,
            row.names = F)

# Formatted PTX3 data
write.table(ptx3CleanData,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_ptx3_data_for_analysis.txt"),
            sep = "\t",
            quote = F,
            row.names = F)

# Adverse events data
write.table(aeData,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_ae_data_with_categories_for_analysis.txt"),
            sep = "\t",
            quote = F,
            row.names = F)
write.table(aeDataFull,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_ae_data_including_subjective_with_categories_for_analysis.txt"),
            sep = "\t",
            quote = F,
            row.names = F)


# Features / mapping info (filtered -> no control probes left)
write.xlsx(featsToKeepNoDuplik,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_features_without_duplicates_for_analysis.xlsx"),
            quote = F,
            rowNames = F)
write.xlsx(featsToKeepNoDuplikProtCoding,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_features_protein_coding_without_duplicates_for_analysis.xlsx"),
            quote = F,
            rowNames = F)

# Summary of number of features across filtering
write.table(summaryNbFeats,
            file.path(projectPath,
                      dataPath,
                      humanPath,
                      "biovacsafe_human_features_filtering_summary_table.csv"),
            sep = ";",
            quote = F,
            row.names = F)

# Reactoscore values
write.csv(individualReactoscores,
          file.path(projectPath,
                    dataPath,
                    humanPath,
                    "biovacsafe_human_ae_patients_reactoscore_per_category.csv"),
          quote = F,
          row.names = F)
write.csv(individualReactoscoresFull,
          file.path(projectPath,
                    dataPath,
                    humanPath,
                    "biovacsafe_human_ae_patients_reactoscore_including_subjective_per_category.csv"),
          quote = F,
          row.names = F)
write.csv(aggregReactoscorePerVaccine,
          file.path(projectPath,
                    dataPath,
                    humanPath,
                    "biovacsafe_human_ae_vaccine_reactoscore_per_category.csv"),
          quote = F,
          row.names = F)
write.csv(aggregReactoscorePerVaccineFull,
          file.path(projectPath,
                    dataPath,
                    humanPath,
                    "biovacsafe_human_ae_vaccine_reactoscore_including_subjective_per_category.csv"),
          quote = F,
          row.names = F)




#################################
### LOW COUNT GENES FILTERING ###
#################################

### Create appropriate directory for QCs and corrections
dir.create(file.path(projectPath,
                     outputPath,
                     humanPath),
           showWarnings = F)

dir.create(file.path(projectPath,
                     outputPath,
                     humanPath,
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
                humanPath,
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
  abline(v = lowCountThreshold[["human"]][[tiss]], col = "red")
  dev.off()
  # Cleaning
  rm(tmpData, tmpMedian)
}
rm(tiss)


### Perform filtering for each tissue type
# Criterion 1: keep only genes with > X counts in at least K samples, with K = nb samples in smallest tissue-TP-treatment group
# => see https://f1000research.com/articles/5-1384"
# Criterion 2: keep only genes with sum of counts across all samples >= X*K with X the low count threshold and K the nb of samples in smallest tissue-TP-treatment group
# => HERE WE USE CRITERION 2
genesToKeepPerTissue <- lapply(tissues, function(tiss) {
  # Number of samples per treatment - timepoint group
  tmp <- metadata %>%
    filter(tissue == tiss) %>%
    group_by(treatment, timepoint) %>%
    dplyr::summarise(n = n())
  K <- floor(min(tmp$n))
  # List of genes passing the filtering criterion
  genesToKeep1 <- rownames(dataAllTissues[[tiss]])[rowSums(dataAllTissues[[tiss]] > lowCountThreshold[["human"]][[tiss]]) >= K]
  genesToKeep2 <- rownames(dataAllTissues[[tiss]])[rowSums(dataAllTissues[[tiss]]) > lowCountThreshold[["human"]][[tiss]]*K]
  res <- list(lowCountThreshold = lowCountThreshold[["human"]][[tiss]],
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
                      humanPath,
                      "QC",
                      paste0("biovacsafe_human_count_data_filtering_low_count_genes_",
                             tiss,
                             ".csv")),
            quote = F,
            row.names = T)
}
rm(tiss)




######################
### QUALITY CHECKS ###
######################

### Retrieve color information for figures
treatmentColors <- colorsPerTreatment$human$colors
names(treatmentColors) <- colorsPerTreatment$human$treatments


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
                   humanPath,
                   "QC",
                   paste0("gene_expression_boxplot_normalised_",
                          tiss,
                          ".png")),
         plot = gp,
         width = 30,
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
                                humanPath,
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
