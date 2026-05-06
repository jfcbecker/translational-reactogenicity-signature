###########################################################
### SIBYLLE PAPER - ANALYSIS OF REACTOGENICITY READOUTS ###
###########################################################


### Libraries
library(tidyverse)
library(GEOquery)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(patchwork)


### Parameters and paths are defined in the YAML file
params <- yaml::read_yaml(file = file.path(this.path::here(),
                                           "sibylle_mouse_human_data_preprocessing_parameters.yaml"))
suppressWarnings(params %>% attach) # To enable accessing the parameters by their name without using params$...
tissues <- list(mouse = tissuesOfInterest$mouse,
                human = tissuesOfInterest$human)
# Set appropriate working directory
setwd(projectPath)


### Choose the model to focus on
modelToUse  <- "MoMuOnHuBl_2classes_lambdaOptim"
alternativeModelToUse <- "MoMuOnHuBl_3classes_lambdaOptim"


### Precise if correlation figures should be done or not (takes some time, if already done no need to redo)
runCorrelFigures <- F


### Precise which version of the reactoscore should be used
filterSubjectiveAEs <- T


### Create folders for correlations between predictions and clinical readouts
dir.create(file.path(projectPath,
                     outputPath,
                     "predictionCorrelations"),
           showWarnings = F)
dir.create(file.path(projectPath,
                     outputPath,
                     "predictionCorrelations",
                     modelToUse),
           showWarnings = F)
dir.create(file.path(projectPath,
                     outputPath,
                     "predictionCorrelations",
                     modelToUse,
                     "per_vaccine"),
           showWarnings = F)
dir.create(file.path(projectPath,
                     outputPath,
                     "predictionCorrelations",
                     modelToUse,
                     "per_individual"),
           showWarnings = F)
dir.create(file.path(projectPath,
                     outputPath,
                     "predictionCorrelations",
                     "reactogenicityReadoutsVisualisation"),
           showWarnings = F)





#######################
### LOAD HUMAN DATA ###
#######################

### Metadata
humanMetadata <- read.table(file.path(projectPath,
                                      dataPath,
                                      humanPath,
                                      "biovacsafe_human_metadata_for_analysis.txt"),
                            sep = "\t",
                            header = T,
                            stringsAsFactors = F) %>%
  # Replace placebo by saline if need be
  mutate(treatment = if_else(treatment == "Placebo", "Saline", treatment))


### Laboratory data

# CRP
crpData <- read.table(file.path(projectPath,
                                dataPath,
                                humanPath,
                                "biovacsafe_human_crp_data_for_analysis.txt"),
                      sep = "\t",
                      header = T,
                      stringsAsFactors = F) %>%
  # Keep only measurements with timeOfDay != "UNSCHEDULED" (these are empty measurements)
  filter(timeOfDay != "UNSCHEDULED 1") %>%
  # Replace placebo by saline if need be
  mutate(treatment = if_else(treatment == "Placebo", "Saline", treatment))

# PTX3
ptx3Data <- read.table(file.path(projectPath,
                                dataPath,
                                humanPath,
                                "biovacsafe_human_ptx3_data_for_analysis.txt"),
                      sep = "\t",
                      header = T,
                      stringsAsFactors = F) %>%
  # At D0 and D1, keep only measurements at 08:00 (for other days, only one measurement and not timeOfDay)
  filter(is.na(timeOfDay)
         | timeOfDay == "08:00") %>%
  # Define variable 'participant' as in other datasats
  mutate(participant = str_replace_all(USUBJID, "-", "")) %>%
  # Replace placebo by saline if need be
  mutate(treatment = if_else(treatment == "Placebo", "Saline", treatment))


### Reactoscore values

# Per patient
individualReactoscores <- read.csv(file.path(projectPath,
                                               dataPath,
                                               humanPath,
                                               paste0("biovacsafe_human_ae_patients_reactoscore_",
                                                      ifelse(filterSubjectiveAEs == F, "including_subjective_", ""),
                                                      "per_category.csv")),
                                   header = T,
                                   stringsAsFactors = F) %>%
  # Replace placebo by saline if need be
  mutate(treatment = if_else(treatment == "Placebo", "Saline", treatment))

# Aggregated per vaccine
aggregReactoscorePerVaccine <- read.csv(file.path(projectPath,
                                                  dataPath,
                                                  humanPath,
                                                  paste0("biovacsafe_human_ae_vaccine_reactoscore_",
                                                         ifelse(filterSubjectiveAEs == F, "including_subjective_", ""),
                                                         "per_category.csv")),
                                        header = T,
                                        stringsAsFactors = F) %>%
  # Replace placebo by saline if need be
  mutate(treatment = if_else(treatment == "Placebo", "Saline", treatment))


### Model predictions (latent variables)

# Load the latent variables from the model of interest
tmp <- lapply(c(modelToUse, alternativeModelToUse), function(modl) {
  res <- read.table(file.path(projectPath,
                              dataPath,
                              "ModelEvaluation",
                              paste0("LatentVariable_",
                                     modl,
                                     ".txt")),
                    sep = " ",
                    fill = T,
                    stringsAsFactors = F)[,(1:8)] %>%
    janitor::row_to_names(row_number = 1) %>%
    # Define treatment as in other datasets
    mutate(treatment = case_when(
      vaccine == "Engerix" ~ "Engerix B",
      vaccine == "Varilix" ~ "Varilrix",
      vaccine == "Placebo" ~ "Saline",
      TRUE ~ vaccine
    )) %>%
    mutate(across(all_of(c("weighted_f1",
                           "lambda",
                           "latent")),
                  .fn = as.numeric)) %>%
    mutate(sample = str_replace_all(sample,
                                    "Varilix",
                                    "Varilrix")) %>%
    mutate(sample = str_replace_all(sample,
                                    "Engerix",
                                    "Engerix B"))
  return(res)
})
names(tmp) <- c(modelToUse, alternativeModelToUse)
reactoPredictions <- tmp[[modelToUse]]
alternativeReactoPredictions <- tmp[[alternativeModelToUse]]
rm(tmp)


  

##########################################################
### CALCULATE AGGREGATED READOUTS AT THE VACCINE LEVEL ###
##########################################################

### Mean CRP and PTX3 levels per vaccine and per timepoint

# CRP
aggregCRPPerVaccineAndTP <- crpData %>%
  group_by(treatment, day) %>%
  summarise(meanLabTestValue = mean(labTestValue, na.rm = T)) %>%
  ungroup()

# PTX3
aggregPTX3PerVaccineAndTP <- ptx3Data %>%
  group_by(treatment, day) %>%
  summarise(meanLabTestValue = mean(labTestValue, na.rm = T)) %>%
  ungroup()


### Mean predictions per vaccine and per timepoints (train + test)
aggregReactoPredictions <- reactoPredictions %>%
  group_by(treatment, tpTrain, tpTest) %>%
  summarise(meanLatentVariable = mean(latent, na.rm = T)) %>%
  ungroup()




######################################################
### COMBINE REACTOGENICITY READOUTS FOR FUTURE USE ###
######################################################

### Per patient

dataForCorrel <- humanMetadata %>%
  # Key columns of metadata
  dplyr::select(participant,
                treatment,
                day,
                sampleID) %>%
  # Define timepoint in hours rather than in days
  mutate(timepoint = paste0(24*day, "h"))

# Join with CRP
tomerge <- crpData %>%
  dplyr::select(participant, day, labTestValue) %>%
  mutate(timepoint = paste0(24*day, "h")) %>%
  pivot_wider(id_cols = "participant",
              names_from = "timepoint",
              values_from = "labTestValue",
              names_prefix = "CRP_")
dataForCorrel <- dataForCorrel %>%
  left_join(tomerge,
            by = "participant")
rm(tomerge)

# Join with PTX3
tomerge <- ptx3Data %>%
  dplyr::select(participant, day, labTestValue) %>%
  mutate(timepoint = paste0(24*day, "h")) %>%
  pivot_wider(id_cols = "participant",
              names_from = "timepoint",
              values_from = "labTestValue",
              names_prefix = "PTX3_")
dataForCorrel <- dataForCorrel %>%
  left_join(tomerge,
            by = "participant")
rm(tomerge)

# Join with reactoscores
aeCats <- unique(individualReactoscores$aeCategory)
tomerge <- NULL
for (aeC in aeCats)
{
  x <- individualReactoscores %>%
    filter(aeCategory == aeC) %>%
    dplyr::select(participant, reactoscore) %>%
    rename(!!paste0("reactoscore_", aeC) := "reactoscore")
  if (is.null(tomerge))
  {
    tomerge <- x
  } else
  {
    tomerge <- tomerge %>%
      left_join(x,
                by = "participant")
  }
  rm(x)
}
rm(aeC, aeCats)
dataForCorrel <- dataForCorrel %>%
  left_join(tomerge,
            by = "participant")
rm(tomerge)

# Join with reactogenicity predictions
tmp <- humanMetadata %>%
  dplyr::select(treatment, treatmentAcronym) %>%
  distinct()
tomerge <- reactoPredictions %>%
  dplyr::select(-name) %>%
  left_join(tmp, by = "treatment") %>%
  mutate(sampleID = str_replace(sample,
                                treatment,
                                treatmentAcronym))
rm(tmp)
dataForCorrel <- tomerge %>%
  left_join(dataForCorrel,
            by = c("treatment", "sampleID"))
rm(tomerge)

# Calculate log2 of clinical readouts for easier use in correlation
clnms <- colnames(dataForCorrel)
clnms <- clnms[!(clnms %in% c("treatment", "tpTrain", "tpTest",
                              "sample", "vaccine", "treatmentAcronym",
                              "sampleID", "participant", "day",
                              "timepoint"))]
for (cl in clnms) 
{
  if (grepl("reactoscore", cl)) # Case of reactoscore, add 1 before log2 to avoid log2(0)
  {
    dataForCorrel <- dataForCorrel %>%
      mutate(across(all_of(cl),
                    .fn = ~ log2(.x +1),
                    .names = "log2_{.col}_plus1"))
  } else
  {
    dataForCorrel <- dataForCorrel %>%
      mutate(across(all_of(cl),
                    .fn = log2,
                    .names = "log2_{.col}"))
  }
}
rm(cl)


### Per vaccine

# Start from reactogenicity predictions
aggregDataForCorrel <- aggregReactoPredictions

# Join with CRP
tomerge <- aggregCRPPerVaccineAndTP %>%
  mutate(timepoint = paste0(24*day, "h")) %>%
  pivot_wider(id_cols = "treatment",
              names_from = "timepoint",
              values_from = "meanLabTestValue",
              names_prefix = "mean_CRP_")
aggregDataForCorrel <- aggregDataForCorrel %>%
  left_join(tomerge,
            by = "treatment")
rm(tomerge)

# Join with PTX3
tomerge <- aggregPTX3PerVaccineAndTP %>%
  mutate(timepoint = paste0(24*day, "h")) %>%
  pivot_wider(id_cols = "treatment",
              names_from = "timepoint",
              values_from = "meanLabTestValue",
              names_prefix = "mean_PTX3_")
aggregDataForCorrel <- aggregDataForCorrel %>%
  left_join(tomerge,
            by = "treatment")
rm(tomerge)

# Join with reactoscore
aeCats <- unique(aggregReactoscorePerVaccine$aeCategory)
tomerge <- NULL
for (aeC in aeCats)
{
  x <- aggregReactoscorePerVaccine %>%
    filter(aeCategory == aeC) %>%
    dplyr::select(treatment, reactosum, meanPatientReactoscore) %>%
    rename(!!paste0("meanPatientReactoscore_", aeC) := "meanPatientReactoscore",
           !!paste0("reactosum_", aeC) := "reactosum")
  if (is.null(tomerge))
  {
    tomerge <- x
  } else
  {
    tomerge <- tomerge %>%
      left_join(x,
                by = "treatment")
  }
  rm(x)
}
rm(aeC, aeCats)
aggregDataForCorrel <- aggregDataForCorrel %>%
  left_join(tomerge,
            by = "treatment")
rm(tomerge)

# Calculate log2 of clinical readouts for easier use in correlation
clnms <- colnames(aggregDataForCorrel)
clnms <- clnms[!(clnms %in% c("treatment", "tpTrain", "tpTest"))]
for (cl in clnms) 
{
  aggregDataForCorrel <- aggregDataForCorrel %>%
    mutate(across(all_of(cl),
                  .fn = log2,
                  .names = "log2_{.col}"))
}
rm(cl)


  
  
#################################################
### DISTRIBUTIONS OF REACTOGENICITY READOUTS ###
#################################################

### Parameters
mycols <- colorsPerTreatment$human$colors
names(mycols) <- colorsPerTreatment$human$treatments
names(mycols)[names(mycols) == "Placebo"] <- "Saline"


### Boxplots of reactoscores
toplot <- individualReactoscores %>%
  mutate(treatment = relevel(factor(treatment),
                             ref = "Saline"))
figName <- file.path(projectPath,
                     outputPath,
                     "predictionCorrelations",
                     "reactogenicityReadoutsVisualisation",
                     paste0("distribution_of_individual_reactoscores_",
                            ifelse(filterSubjectiveAEs == F, "including_subjective_", ""),
                            "per_vaccine.png"))
gp <- ggplot(data = toplot) +
  theme_bw(base_size = 16) +
  geom_boxplot(aes(x = treatment,
                   y = log2(reactoscore + 1),
                   fill = treatment),
               colour = "black") +
  scale_fill_manual(values = mycols) +
  labs(x = "",
       y = "Log2(reactoscore + 1)",
       fill = "") +
  guides(fill = "none")
ggsave(filename = figName,
       plot = gp,
       width = 8,
       height = 6)
rm(toplot, figName, gp)


### Individual CRP levels per vaccine and per timepoint
tps <- unique(crpData$day)*24
tps <- tps[tps >= 0]
tps <- tps[order(tps)]
tps <- paste0(tps, "h")
toplot <- crpData %>%
  mutate(treatment = relevel(factor(treatment),
                             ref = "Saline")) %>%
  filter(day >= 0) %>%
  mutate(timepoint = factor(paste0(as.character(day*24),
                                   "h"),
                            levels = tps))
figName <- file.path(projectPath,
                     outputPath,
                     "predictionCorrelations",
                     "reactogenicityReadoutsVisualisation",
                     "distribution_of_CRP_per_vaccine_per_timepoint.png")
gp <- ggplot(data = toplot) +
  theme_bw(base_size = 16) +
  geom_boxplot(aes(x = treatment,
                   y = log2(labTestValue),
                   fill = treatment),
               colour = "black") +
  scale_fill_manual(values = mycols) +
  labs(x = "",
       y = "Log2(CRP)",
       fill = "") +
  facet_wrap(~ timepoint) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))
ggsave(filename = figName,
       plot = gp,
       width = 10,
       height = 8)
rm(toplot, figName, gp, tps)


### Individual PTX3 levels per vaccine and per timepoint
tps <- unique(crpData$day)*24
tps <- tps[tps >= 0]
tps <- tps[order(tps)]
tps <- paste0(tps, "h")
toplot <- ptx3Data %>%
  mutate(treatment = relevel(factor(treatment),
                             ref = "Saline")) %>%
  filter(day >= 0) %>%
  mutate(timepoint = factor(paste0(as.character(day*24),
                                   "h"),
                            levels = tps))
figName <- file.path(projectPath,
                     outputPath,
                     "predictionCorrelations",
                     "reactogenicityReadoutsVisualisation",
                     "distribution_of_PTX3_per_vaccine_per_timepoint.png")
gp <- ggplot(data = toplot) +
  theme_bw(base_size = 16) +
  geom_boxplot(aes(x = treatment,
                   y = log2(labTestValue),
                   fill = treatment),
               colour = "black") +
  scale_fill_manual(values = mycols) +
  labs(x = "",
       y = "Log2(PTX3)",
       fill = "") +
  facet_wrap(~ timepoint) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))
ggsave(filename = figName,
       plot = gp,
       width = 10,
       height = 8)
rm(toplot, figName, gp, tps)


### Individual reactogenicity predictions per vaccine and per timepoint, for the chosen model(s)
tmp <- unique(reactoPredictions$tpTest)
tps <- humanMetadata %>%
  filter(timepoint %in% tmp) %>%
  dplyr::select(day, timepoint) %>%
  distinct() %>%
  arrange(day) %>%
  dplyr::select(timepoint) %>%
  unlist() %>%
  unname()
makeFigs <- lapply(unique(reactoPredictions$tpTrain), function(tp) {
  toplot <- reactoPredictions %>%
    filter(tpTrain == tp) %>%
    mutate(tpTest = factor(tpTest,
                           levels = tps))
  figName <- file.path(projectPath,
                       outputPath,
                       "predictionCorrelations",
                       "reactogenicityReadoutsVisualisation",
                       paste0("distribution_of_reacto_predictions_per_vaccine_per_timepoint_tpTrain",
                              tp,
                              ".png"))
  gp <- ggplot(data = toplot) +
    theme_bw(base_size = 16) +
    geom_boxplot(aes(x = treatment,
                     y = latent,
                     fill = treatment),
                 colour = "black") +
    scale_fill_manual(values = mycols) +
    labs(x = "",
         y = "Latent variable",
         fill = "") +
    facet_wrap(~ tpTest) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1))
  ggsave(filename = figName,
         plot = gp,
         width = 12,
         height = 6)
  rm(toplot, figName, gp)
})
rm(makeFigs, tmp, tps)


### Cleaning
rm(mycols)





####################################################################
### CORRELATION BETWWEEN REACTOGENICITY READOUTS AND PREDICTIONS ###
####################################################################

# Run only if necessary
if (runCorrelFigures)
{
  ### At the individual level
  
  # Parameters
  timepoints <- dataForCorrel %>%
    filter(timepoint >= 0) %>%
    dplyr::select(timepoint) %>%
    unlist() %>%
    unname() %>%
    unique()
  aeCats <- unique(aggregReactoscorePerVaccine$aeCategory)
  clinReadouts <- c(paste0("CRP_", timepoints),
                    paste0("PTX3_", timepoints),
                    paste0("reactoscore_", aeCats))
  clinReadouts <- c(clinReadouts,
                    paste0("log2_",
                           clinReadouts[!grepl("reactoscore", clinReadouts)]),
                    paste0("log2_",
                           clinReadouts[grepl("reactoscore", clinReadouts)],
                           "_plus1"))
  responseVar <- c("latent",
                   "log2_latent")
  separateVaccines <- c("all",
                        unique(dataForCorrel$treatment))
  excludeZeroReactoscores <- c(F, T)
  
  # Loop on response variables
  predClinCorrel <- lapply(responseVar, function(resp) {
    # Loop on clinical readouts
    figDoer <- lapply(clinReadouts, function(clinR) {
      
      # If looking at reactoscores, change name if figure if filtering subjective AEs or not
      figSuff <- case_when(
        str_detect(clinR, "reactoscore")
        & filterSubjectiveAEs == F ~ paste0(clinR, "_including_subjective"),
        str_detect(clinR, "reactoscore")
        & filterSubjectiveAEs == T ~ clinR,
        TRUE ~ clinR
      )
      
      # Different versions
      # => all vaccines together or vaccine per vaccine
      # => excluding or not zero values for reactoscore
      vers <- lapply(separateVaccines, function(sepVacc) {
        vers2 <- lapply(excludeZeroReactoscores, function(excl) {
          if (sepVacc == "all")
          {
            touse <- dataForCorrel
            # Figure name and and path
            figName <- file.path(projectPath,
                                 outputPath,
                                 "predictionCorrelations",
                                 modelToUse,
                                 "per_individual",
                                 paste0("INDIVIDUAL_ALLVACC_human_blood",
                                        ifelse(resp == "latent", "_", "_log2_"),
                                        "reactoG_predictions_vs_",
                                        figSuff,
                                        ".png"))
          } else
          {
            touse <- dataForCorrel %>%
              filter(treatment == sepVacc)
            # Figure name and and path
            figName <- file.path(projectPath,
                                 outputPath,
                                 "predictionCorrelations",
                                 modelToUse,
                                 "per_individual",
                                 paste0("INDIVIDUAL_",
                                        gsub(" ", "_", sepVacc),
                                        "_human_blood",
                                        ifelse(resp == "latent", "_", "_log2_"),
                                        "reactoG_predictions_vs_",
                                        figSuff,
                                        ".png"))
          }
          if (excl == T & grepl("reactoscore", clinR))
          {
            touse <- touse %>%
              filter(get(clinR) > 0)
            figName <- file.path(projectPath,
                                 outputPath,
                                 "predictionCorrelations",
                                 modelToUse,
                                 "per_individual",
                                 paste0("INDIVIDUAL_",
                                        ifelse(sepVacc == "all",
                                               "ALLVACC",
                                               gsub(" ", "_", sepVacc)),
                                        "_human_blood",
                                        ifelse(resp == "latent", "_", "_log2_"),
                                        "reactoG_predictions_vs_",
                                        figSuff,
                                        "_without_zeros",
                                        ".png"))
          }
          check <- touse %>%
            filter(!is.na(get(clinR))
                   & !is.na(get(resp)))
          if (nrow(check) > 0 & length(unique(check$participant)) > 2)
          {
            # Calculate the correlation coefficients (Spearman and Pearson)
            correlCoefs <- lapply(unique(touse$tpTrain), function(tpModel) {
              doer <- lapply(unique(touse$tpTest), function(tpObs) {
                tmp <- touse %>%
                  filter(tpTrain == tpModel
                         & tpTest == tpObs)
                x <- tmp %>%
                  dplyr::select(all_of(resp)) %>%
                  unlist() %>%
                  unname()
                y <- tmp %>%
                  dplyr::select(all_of(clinR)) %>%
                  unlist() %>%
                  unname()
                pcc <- cor.test(x,
                                y,
                                method = "pearson")$estimate
                scc <- cor.test(x,
                                y,
                                method = "spearman")$estimate
                res <- data.frame(tpTrain = tpModel,
                                  tpTest = tpObs,
                                  PCC = pcc,
                                  SCC = scc,
                                  stringsAsFactors = F)
                return(res)
              })
              doer <- do.call("bind_rows", doer)
              return(doer)
            })
            correlCoefs <- do.call("bind_rows", correlCoefs)
            
            # Data to plot
            toplot <- touse %>%
              left_join(correlCoefs,
                        by = c("tpTrain", "tpTest")) %>%
              # Position for correlation coefficients
              group_by(tpTrain, tpTest) %>%
              mutate(xPos = mean(get(resp), na.rm = T),
                     yPos = 0.85*max(get(clinR), na.rm = T)) %>%
              # Order of vaccines
              mutate(treatment = factor(treatment,
                                        levels = c("Saline",
                                                   "Engerix B",
                                                   "Agrippal",
                                                   "Varilrix",
                                                   "Stamaril",
                                                   "Fluad")))
            
            # Make the figures
            mycols <- colorsPerTreatment$human$colors
            names(mycols) <- colorsPerTreatment$human$treatments
            names(mycols)[names(mycols) == "Placebo"] <- "Saline"
            gp <- ggplot(data = toplot) +
              theme_bw(base_size = 28) +
              geom_point(aes(x = get(resp),
                             y = get(clinR),
                             colour = treatment),
                         size = 5) +
              geom_text(aes(x = xPos,
                            y = yPos,
                            label = paste0("SCC = ",
                                           round(SCC, digits = 2),
                                           "\n",
                                           "PCC = ",
                                           round(PCC, digits = 2))),
                        size = 8) +
              scale_colour_manual(values = mycols) +
              labs(x = ifelse(resp == "latent",
                              "Latent variable",
                              "Log2(latent variable)"),
                   y = Hmisc::capitalize(gsub("_", " ", clinR)),
                   colour = "") +
              facet_wrap(~ tpTrain + tpTest,
                         scales = "free_x",
                         labeller = labeller(tpTrain = ~ paste("MoMu: ", .),
                                             tpTest = ~ paste("HuBl: ", .),
                                             .multi_line = FALSE))
            ggsave(filename = figName,
                   plot = gp,
                   width = 24,
                   height = 18)
            
            # Output to return
            outp <- correlCoefs %>%
              mutate(clinicalReadout = clinR,
                     responseVariable = resp,
                     vaccinesIncluded = sepVacc,
                     excludeZeroValues = ifelse(grepl("reactoscore", clinR),
                                                excl,
                                                NA))
            return(outp)
          }
        })
        vers2 <- do.call("bind_rows", vers2)
        return(vers2)
      })
      vers <- do.call("bind_rows", vers)
      return(vers)
    })
    figDoer <- do.call("bind_rows", figDoer)
    return(figDoer)
  })
  predClinCorrel <- do.call("bind_rows", predClinCorrel)
  write.csv(predClinCorrel,
            file = file.path(projectPath,
                             outputPath,
                             "predictionCorrelations",
                             modelToUse,
                             paste0("INDIVIDUAL_correlations_human_blood_reactoG_predictions_vs_clinical_readouts",
                                    ifelse(filterSubjectiveAEs == F, "_including_subjective", ""),
                                    ".csv")),
            quote = F,
            row.names = F)
  rm(timepoints,
     aeCats,
     clinReadouts,
     responseVar,
     separateVaccines,
     excludeZeroReactoscores)
  
  
  ### At the vaccine level
  
  # Parameters
  timepoints <- dataForCorrel %>%
    filter(timepoint >= 0) %>%
    dplyr::select(timepoint) %>%
    unlist() %>%
    unname() %>%
    unique()
  aeCats <- unique(aggregReactoscorePerVaccine$aeCategory)
  clinReadouts <- c(paste0("mean_CRP_", timepoints),
                    paste0("mean_PTX3_", timepoints),
                    paste0("meanPatientReactoscore_", aeCats),
                    paste0("reactosum_", aeCats))
  clinReadouts <- c(clinReadouts,
                    paste0("log2_", clinReadouts))
  responseVar <- c("meanLatentVariable",
                   "log2_meanLatentVariable")
  
  # Loop on response variables
  aggregPredClinCorrel <- lapply(responseVar, function(resp) {
    # Loop on clinical readouts
    figDoer <- lapply(clinReadouts, function(clinR) {
      
      # If looking at reactoscores, change name if figure if filtering subjective AEs or not
      figSuff <- case_when(
        str_detect(clinR, "Reactoscore")
        & filterSubjectiveAEs == F ~ paste0(clinR, "_including_subjective"),
        str_detect(clinR, "Reactoscore")
        & filterSubjectiveAEs == T ~ clinR,
        TRUE ~ clinR
      )
      
      # Figure name and path
      figName <- file.path(projectPath,
                           outputPath,
                           "predictionCorrelations",
                           modelToUse,
                           "per_vaccine",
                           paste0("VACCINE_AGGREG_human_blood",
                                  ifelse(resp == "meanLatentVariable", "_", "_log2_"),
                                  "reactoG_predictions_vs_",
                                  figSuff,
                                  ".png"))
      
      # Calculate the correlation coefficients (Spearman and Pearson)
      correlCoefs <- lapply(unique(aggregDataForCorrel$tpTrain), function(tpModel) {
        doer <- lapply(unique(aggregDataForCorrel$tpTest), function(tpObs) {
          tmp <- aggregDataForCorrel %>%
            filter(tpTrain == tpModel
                   & tpTest == tpObs)
          x <- tmp %>%
            dplyr::select(all_of(resp)) %>%
            unlist() %>%
            unname()
          y <- tmp %>%
            dplyr::select(all_of(clinR)) %>%
            unlist() %>%
            unname()
          pcc <- cor.test(x,
                          y,
                          method = "pearson")$estimate
          scc <- cor.test(x,
                          y,
                          method = "spearman")$estimate
          res <- data.frame(tpTrain = tpModel,
                            tpTest = tpObs,
                            PCC = pcc,
                            SCC = scc,
                            stringsAsFactors = F)
          return(res)
        })
        doer <- do.call("bind_rows", doer)
        return(doer)
      })
      correlCoefs <- do.call("bind_rows", correlCoefs)
      
      # Data to plot
      toplot <- aggregDataForCorrel %>%
        left_join(correlCoefs,
                  by = c("tpTrain", "tpTest")) %>%
        # Position for correlation coefficients
        group_by(tpTrain, tpTest) %>%
        mutate(xPos = mean(get(resp), na.rm = T),
               yPos = 0.8*max(get(clinR), na.rm = T)) %>%
        # Order of vaccines
        mutate(treatment = factor(treatment,
                                  levels = c("Saline",
                                             "Engerix B",
                                             "Agrippal",
                                             "Varilrix",
                                             "Stamaril",
                                             "Fluad")))
      
      # Make the figures
      mycols <- colorsPerTreatment$human$colors
      names(mycols) <- colorsPerTreatment$human$treatments
      names(mycols)[names(mycols) == "Placebo"] <- "Saline"
      gp <- ggplot(data = toplot) +
        theme_bw(base_size = 28) +
        geom_point(aes(x = get(resp),
                       y = get(clinR),
                       colour = treatment),
                   size = 5) +
        geom_text(aes(x = xPos,
                      y = yPos,
                      label = paste0("SCC = ",
                                     round(SCC, digits = 2),
                                     "\n",
                                     "PCC = ",
                                     round(PCC, digits = 2))),
                  size = 8) +
        scale_colour_manual(values = mycols) +
        labs(x = ifelse(resp == "meanLatentVariable",
                        "Mean latent variable",
                        "Log2(mean latent variable)"),
             y = Hmisc::capitalize(gsub("_", " ", clinR)),
             colour = "") +
        facet_wrap(~ tpTrain + tpTest,
                   scales = "free_x",
                   labeller = labeller(tpTrain = ~ paste("MoMu: ", .),
                                       tpTest = ~ paste("HuBl: ", .),
                                       .multi_line = FALSE))
      ggsave(filename = figName,
             plot = gp,
             width = 24,
             height = 18)
      
      # Save correlations values
      outp <- correlCoefs %>%
        mutate(clinicalReadout = clinR,
               responseVariable = resp)
      return(outp)
    })
    figDoer <- do.call("bind_rows", figDoer)
    return(figDoer)
  })
  aggregPredClinCorrel <- do.call("bind_rows", aggregPredClinCorrel)
  write.csv(aggregPredClinCorrel,
            file = file.path(projectPath,
                             outputPath,
                             "predictionCorrelations",
                             modelToUse,
                             paste0("VACCINE_AGGREG_correlations_human_blood_reactoG_predictions_vs_clinical_readouts",
                                    ifelse(filterSubjectiveAEs == F, "_including_subjective", ""),
                                    ".csv")),
            quote = F,
            row.names = F)
  rm(timepoints,
     aeCats,
     clinReadouts,
     responseVar)
} else
{
  # Simply read the values of correlation coefficients
  predClinCorrel <- read.csv(file.path(projectPath,
                                      outputPath,
                                      "predictionCorrelations",
                                      modelToUse,
                                      paste0("INDIVIDUAL_correlations_human_blood_reactoG_predictions_vs_clinical_readouts",
                                             ifelse(filterSubjectiveAEs == F, "_including_subjective", ""),
                                             ".csv")),
                             header = T,
                             stringsAsFactors = F)
  aggregPredClinCorrel <- read.csv(file.path(projectPath,
                                             outputPath,
                                             "predictionCorrelations",
                                             modelToUse,
                                             paste0("VACCINE_AGGREG_correlations_human_blood_reactoG_predictions_vs_clinical_readouts",
                                                    ifelse(filterSubjectiveAEs == F, "_including_subjective", ""),
                                                    ".csv")),
                                   header = T,
                                   stringsAsFactors = F)
}


### Look at comparisons showing the highest correlation level
highCorrel <- predClinCorrel %>%
  filter(SCC >= 0.6
         | PCC >= 0.6) %>%
  filter(tpTrain != "72h"
         & vaccinesIncluded != "Saline") %>%
  arrange(tpTrain,
          desc(PCC),
          desc(SCC))
highAggregCorrel <- aggregPredClinCorrel %>%
  filter(SCC >= 0.6
         | PCC >= 0.6) %>%
  filter(tpTrain != "72h") %>%
  arrange(tpTrain,
          desc(PCC),
          desc(SCC))




###############################
### FIGURES FOR FINAL PAPER ###
###############################

### Figure 4 = boxplots of latent variable, reactoscore & CRP

# Palette
mycols <- colorsPerTreatment$human$colors
names(mycols) <- colorsPerTreatment$human$treatments
names(mycols)[names(mycols) == "Placebo"] <- "Saline"

# Boxplots of latent variable (showing all train and test timepoints)
toplot <- reactoPredictions %>%
  mutate(xlabel = "Mouse muscle",
         ylabel = "Human blood (2 classes)") %>%
  mutate(treatment = factor(treatment,
                            levels = c("Saline",
                                       "Engerix B",
                                       "Agrippal",
                                       "Varilrix",
                                       "Stamaril",
                                       "Fluad")))
bpLatent <- ggplot(data = toplot) +
  theme_bw(base_size = 16) +
  geom_boxplot(aes(x = treatment,
                   y = latent,
                   fill = treatment),
               colour = "black",
               outliers = F) +
  scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
  scale_fill_manual(values = mycols) +
  ggh4x::facet_nested(xlabel + tpTrain ~ ylabel + tpTest, scales = "free_y",  independent = "y") +
  labs(x = "",
       y = "Latent variable",
       fill = "") +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle=45,
                                   hjust=1))
rm(toplot)

# Boxplot of reactoscore
toplot <- individualReactoscores %>%
  mutate(treatment = factor(treatment,
                            levels = c("Saline",
                                       "Engerix B",
                                       "Agrippal",
                                       "Varilrix",
                                       "Stamaril",
                                       "Fluad")))
bpReacto <- ggplot(data = toplot) +
  theme_bw(base_size = 16) +
  geom_boxplot(aes(x = treatment,
                   y = log2(reactoscore + 1),
                   fill = treatment),
               colour = "black") +
  scale_fill_manual(values = mycols) +
  labs(x = "",
       y = "Log2(reactoscore + 1)",
       fill = "") +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle=45,
                                   hjust=1))
rm(toplot)

# Boxplots of CRP
tps <- unique(crpData$day)*24
tps <- tps[tps >= 0]
tps <- tps[order(tps)]
tps <- paste0(tps, "h")
toplot <- crpData %>%
  mutate(treatment = relevel(factor(treatment),
                             ref = "Saline")) %>%
  filter(day >= 1
         & day <= 3) %>%
  mutate(timepoint = factor(paste0(as.character(day*24),
                                   "h"),
                            levels = tps)) %>%
  mutate(treatment = factor(treatment,
                            levels = c("Saline",
                                       "Engerix B",
                                       "Agrippal",
                                       "Varilrix",
                                       "Stamaril",
                                       "Fluad")))
bpCRP <- ggplot(data = toplot) +
  theme_bw(base_size = 16) +
  geom_boxplot(aes(x = treatment,
                   y = log2(labTestValue),
                   fill = treatment),
               colour = "black") +
  scale_fill_manual(values = mycols) +
  labs(x = "",
       y = "Log2(CRP)",
       fill = "") +
  ggh4x::facet_nested(~ timepoint) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle=45,
                                   hjust=1))
rm(toplot)

# Arrange figures together
figbc <- ggarrange(plotlist = list(bpCRP, bpReacto),
                   nrow = 1,
                   ncol = 2,
                   widths = c(3, 1),
                   font.label = list(size = 24),
                   legend = "none",
                   labels = c("b", "c"))
figabc <- ggarrange(plotlist = list(bpLatent, figbc),
                    nrow = 2,
                    ncol = 1,
                    heights = c(2, 1),
                    font.label = list(size = 24),
                    legend = "top",
                    common.legend = T,
                    labels = c("a", ""))
figName <- file.path(projectPath,
                         outputPath,
                         "predictionCorrelations",
                         modelToUse,
                         paste0("Figure_4",
                                ifelse(filterSubjectiveAEs == F, "_including_subjective", ""),
                                ".pdf"))
ggsave(filename = figName,
       plot = figabc,
       width = 14,
       height = 12)
rm(figbc, bpLatent, bpReacto, bpCRP, figabc, figName)



### Supp figure 5 : heatmap of pvalues for latent variable (2 classes and 3 classes), CRP and reactoscore

# Palette
pvalueBreaks <- c(0,1e-3,1e-2,0.05,1e-1,1)
pvaluePalette <- RColorBrewer::brewer.pal(length(pvalueBreaks)-1, "RdYlBu")
pvaluePalette <- pvaluePalette[length(pvaluePalette):1]
names(pvaluePalette) <- levels(cut(c(0, 1), breaks = pvalueBreaks))

# Heatmap of pvalues - latent variables - 2 classes model
pvals <- lapply(unique(reactoPredictions$tpTrain), function(tpTrn) {
  res <- lapply(unique(reactoPredictions$tpTest), function(tpTst) {
    x <- reactoPredictions %>%
      filter(tpTest == tpTst
             & tpTrain == tpTrn) %>%
      mutate(treatment = factor(treatment,
                                 levels = c("Saline",
                                            "Engerix B",
                                            "Agrippal",
                                            "Varilrix",
                                            "Stamaril",
                                            "Fluad")))
    tmp <- pairwise.t.test(x$latent,
                           x$treatment,
                           na.rm = TRUE,
                           p.adjust.method = "BH")
    pvals <- reshape2::melt(tmp$p.value) %>%
      rename("treatment1" = "Var1",
             "treatment2" = "Var2",
             "adjPval" = "value") %>%
      filter(!is.na(adjPval)) %>%
      mutate(adjPvalCat = cut(adjPval,
                              breaks = pvalueBreaks)) %>%
      mutate(tpTest = tpTst,
             tpTrain = tpTrn)
    return(pvals)
  })
  res <- do.call("bind_rows", res)
  return(res)
})
pvals <- do.call("bind_rows", pvals) %>%
  mutate(xlabel = "Mouse muscle",
         ylabel = "Human blood (2 classes)")
phmLatent2 <- ggplot(data = pvals) +
  theme_bw(base_size = 16) +
  geom_tile(aes(x = treatment1,
                y = treatment2,
                fill = adjPvalCat)) +
  ggh4x::facet_nested(xlabel + tpTrain ~ ylabel + tpTest, scales = "free_y",  independent = "y") +
  labs(x = "",
       y = "",
       fill = "p-value",
       title = "Latent variable") +
  scale_fill_manual(values = pvaluePalette,
                    drop = F) +
  theme(axis.text.x = element_text(angle=45,
                                   hjust=1))
rm(pvals)


# Heatmap of pvalues - latent variables - 3 classes model
pvals <- lapply(unique(alternativeReactoPredictions$tpTrain), function(tpTrn) {
  res <- lapply(unique(alternativeReactoPredictions$tpTest), function(tpTst) {
    x <- alternativeReactoPredictions %>%
      filter(tpTest == tpTst
             & tpTrain == tpTrn) %>%
      mutate(treatment = factor(treatment,
                                levels = c("Saline",
                                           "Engerix B",
                                           "Agrippal",
                                           "Varilrix",
                                           "Stamaril",
                                           "Fluad")))
    tmp <- pairwise.t.test(x$latent,
                           x$treatment,
                           na.rm = TRUE,
                           p.adjust.method = "BH")
    pvals <- reshape2::melt(tmp$p.value) %>%
      rename("treatment1" = "Var1",
             "treatment2" = "Var2",
             "adjPval" = "value") %>%
      filter(!is.na(adjPval)) %>%
      mutate(adjPvalCat = cut(adjPval,
                              breaks = pvalueBreaks)) %>%
      mutate(tpTest = tpTst,
             tpTrain = tpTrn)
    return(pvals)
  })
  res <- do.call("bind_rows", res)
  return(res)
})
pvals <- do.call("bind_rows", pvals) %>%
  mutate(xlabel = "Mouse muscle",
         ylabel = "Human blood (3 classes)")
phmLatent3 <- ggplot(data = pvals) +
  theme_bw(base_size = 16) +
  geom_tile(aes(x = treatment1,
                y = treatment2,
                fill = adjPvalCat)) +
  ggh4x::facet_nested(xlabel + tpTrain ~ ylabel + tpTest, scales = "free_y",  independent = "y") +
  labs(x = "",
       y = "",
       fill = "p-value",
       title = "") +
  scale_fill_manual(values = pvaluePalette,
                    drop = F) +
  theme(axis.text.x = element_text(angle=45,
                                   hjust=1))
rm(pvals)

# Heatmap of pvalues - CRP
tps <- unique(crpData$day)*24
tps <- tps[tps >= 0]
tps <- tps[order(tps)]
tps <- paste0(tps, "h")
tempCrp <- crpData %>%
  filter(day >= 1
         & day <= 3) %>%
  mutate(timepoint = factor(paste0(as.character(day*24),
                                   "h"),
                            levels = tps)) %>%
  mutate(treatment = factor(treatment,
                            levels = c("Saline",
                                       "Engerix B",
                                       "Agrippal",
                                       "Varilrix",
                                       "Stamaril",
                                       "Fluad")))
pvals <- lapply(unique(tempCrp$timepoint), function(tp) {
  x <- tempCrp %>%
    filter(timepoint == tp)
  tmp <- pairwise.t.test(x$labTestValue,
                         x$treatment,
                         na.rm = TRUE,
                         p.adjust.method = "BH")
  pvals <- reshape2::melt(tmp$p.value) %>%
    rename("treatment1" = "Var1",
           "treatment2" = "Var2",
           "adjPval" = "value") %>%
    filter(!is.na(adjPval)) %>%
    mutate(adjPvalCat = cut(adjPval,
                            breaks = pvalueBreaks)) %>%
    mutate(timepoint = tp,
           xlabel = tp)
  return(pvals)
})
pvals <- do.call("bind_rows", pvals)
phmCRP <- ggplot(data = pvals) +
  theme_bw(base_size = 16) +
  geom_tile(aes(x = treatment1,
                y = treatment2,
                fill = adjPvalCat)) +
  ggh4x::facet_nested(~ xlabel) +
  labs(x = "",
       y = "",
       fill = "p-value",
       title = "CRP") +
  scale_fill_manual(values = pvaluePalette,
                    drop = F) +
  theme(axis.text.x = element_text(angle=45,
                                   hjust=1))
rm(tps, tempCrp, pvals)

# Heatmap of pvalues - reactoscore
individualReactoscores <- individualReactoscores %>%
  mutate(treatment = factor(treatment,
                            levels = c("Saline",
                                       "Engerix B",
                                       "Agrippal",
                                       "Varilrix",
                                       "Stamaril",
                                       "Fluad")))
tmp <- pairwise.t.test(individualReactoscores$reactoscore,
                       individualReactoscores$treatment,
                       na.rm = TRUE,
                       p.adjust.method = "BH")
pvals <- reshape2::melt(tmp$p.value) %>%
  rename("treatment1" = "Var1",
         "treatment2" = "Var2",
         "adjPval" = "value") %>%
  filter(!is.na(adjPval)) %>%
  mutate(adjPvalCat = cut(adjPval,
                          breaks = pvalueBreaks))
phmReacto <- ggplot(data = pvals) +
  theme_bw(base_size = 16) +
  geom_tile(aes(x = treatment1,
                y = treatment2,
                fill = adjPvalCat)) +
  labs(x = "",
       y = "",
       fill = "p-value",
       title = "Reactoscore") +
  scale_fill_manual(values = pvaluePalette,
                    drop = F) +
  theme(axis.text.x = element_text(angle=45,
                                   hjust=1))
rm(tmp, pvals)

# Arrange the figures together
figab <- ggarrange(plotlist = list(phmLatent2, phmLatent3),
                   nrow = 1,
                   ncol = 2,
                   font.label = list(size = 24),
                   legend = "top",
                   common.legend = T,
                   labels = c("a", "b"))
figcd <- ggarrange(plotlist = list(phmCRP, phmReacto),
                   nrow = 1,
                   ncol = 2,
                   widths = c(3, 1),
                   font.label = list(size = 24),
                   legend = "none",
                   labels = c("c", "d"))
figabcd <- ggarrange(plotlist = list(figab, figcd),
                   nrow = 2,
                   ncol = 1,
                   heights = c(1.5,1),
                   font.label = list(size = 24),
                   legend = "top",
                   common.legend = T,
                   labels = c("", ""))
suppfigName <- file.path(projectPath,
                      outputPath,
                      "predictionCorrelations",
                      modelToUse,
                      paste0("Supplementary_Figure_5",
                             ifelse(filterSubjectiveAEs == F, "_including_subjective", ""),
                             ".pdf"))
ggsave(filename = suppfigName,
       plot = figabcd,
       width = 18,
       height = 10)
rm(figab, figcd, figabcd, suppfigName)





### OLD CODE FOR FIGURES
# # Fig 4: reactoscore
# myReadouts <- c("meanPatientReactoscore_all")
# myResponse <- "meanLatentVariable"
# 
# listFig4 <- lapply(myReadouts, function(readOut) {
#   toplot1 <- aggregPredClinCorrel %>%
#     filter(tpTrain == myTrainTP
#            & clinicalReadout == readOut
#            & responseVariable == myResponse) 
#   toplot2 <- aggregDataForCorrel %>%
#     filter(tpTrain == myTrainTP) %>%
#     dplyr::select(tpTrain, tpTest, treatment,
#                   all_of(c(readOut, myResponse))) %>%
#     # Order of vaccines
#     mutate(treatment = factor(treatment,
#                               levels = c("Placebo",
#                                          "Engerix B",
#                                          "Agrippal",
#                                          "Varilrix",
#                                          "Stamaril",
#                                          "Fluad")))
#   toplot <- toplot2 %>%
#     left_join(toplot1,
#               by = c("tpTrain", "tpTest")) %>%
#     # Position for correlation coefficients
#     group_by(tpTrain, tpTest) %>%
#     mutate(xPos = 1.05*min(get(myResponse), na.rm = T),
#            yPos = 0.9*max(get(readOut), na.rm = T))
#   rm(toplot1, toplot2)
#   gp <- ggplot(data = toplot) +
#     theme_bw(base_size = 24) +
#     geom_point(aes(x = get(myResponse),
#                    y = get(readOut),
#                    colour = treatment),
#                size = 5) +
#     geom_text(aes(x = xPos,
#                   y = yPos,
#                   label = paste0("SCC = ",
#                                  round(SCC, digits = 2),
#                                  "\n",
#                                  "PCC = ",
#                                  round(PCC, digits = 2))),
#               size = 6) +
#     scale_colour_manual(values = mycols) +
#     labs(x = ifelse(myResponse == "meanLatentVariable",
#                     "Mean latent variable",
#                     "Log2(mean latent variable)"),
#          y = case_when(
#            readOut == "meanPatientReactoscore_all" ~ "Mean reactoscore",
#            readOut == "meanPatientReactoscore_systemic" ~ "Mean reactoscore (systemic)",
#            readOut == "meanPatientReactoscore_local" ~ "Mean reactoscore (local)"
#            ),
#          colour = "") +
#     facet_wrap(~ tpTrain + tpTest,
#                scales = "free_x",
#                labeller = labeller(tpTrain = ~ paste("MoMu: ", .),
#                                    tpTest = ~ paste("HuBl: ", .),
#                                    .multi_line = FALSE))
#   return(gp)
# })
# fig4 <- ggpubr::ggarrange(plotlist = listFig4,
#                           ncol = 1,
#                           nrow = length(listFig4),
#                           font.label = list(size = 24),
#                           legend = "right",
#                           common.legend = T)
# fig4Name <- file.path(projectPath,
#                      outputPath,
#                      "predictionCorrelations",
#                      modelToUse,
#                      paste0("Figure_4",
#                             ifelse(filterSubjectiveAEs == F, "_including_subjective", ""),
#                             ".pdf"))
# ggsave(filename = fig4Name,
#        plot = fig4,
#        width = 24,
#        height = 6)
# rm(myReadouts,
#    myResponse,
#    listFig4,
#    fig4,
#    fig4Name)
# 
# # Fig 5: CRP
# myReadouts <- c("mean_CRP_24h",
#                 "mean_CRP_48h",
#                 "mean_CRP_72h")
# myResponse <- "meanLatentVariable"
# 
# listFig5 <- lapply(myReadouts, function(readOut) {
#   toplot1 <- aggregPredClinCorrel %>%
#     filter(tpTrain == myTrainTP
#            & clinicalReadout == readOut
#            & responseVariable == myResponse) 
#   toplot2 <- aggregDataForCorrel %>%
#     filter(tpTrain == myTrainTP) %>%
#     dplyr::select(tpTrain, tpTest, treatment,
#                   all_of(c(readOut, myResponse))) %>%
#     # Order of vaccines
#     mutate(treatment = factor(treatment,
#                               levels = c("Placebo",
#                                          "Engerix B",
#                                          "Agrippal",
#                                          "Varilrix",
#                                          "Stamaril",
#                                          "Fluad")))
#   toplot <- toplot2 %>%
#     left_join(toplot1,
#               by = c("tpTrain", "tpTest")) %>%
#     # Position for correlation coefficients
#     group_by(tpTrain, tpTest) %>%
#     mutate(xPos = 1.05*min(get(myResponse), na.rm = T),
#            yPos = 0.8*max(get(readOut), na.rm = T))
#   rm(toplot1, toplot2)
#   gp <- ggplot(data = toplot) +
#     theme_bw(base_size = 24) +
#     geom_point(aes(x = get(myResponse),
#                    y = get(readOut),
#                    colour = treatment),
#                size = 5) +
#     geom_text(aes(x = xPos,
#                   y = yPos,
#                   label = paste0("SCC = ",
#                                  round(SCC, digits = 2),
#                                  "\n",
#                                  "PCC = ",
#                                  round(PCC, digits = 2))),
#               size = 6) +
#     scale_colour_manual(values = mycols) +
#     labs(x = ifelse(myResponse == "meanLatentVariable",
#                     "Mean latent variable",
#                     "Log2(mean latent variable)"),
#          y = Hmisc::capitalize(gsub("_", " ", readOut)),
#          colour = "") +
#     facet_wrap(~ tpTrain + tpTest,
#                scales = "free_x",
#                labeller = labeller(tpTrain = ~ paste("MoMu: ", .),
#                                    tpTest = ~ paste("HuBl: ", .),
#                                    .multi_line = FALSE))
#   return(gp)
# })
# fig5 <- ggpubr::ggarrange(plotlist = listFig5,
#                           ncol = 1,
#                           nrow = length(listFig5),
#                           labels = c("a)", "b)", "c)"),
#                           font.label = list(size = 24),
#                           legend = "right",
#                           common.legend = T)
# fig5Name <- file.path(projectPath,
#                       outputPath,
#                       "predictionCorrelations",
#                       modelToUse,
#                       "Figure_5.pdf")
# ggsave(filename = fig5Name,
#        plot = fig5,
#        width = 24,
#        height = 18)
# rm(myReadouts,
#    myResponse,
#    listFig5,
#    fig5,
#    fig5Name)
# 
# # Fig 6: PTX (TBD)
# myReadouts <- c("mean_PTX3_24h",
#                 "mean_PTX3_48h",
#                 "mean_PTX3_72h")
# myResponse <- "meanLatentVariable"
# 
# listfig6 <- lapply(myReadouts, function(readOut) {
#   toplot1 <- aggregPredClinCorrel %>%
#     filter(tpTrain == myTrainTP
#            & clinicalReadout == readOut
#            & responseVariable == myResponse) 
#   toplot2 <- aggregDataForCorrel %>%
#     filter(tpTrain == myTrainTP) %>%
#     dplyr::select(tpTrain, tpTest, treatment,
#                   all_of(c(readOut, myResponse))) %>%
#     # Order of vaccines
#     mutate(treatment = factor(treatment,
#                               levels = c("Placebo",
#                                          "Engerix B",
#                                          "Agrippal",
#                                          "Varilrix",
#                                          "Stamaril",
#                                          "Fluad")))
#   toplot <- toplot2 %>%
#     left_join(toplot1,
#               by = c("tpTrain", "tpTest")) %>%
#     # Position for correlation coefficients
#     group_by(tpTrain, tpTest) %>%
#     mutate(xPos = 1.05*min(get(myResponse), na.rm = T),
#            yPos = 0.95*max(get(readOut), na.rm = T))
#   rm(toplot1, toplot2)
#   gp <- ggplot(data = toplot) +
#     theme_bw(base_size = 24) +
#     geom_point(aes(x = get(myResponse),
#                    y = get(readOut),
#                    colour = treatment),
#                size = 5) +
#     geom_text(aes(x = xPos,
#                   y = yPos,
#                   label = paste0("SCC = ",
#                                  round(SCC, digits = 2),
#                                  "\n",
#                                  "PCC = ",
#                                  round(PCC, digits = 2))),
#               size = 6) +
#     scale_colour_manual(values = mycols) +
#     labs(x = ifelse(myResponse == "meanLatentVariable",
#                     "Mean latent variable",
#                     "Log2(mean latent variable)"),
#          y = Hmisc::capitalize(gsub("_", " ", readOut)),
#          colour = "") +
#     facet_wrap(~ tpTrain + tpTest,
#                scales = "free_x",
#                labeller = labeller(tpTrain = ~ paste("MoMu: ", .),
#                                    tpTest = ~ paste("HuBl: ", .),
#                                    .multi_line = FALSE))
#   return(gp)
# })
# fig6 <- ggpubr::ggarrange(plotlist = listfig6,
#                           ncol = 1,
#                           nrow = length(listfig6),
#                           labels = c("a)", "b)", "c)"),
#                           font.label = list(size = 24),
#                           legend = "right",
#                           common.legend = T)
# fig6Name <- file.path(projectPath,
#                       outputPath,
#                       "predictionCorrelations",
#                       modelToUse,
#                       "Figure_6.pdf")
# ggsave(filename = fig6Name,
#        plot = fig6,
#        width = 24,
#        height = 18)
# rm(myReadouts,
#    myResponse,
#    listfig6,
#    fig6,
#    fig6Name)


# ### Supplementary figure for the paper -> reactoscore + CRP
# 
# # Aggregated correlations at the vaccine level
# # Show mouse muscle 24h, human blood 24-48-72h, latent variable VS
# # -> reactoscore
# # -> CRP 24-48h (remove 72h)
# # -> PTX3 ? NO
# 
# # Common parameters
# mycols <- colorsPerTreatment$human$colors
# names(mycols) <- colorsPerTreatment$human$treatments
# names(mycols)[names(mycols) == "Placebo"] <- "Saline"
# myTrainTP <- "72h"
# myTestTP <- c("24h", "48h")
# 
# myReadouts <- c("meanPatientReactoscore_all",
#                 "mean_CRP_24h",
#                 "mean_CRP_48h",
#                 "mean_CRP_72h")
# myResponse <- "meanLatentVariable"
# 
# dataForFig4 <- lapply(myReadouts, function(readOut) {
#   toplot1 <- aggregPredClinCorrel %>%
#     filter(tpTrain == myTrainTP
#            & tpTest %in% myTestTP
#            & clinicalReadout == readOut
#            & responseVariable == myResponse)
#   toplot2 <- aggregDataForCorrel %>%
#     filter(tpTrain == myTrainTP
#            & tpTest %in% myTestTP) %>%
#     dplyr::select(tpTrain, tpTest, treatment,
#                   all_of(c(readOut, myResponse))) %>%
#     # Order of vaccines
#     mutate(treatment = factor(treatment,
#                               levels = c("Saline",
#                                          "Engerix B",
#                                          "Agrippal",
#                                          "Varilrix",
#                                          "Stamaril",
#                                          "Fluad")))
#   toplot <- toplot2 %>%
#     left_join(toplot1,
#               by = c("tpTrain", "tpTest")) %>%
#     # Rename columns for easier call in facet plot
#     rename("response" = myResponse,
#            "readout" = readOut) %>%
#     # Define nested labels for facet figure
#     mutate(xlabel = if_else(readOut == "meanPatientReactoscore_all",
#                             "Reactoscore",
#                             Hmisc::capitalize(gsub("_", " ", gsub("mean_", "", readOut))))) %>%
#     mutate(ylabel = "Human blood")
#   rm(toplot1, toplot2)
#   
#   return(toplot)
# })
# dataForFig4 <- do.call("bind_rows", dataForFig4) %>%
#   mutate(xlabel = relevel(as.factor(xlabel), ref = "Reactoscore"))
# # Part 1: reactoscore/CRP vs latent variable
# gp <- ggplot(data = dataForFig4) +
#   theme_bw(base_size = 24) +
#   geom_point(aes(x = readout,
#                  y = response,
#                  colour = treatment),
#              size = 5) +
#   scale_colour_manual(values = mycols) +
#   labs(y = ifelse(myResponse == "meanLatentVariable",
#                   "Mean latent variable",
#                   "Log2(mean latent variable)"),
#        x = "Mean reactoscore / CRP",
#        colour = "") +
#   ggh4x::facet_nested(tpTest ~ xlabel,
#                       scales = "free", 
#                       independent = "all")
# # Part 2: reactoscore VS CRP
# x <- dataForFig4 %>%
#   filter(clinicalReadout == "meanPatientReactoscore_all") %>%
#   dplyr::select(-all_of(c("response",
#                           "responseVariable",
#                           "PCC",
#                           "SCC",
#                           "ylabel",
#                           "tpTrain",
#                           "tpTest"))) %>%
#   rename("response" = "readout",
#          "responseVariable" = "clinicalReadout",
#          "ylabel" = "xlabel") %>%
#   distinct()
# y <- dataForFig4 %>%
#   filter(clinicalReadout != "meanPatientReactoscore_all") %>%
#   dplyr::select(-all_of(c("response",
#                           "responseVariable",
#                           "PCC",
#                           "SCC",
#                           "ylabel",
#                           "tpTrain",
#                           "tpTest"))) %>%
#   distinct()
# dataForFig4_v2 <- y %>%
#   left_join(x,
#             by = "treatment") %>%
#   mutate(xlabel = factor(as.character(xlabel)))
# rm(x, y)
# gp2 <- ggplot(data = dataForFig4_v2) +
#   theme_bw(base_size = 24) +
#   geom_point(aes(x = readout,
#                  y = response,
#                  colour = treatment),
#              size = 5) +
#   scale_colour_manual(values = mycols) +
#   labs(y = "Mean reactoscore",
#        x = "Mean CRP",
#        colour = "") +
#   ggh4x::facet_nested(~ xlabel,
#                       scales = "free", 
#                       independent = "all")
# gpblank <- ggplot() + geom_blank() + theme_void()
# gp2 <- ggarrange(plotlist = list(gpblank, gp2),
#                  nrow = 1,
#                  ncol = 2,
#                  widths = c(1, 3),
#                  font.label = list(size = 24),
#                  legend = "none",
#                  labels = c("", "b"))
# fig4 <- ggarrange(plotlist = list(gp, gp2),
#                   ncol = 1,
#                   nrow = 2,
#                   heights = c(2, 1),
#                   font.label = list(size = 24),
#                   legend = "right",
#                   common.legend = T,
#                   labels = c("a", ""))
# fig4Name <- file.path(projectPath,
#                       outputPath,
#                       "predictionCorrelations",
#                       modelToUse,
#                       paste0("OLD_Figure_4",
#                              ifelse(filterSubjectiveAEs == F, "_including_subjective", ""),
#                              ".pdf"))
# ggsave(filename = fig4Name,
#        plot = fig4,
#        width = 24,
#        height = 14)
# rm(myReadouts,
#    myResponse,
#    dataForFig4,
#    gp,
#    gp2,
#    fig4,
#    fig4Name,
#    mycols,
#    myTrainTP,
#    myTestTP)





### Cleaning
rm(list = ls())
gc()


