################################################################
### SIBYLLE PAPER - REDO DIFFERENTIAL ANALYSIS ON MOUSE DATA ###
################################################################


### Libraries
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(edgeR)


### Parameters and paths are defined in the YAML file
params <- yaml::read_yaml(file = file.path(this.path::here(),
                                           "sibylle_mouse_human_data_preprocessing_parameters.yaml"))
suppressWarnings(params %>% attach) # To enable accessing the parameters by their name without using params$...
tissues <- list(mouse = tissuesOfInterest$mouse,
                human = tissuesOfInterest$human)
# Set appropriate working directory
setwd(projectPath)


### Should low variance filtering be applied?
lowVarFiltering <- F


### Create appropriate directory
dir.create(file.path(projectPath,
                     outputPath,
                     mousePath,
                     "DEG"),
           showWarnings = F)


### Load mouse data

# Metadata
mouseMeta <- read.table(file.path(projectPath,
                                  dataPath,
                                  mousePath,
                                  "biovacsafe_mouse_metadata_for_analysis.txt"),
                        sep = "\t",
                        header = T,
                        stringsAsFactors = F) %>%
  mutate(treatment_timepoint = paste0(gsub(" ", "", treatment),
                                      "_",
                                      timepoint),
         reactoClass_tp = paste0(reactoClass,
                                 "_",
                                 timepoint))

# Count data in muscle
mouseCountMuscle <- read.table(file.path(projectPath,
                                         dataPath,
                                         mousePath,
                                         "biovacsafe_mouse_muscle_count_data_for_analysis.txt"),
                               sep = "\t",
                               header = T,
                               stringsAsFactors = F)
mouseCountBlood <- read.table(file.path(projectPath,
                                        dataPath,
                                        mousePath,
                                        "biovacsafe_mouse_blood_count_data_for_analysis.txt"),
                              sep = "\t",
                              header = T,
                              stringsAsFactors = F)
mouseCountAll <- list(muscle = mouseCountMuscle,
                      blood = mouseCountBlood)


### Filtering of low count genes if need be
if (lowVarFiltering == T)
{
  # Criterion: keep only genes with > X counts in at least K samples, with K = nb samples in smallest tissue-TP-treatment group
  # => see https://f1000research.com/articles/5-1384"
  genesToKeepPerTissue <- lapply(c("muscle", "blood"), function(tiss) {
    curData <- mouseCountAll[[tiss]]
    sds <- matrixStats::rowSds(as.matrix(curData))
    sds <- sds[order(sds, decreasing = T)]
    genesToKeep <- names(sds[(1:5000)])
    return(genesToKeep)
  })
  names(genesToKeepPerTissue) <- c("muscle", "blood")
  
  # Filter low count genes
  mouseCountAllFiltered <- lapply(names(mouseCountAll), function(tiss) {
    # Genes to keep for considered tissue
    genesToKeep <- genesToKeepPerTissue[[tiss]]
    # Samples for considered tissue
    sampTissue <- mouseMeta %>%
      filter(tissue == tiss) %>%
      dplyr::select(sampleID) %>%
      unlist() %>%
      unname()
    # Filtered data
    dat <- mouseCountAll[[tiss]][genesToKeep, sampTissue]
    return(dat)
  })
  names(mouseCountAllFiltered) <- names(mouseCountAll)
} else
{
  mouseCountAllFiltered <- mouseCountAll
}


### Do differential analysis - comparing vaccines

doDEG <- lapply(c("muscle" ,"blood"), function(tiss) {
  
  # Arrange data and metadata in the same order
  mouseCount <- mouseCountAllFiltered[[tiss]]
  curMouseMeta <- mouseMeta %>%
    filter(sampleID %in% colnames(mouseCount))
  mouseCount <- mouseCount[,curMouseMeta$sampleID]
  
  # Design matrix
  design <- model.matrix(~ 0 + treatment_timepoint,
                         curMouseMeta)
  
  # Contrasts = (treatment_tp - treatment_0h) - (saline_tp - saline_0h)
  tps <- unique(curMouseMeta$timepoint)
  tps <- tps[tps != "0h"]
  treats <- gsub(" ", "", unique(curMouseMeta$treatment))
  treats <- treats[treats != "Saline"]
  contrast0 <- unlist(lapply(tps, function(t) {
    unlist(lapply(treats, function(g) {
      paste0("(treatment_timepoint", g, "_", t,
             "-treatment_timepoint", g, "_0h)",
             "-(",
             "treatment_timepointSaline_", t,
             "-treatment_timepointSaline_0h)")
    }))
  }))
  contrast <- makeContrasts(contrasts = contrast0,
                            levels = design)
  nb.contrast <- dim(contrast)[2]
  rm(tps, treats)
  
  # Limma fit with all contrasts
  fit.ebayes <- eBayes(contrasts.fit(lmFit(mouseCount,
                                           design = design),
                                     contrasts = contrast),
                       robust = T,
                       trend = T)
  
  # Table of DEG results, with treatment and timepoint info
  toptableRes <- lapply(colnames(fit.ebayes$coefficients), function(ctr) {
    # Retrieve treatment and timepoint info
    tmp <- unlist(str_split(gsub("treatment_timepoint", "", ctr), pattern = "\\)-\\("))
    tmp <- gsub("\\(", "", gsub("\\(", "", tmp))
    tmp <- strsplit(tmp, split = "-")[[1]]
    treat <- strsplit(tmp, split = "_")[[1]][1]
    tp <- strsplit(tmp, split = "_")[[1]][2]
    rm(tmp)
    res <- topTable(fit.ebayes,
                    coef = ctr,
                    adjust = "BH",
                    n = Inf,
                    confint = T)
    res <- res %>%
      mutate(contrastFormula = ctr,
             treatment = treat,
             timepoint = tp) %>%
      rownames_to_column(var = "geneID")
    return(res)
  })
  names(toptableRes) <- colnames(fit.ebayes$coefficients)
  toptableResDf <- do.call("bind_rows", toptableRes) %>%
    # Add tissue info
    mutate(tissueType = tiss)
  write.xlsx(toptableResDf,
             file.path(projectPath,
                       outputPath,
                       mousePath,
                       "DEG",
                       paste0("biovacsafe_mouse_differential_analysis_vaccine_",
                              tiss,
                              "_results_all_genes.xlsx")),
            quote = F,
            rowNames = F)
  
  
  ### Extract DEG = genes with pvalue < 0.01 (same criterion as in McKay original paper)
  degRes <- toptableResDf %>%
    filter(adj.P.Val < 0.01)
  write.xlsx(degRes,
            file.path(projectPath,
                      outputPath,
                      mousePath,
                      "DEG",
                      paste0("biovacsafe_mouse_differential_analysis_vaccine_",
                             tiss,
                             "_results_only_DEG.xlsx")),
            quote = F,
            rowNames = F)
})
rm(doDEG)


### Do differential analysis - comparing reacto classes (low VS high)

doDEG <- lapply(c("muscle" ,"blood"), function(tiss) {
  # Arrange data and metadata in the same order
  mouseCount <- mouseCountAllFiltered[[tiss]]
  curMouseMeta <- mouseMeta %>%
    filter(sampleID %in% colnames(mouseCount))
  mouseCount <- mouseCount[,curMouseMeta$sampleID]
  
  # Design matrix
  design <- model.matrix(~ 0 + reactoClass_tp,
                         curMouseMeta)
  
  # Contrasts
  tps <- unique(curMouseMeta$timepoint)
  tps <- tps[tps != "0h"]
  treats <- gsub(" ", "", unique(curMouseMeta$reactoClass))
  treats <- treats[treats != "low"]
  contrast0 <- unlist(lapply(tps, function(t) {
    unlist(lapply(treats, function(g) {
      paste0("(reactoClass_tp", g, "_", t,
             "-reactoClass_tp", g, "_0h)",
             "-(",
             "reactoClass_tplow_", t,
             "-reactoClass_tplow_0h)")
    }))
  }))
  contrast <- makeContrasts(contrasts = contrast0,
                            levels = design)
  nb.contrast <- dim(contrast)[2]
  rm(tps)
  
  # Limma fit with all contrasts
  fit.ebayes <- eBayes(contrasts.fit(lmFit(mouseCount,
                                           design = design),
                                     contrasts = contrast),
                       robust = T,
                       trend = F)
  
  # Table of DEG results, with reactoclass and timepoint info
  toptableRes <- lapply(colnames(fit.ebayes$coefficients), function(ctr) {
    # Retrieve reactoclass and timepoint info
    tmp <- unlist(str_split(gsub("reactoClass_tp", "", ctr), pattern = "\\)-\\("))
    tmp <- gsub("\\(", "", gsub("\\(", "", tmp))
    tmp <- strsplit(tmp, split = "-")[[1]]
    reactoclass <- strsplit(tmp, split = "_")[[1]][1]
    tp <- strsplit(tmp, split = "_")[[1]][2]
    rm(tmp)
    res <- topTable(fit.ebayes,
                    coef = ctr,
                    adjust = "BH",
                    n = Inf,
                    confint = T)
    res <- res %>%
      mutate(contrastFormula = ctr,
             reactoClass = reactoclass,
             timepoint = tp) %>%
      rownames_to_column(var = "geneID")
    return(res)
  })
  names(toptableRes) <- colnames(fit.ebayes$coefficients)
  toptableResDf <- do.call("bind_rows", toptableRes) %>%
    # Add tissue info
    mutate(tissueType = tiss) %>%
    # Limit to the comparison high vs low
    filter(!str_detect(contrastFormula, "medium"))
  write.xlsx(toptableResDf,
             file.path(projectPath,
                       outputPath,
                       mousePath,
                       "DEG",
                       paste0("biovacsafe_mouse_differential_analysis_reactoclass_",
                              tiss,
                              "_results_all_genes.xlsx")),
             quote = F,
             rowNames = F)
  
  
  ### Extract DEG = genes with pvalue < 0.01 (same criterion as in McKay original paper)
  degRes <- toptableResDf %>%
    filter(adj.P.Val < 0.01)
  write.xlsx(degRes,
             file.path(projectPath,
                       outputPath,
                       mousePath,
                       "DEG",
                       paste0("biovacsafe_mouse_differential_analysis_reactoclass_",
                              tiss,
                              "_results_only_DEG.xlsx")),
             quote = F,
             rowNames = F)
})
rm(doDEG)

# Quick check to validate pvalues
tiss <- "muscle"
gene <- "A_55_P1965467"
todo <- lapply(c("0h", "168h"), function(tp) {
  selSampHigh <- mouseMeta %>%
    filter(tissue == tiss
           & timepoint == tp
           & reactoClass == "high")
  selSampLow <- mouseMeta %>%
    filter(tissue == tiss
           & timepoint == tp
           & reactoClass == "low")
  selCountHigh <- t(mouseCountAllFiltered[[tiss]][gene,selSampHigh$sampleID])
  selCountLow <- t(mouseCountAllFiltered[[tiss]][gene,selSampLow$sampleID])
  toplot1 <- data.frame(geneID = gene,
                        reactoClass = "high",
                        timepoint = tp,
                        sampleID = rownames(selCountHigh),
                        geneExpr = selCountHigh)
  toplot2 <- data.frame(geneID = gene,
                        reactoClass = "low",
                        timepoint = tp,
                        sampleID = rownames(selCountLow),
                        geneExpr = selCountLow)
  toplot <- toplot1 %>%
    bind_rows(toplot2) %>%
    rename("geneExpr" = gene)
  rm(toplot1, toplot2)
  return(toplot)
})
toplot <- do.call("bind_rows", todo)
rm(todo)
# forstat <- toplot %>%
#   mutate(newSampleID = str_remove(str_remove(sampleID, "_0h"), "_4h")) %>%
#   pivot_wider(id_cols = all_of(c("newSampleID", "geneID", "reactoClass")),
#               names_from = timepoint,
#               values_from = geneExpr,
#               names_prefix = "geneExpr_") %>%
#   mutate(geneExprDiff = geneExpr_4h - geneExpr_0h)
tmp1 <- sapply(toplot$geneExpr[toplot$timepoint == "4h" & toplot$reactoClass == "high"],
              "-",
              toplot$geneExpr[toplot$timepoint == "0h" & toplot$reactoClass == "high"])
tmp1 <- as.vector(tmp1)
tmp2 <- sapply(toplot$geneExpr[toplot$timepoint == "4h" & toplot$reactoClass == "low"],
               "-",
               toplot$geneExpr[toplot$timepoint == "0h" & toplot$reactoClass == "low"])
tmp2 <- as.vector(tmp2)
forstat <- data.frame(reactoClass = c(rep("high", length(tmp1)),
                                      rep("low", length(tmp2))),
                      geneExprDiff = c(tmp1,
                                       tmp2),
                      stringsAsFactors = F)
rm(tmp1, tmp2)
ttest <- ggpubr::compare_means(geneExprDiff ~ reactoClass,
                               data = forstat,
                               method = "t.test",
                               paired = F,
                               p.adjust.method = "BH")
gp <- ggplot(data = toplot) +
  theme_bw(base_size = 16) +
  geom_boxplot(aes(x = reactoClass,
                   y = geneExpr,
                   fill = timepoint)) +
  labs(x = "Reactogenicity class",
       y = "Gene expression",
       fill = "TP",
       title = gene,
       subtitle = paste0("pval = ", ttest$p.adj))
plot(gp)
rm(tiss,
   gene,
   tp,
   selSampHigh,
   selSampLow,
   selCountLow,
   selCountHigh,
   toplot,
   gp)



### Cleaning
rm(list = ls())
gc()

