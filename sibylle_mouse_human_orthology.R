######################################################################
### SIBYLLE HUMAN STUDY (blood only) - ORTHOLOGY WITH MOUSE STUDY  ###
######################################################################




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


### Create folder for orthology
dir.create(file.path(projectPath,
                     outputPath,
                     "orthology"),
           showWarnings = F)




#########################################
### LOAD FEATURES FOR MOUSE AND HUMAN ###
#########################################

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




#######################################################
### MATCH HUMAN (blood only) GENES with MOUSE GENES ###
#######################################################


### Genes in human (blood only) study -> features not filtered out and protein coding
featsToKeepNoDuplikProtCoding <- read.xlsx(file.path(projectPath,
                                                     dataPath,
                                                     humanPath,
                                                     "biovacsafe_human_features_protein_coding_without_duplicates_for_analysis.xlsx"),
                                           sheet = 1)
humanGenesMapping <- featsToKeepNoDuplikProtCoding
# Add mapping to HGNC (GENE_SYMBOL is supposed to be HGNCSymbol most of the time, but HGNC mapping can bring additional info)
conv <- limma::alias2SymbolUsingNCBI(alias = humanGenesMapping$GENE_SYMBOL,
                                     gene.info.file = file.path(projectPath,
                                                                dataPath,
                                                                humanPath,
                                                                "Homo_sapiens.gene_info"))
colnames(conv) <- c("HGNCID",
                    "HGNCSymbol",
                    "description")
humanGenesMapping <- cbind(humanGenesMapping,
                           conv)
rm(conv)
write.xlsx(humanGenesMapping,
           file.path(projectPath,
                     dataPath,
                     humanPath,
                     "biovacsafe_human_mapping_to_hgnc.xlsx"))


### Load probes information with chromosomal location details

# Mouse
mouseChrLoc <- read.table(file.path(projectPath,
                                    dataPath,
                                    mousePath,
                                    "biovacsafe_mouse_features_chromosomal_location_with_probe_ids.txt"),
                          sep = "\t",
                          stringsAsFactors = F,
                          header = T) %>%
  # Limit to probes included in mouse study
  filter(SPOT_ID %in% mouseFeatures$SPOT_ID) %>%
  # Reorder start and end
  mutate(startCorr = if_else(end < start, end, start),
         endCorr = if_else(end < start, start, end)) %>%
  mutate(strand = if_else(start != startCorr & end != endCorr,
                          "-",
                          "+")) %>%
  dplyr::select(-all_of(c("start", "end"))) %>%
  rename(start = startCorr,
         end = endCorr) %>%
  # Chromosomal location with start and end in the proper order and without the -1 correction needed for LiftOver
  mutate(chrLoc.mouse = paste0(chr, ":", start, "-", end)) %>%
  # Exclude probes with unmapped chromosomal location
  filter(chr != "unmapped")

# Human
humanChrLoc <- read.table(file.path(projectPath,
                                    dataPath,
                                    humanPath,
                                    "biovacsafe_human_features_chromosomal_location_with_probe_ids.txt"),
                          sep = "\t",
                          stringsAsFactors = F,
                          header = T) %>%
  # Limit to probes included in human study
  filter(SPOT_ID %in% featsToKeepNoDuplikProtCoding$SPOT_ID) %>%
  # Reorder start and end
  mutate(startCorr = if_else(end < start, end, start),
         endCorr = if_else(end < start, start, end)) %>%
  mutate(strand = if_else(start != startCorr & end != endCorr,
                          "-",
                          "+")) %>%
  dplyr::select(-all_of(c("start", "end"))) %>%
  rename(start = startCorr,
         end = endCorr) %>%
  # Chromosomal location with start and end in the proper order and without the -1 correction needed for LiftOver
  mutate(chrLoc.human = paste0(chr, ":", start, "-", end)) %>%
  # Exclude probes with unmapped chromosomal location
  filter(chr != "unmapped")


### Load LiftOver (UCSC) match and annotation between mouse and human chromosomal locations

# NB: some chromosomal locations in mouse have multiple matches in human genome

# Load full match
chrLocMatchesFull <- read.table(file.path(projectPath,
                                          dataPath,
                                          "liftover_chromosomal_location_match_and_annotation_from_mouse_ncbi37_to_human_hg19.bed"),
                                sep = "\t",
                                header = F,
                                stringsAsFactors = F)
colnames(chrLocMatchesFull) <- c("chr.human",
                                 "start.human",
                                 "end.human",
                                 "chrLoc.mouse",
                                 "matchIdx")
chrLocMatchesFull <- chrLocMatchesFull %>%
  # Remove duplicated lines (some mouse probes have same chromosomal location, generated duplicated lines in LiftOver output)
  distinct() %>%
  # if multiple human chromosomal locations matching one mouse chromosomal location, keep only the first one
  filter(matchIdx == 1) %>%
  # Reformat chromosomal position information
  mutate(chrLoc.human = paste0(chr.human,
                               ":",
                               start.human,
                               "-",
                               end.human)) %>%
  mutate(positionForLiftOver = str_replace(str_replace(chrLoc.mouse, ":", " "),
                                           "-", " ")) %>%
  mutate(chr.mouse = str_split(positionForLiftOver, " ", simplify = T)[,1],
         start.mouse = as.integer(str_split(positionForLiftOver, " ", simplify = T)[,2]),
         end.mouse = as.integer(str_split(positionForLiftOver, " ", simplify = T)[,3])) %>%
  # Rearrange columns for use
  dplyr::select(chrLoc.mouse,
                chr.mouse,
                start.mouse,
                end.mouse,
                positionForLiftOver,
                chrLoc.human,
                chr.human,
                start.human,
                end.human,
                matchIdx)

# Limit to transcripts of interest in mouse study
chrLocMatches <- chrLocMatchesFull %>%
  mutate(positionForLiftOver = paste(chr.mouse,
                                     start.mouse,
                                     end.mouse)) %>%
  filter(chrLoc.mouse %in% mouseChrLoc$chrLoc.mouse) %>%
  dplyr::select(-all_of(c("positionForLiftOver"))) %>%
  inner_join(mouseChrLoc, by = "chrLoc.mouse") %>%
  dplyr::select(-all_of(c("positionForLiftOver",
                          "chr",
                          "start",
                          "end"))) %>%
  rename(probeName.mouse = SPOT_ID,
         geneSymbol.mouse = GENE_SYMBOL,
         ensemblTranscriptID.mouse = ENSEMBL_ID,
         chromosomalLocation.mouse = CHROMOSOMAL_LOCATION,
         refseq.mouse = REFSEQ,
         strand.mouse = strand)

# NB : nrow(chrLocMatches) < nrow(mouseChrLoc) because
# - some transcripts in mouseChrLoc had chr.mouse == "unmapped" so no match
# - somes transcripts with valid chromosomal positions had no or insufficient match (6055)
# => see liftover_match_mouse_ncbi37_human_hg19_failure_file.txt


### Identify human matches to mouse chromosomal locations in BioVacSafe human data
# => validating based on the exon positions

# Load genome with exons details
gtf <- rtracklayer::import("/sps/bioaster/Data/Genomes/H.sapiens/hg19/GCF_000001405.25_GRCh37.p13_genomic.gtf")
gtf <- as.data.frame(gtf)
# Load correspondence between REFSEQ and chromosome number (from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/)
x <- read.table(file.path(projectPath,
                          dataPath,
                          humanPath,
                          "GCF_000001405.25_GRCh37.p13_assembly_report.txt"),
                sep = "\t",
                header = F,
                stringsAsFactors = F)
colnames(x) <- c("Sequence-Name",
                 "Sequence-Role",
                 "Assigned-Molecule",
                 "Assigned-Molecule-Location-Type",
                 "GenBank-Accn",
                 "Relationship",
                 "RefSeq-Accn",
                 "Assembly-Unit",
                 "Sequence-Length",
                 "UCSC-style-name")
# Extract the exons on human genome (NB: start < end for all lines here, not inverted depending on strand)
allExons <- gtf %>%
  filter(type == "exon")
tmp <- allExons$transcript_id
tmp <- str_split(tmp, pattern = "\\.", simplify = T)[,1]
allExons$transcript_id_without_version <- tmp
rm(tmp)
# Add the chromosome number
tmp <- x %>%
  dplyr::select("RefSeq-Accn", "UCSC-style-name") %>%
  dplyr::rename("seqnames" = "RefSeq-Accn",
                "chromosomeUCSC" = "UCSC-style-name")
allExons <- allExons %>%
  left_join(tmp, by = "seqnames")
rm(tmp)
# Limit to exons with the same transcript ID (RefSeq) as in BioVacSafe human data
exons <- allExons %>%
  filter(transcript_id_without_version %in% featsToKeepNoDuplikProtCoding$REFSEQ)

# For each chromosomal location from BioVacSafe human, find the intersecting exons with same strand as transcript
exonsHumanChrLoc <- lapply((1:nrow(humanChrLoc)), function(k) {
  curHumanProbes <- humanChrLoc[k,]
  curExons <- exons %>%
    filter(chromosomeUCSC == curHumanProbes$chr) %>%
    filter((start >= curHumanProbes$start
            & end <= curHumanProbes$end) # overlap on the left
           | (start >= curHumanProbes$start
              & start <= curHumanProbes$end
              & end >= curHumanProbes$end) # overlap on the right
           | (start >= curHumanProbes$start
              & start <= curHumanProbes$end
              & end <= curHumanProbes$end) # overlap in the middle
           | (start <= curHumanProbes$start
              & end >= curHumanProbes$end)) %>% # exon wider than gene
    filter(as.character(strand) == curHumanProbes$strand) %>%
    dplyr::select(seqnames,
                  start,
                  end,
                  strand,
                  type,
                  transcript_id,
                  transcript_id_without_version,
                  exon_number,
                  chromosomeUCSC) %>%
    rename(exonRefSeq = seqnames,
           exonStart = start,
           exonEnd = end,
           exonStrand = strand,
           exonType = type,
           exonTranscriptID = transcript_id,
           exonTranscriptIDWithoutVersion = transcript_id_without_version,
           exonNumber = exon_number)
  if (nrow(curExons) > 0)
  {
    res <- cbind(curExons, curHumanProbes) %>%
      mutate(nbExonsMatched = nrow(curExons))
    return(res)
  }
})
exonsHumanChrLoc <- do.call("bind_rows", exonsHumanChrLoc)

# For each chromosomal location from mouse-human LiftOver match, find the intersecting exons
exonsChrLocMatches <- lapply((1:nrow(chrLocMatches)), function(k) {
  curLocMatches <- chrLocMatches[k,]
  curExons <- exons %>%
    filter(chromosomeUCSC == curLocMatches$chr.human) %>%
    filter((start >= curLocMatches$start.human
            & end <= curLocMatches$end.human) # overlap on the left
           | (start >= curLocMatches$start.human
              & start <= curLocMatches$end.human
              & end >= curLocMatches$end.human) # overlap on the right
           | (start >= curLocMatches$start.human
              & start <= curLocMatches$end.human
              & end <= curLocMatches$end.human) # overlap in the middle
           | (start <= curLocMatches$start.human
              & end >= curLocMatches$end.human)) %>% # exon wider than gene
    # NB: we do not have the strand information for chromosomal locations matched from LiftOver
    dplyr::select(seqnames,
                  start,
                  end,
                  strand,
                  type,
                  transcript_id,
                  transcript_id_without_version,
                  exon_number,
                  chromosomeUCSC) %>%
    rename(exonRefSeq = seqnames,
           exonStart = start,
           exonEnd = end,
           exonStrand = strand,
           exonType = type,
           exonTranscriptID = transcript_id,
           exonTranscriptIDWithoutVersion = transcript_id_without_version,
           exonNumber = exon_number)
  if (nrow(curExons) > 0)
  {
    res <- cbind(curExons, curLocMatches) %>%
      mutate(nbExonsMatched = nrow(curExons))
    return(res)
  }
})
exonsChrLocMatches <- do.call("bind_rows", exonsChrLocMatches)
save.image(file.path(projectPath,
                     outputPath,
                     "orthology",
                     "mouse_human_orthology_intermediate_workspace.RData"))

# For each chromosomal location matched with BioVacSafe mouse probes, find the BioVacSafe human probe that correspond
# load(paste0(projectPath,
#                     outputPath,
#                     "orthology",,
#             "mouse_human_orthology_intermediate_workspace.RData"))
chrLocMatchesHumanMouse <- lapply((1:nrow(chrLocMatches)), function(k) {
  curLocMatch <- chrLocMatches[k,]
  curExons <- exonsChrLocMatches %>%
    filter(chr.human == curLocMatch$chr.human
           & start.human == curLocMatch$start.human
           & end.human == curLocMatch$end.human)
  curExonsMatch <- exonsHumanChrLoc %>%
    filter(exonRefSeq %in% curExons$exonRefSeq
           & exonStart %in% curExons$exonStart
           & exonEnd %in% curExons$exonEnd
           & exonTranscriptID %in% curExons$exonTranscriptID
           & chromosomeUCSC %in% curExons$chromosomeUCSC)
  if (nrow(curExonsMatch) > 0)
  {
    res <- curExonsMatch %>%
      dplyr::select(SPOT_ID,
                    GENE_SYMBOL,
                    ENSEMBL_ID,
                    CHROMOSOMAL_LOCATION,
                    REFSEQ) %>%
      distinct() %>%
      mutate(nbExonsInLiftOverChromosomalLocation = nrow(curExons),
             nbMatchedExons = nrow(curExonsMatch))
    res <- res %>%
      bind_cols(curLocMatch)
    return(res)
  }
})
chrLocMatchesHumanMouse <- do.call("bind_rows", chrLocMatchesHumanMouse) %>%
  rename(probeName.human = SPOT_ID,
         geneSymbol.human = GENE_SYMBOL,
         ensemblTranscriptID.human = ENSEMBL_ID,
         chromosomalLocation.human = CHROMOSOMAL_LOCATION,
         refseq.human = REFSEQ)
# Arrange columns to have all human fields aside and all mouse fields aside
clnms <- colnames(chrLocMatchesHumanMouse)
hClnms <- clnms[grepl("human", tolower(clnms))]
mClnms <- clnms[grepl("mouse", tolower(clnms))]
chrLocMatchesHumanMouse <- chrLocMatchesHumanMouse[,c(mClnms,
                                                      hClnms,
                                                      clnms[!(clnms %in% c(mClnms,
                                                                           hClnms))])]
rm(clnms, hClnms, mClnms)

# Verify content of match
nrow(chrLocMatchesHumanMouse) # 17560
verif <- chrLocMatchesHumanMouse %>%
  filter(!is.na(probeName.human) & !is.na(probeName.mouse)) %>%
  dplyr::select(probeName.mouse, probeName.human) %>%
  distinct()
nrow(verif) # 17560 => no bizarre match with no probe associated
length(unique(verif$probeName.human)) # 13056 => human probes matched to multiple mouse probes
length(unique(verif$probeName.mouse)) # 14841 => mouse probes matched to multiple human probes
rm(verif)

# Investigate multiple matches and single matches
mouseHumanPairs <- chrLocMatchesHumanMouse %>%
  filter(!is.na(probeName.human) & !is.na(probeName.mouse)) %>%
  dplyr::select(probeName.mouse, probeName.human) %>%
  distinct() # 17560 rows
humanToManyMouse <- mouseHumanPairs %>%
  group_by(probeName.human) %>%
  summarise(nMouseProbesMatching = n()) %>%
  filter(nMouseProbesMatching > 1)
mouseToManyHuman <- mouseHumanPairs %>%
  group_by(probeName.mouse) %>%
  summarise(nHumanProbesMatching = n()) %>%
  filter(nHumanProbesMatching > 1)
multMatches <- chrLocMatchesHumanMouse %>%
  filter(probeName.mouse %in% mouseToManyHuman$probeName.mouse
         | probeName.human %in% humanToManyMouse$probeName.human) %>%
  mutate(multipleMatchType = case_when(
    probeName.mouse %in% mouseToManyHuman$probeName.mouse
    & probeName.human %in% humanToManyMouse$probeName.human ~ "multiple match both in mouse and human",
    probeName.mouse %in% mouseToManyHuman$probeName.mouse
    & !(probeName.human %in% humanToManyMouse$probeName.human) ~ "mouse probe matching multiple human probes",
    !(probeName.mouse %in% mouseToManyHuman$probeName.mouse)
    & probeName.human %in% humanToManyMouse$probeName.human ~ "human probe matching multiple mouse probes"
  ))
plyr::count(multMatches$multipleMatchType)
singleMatches <- chrLocMatchesHumanMouse %>%
  filter(!(probeName.mouse %in% mouseToManyHuman$probeName.mouse)
         & !(probeName.human %in% humanToManyMouse$probeName.human))

# Summary table of numbers // match content verification
myTab <- matrix("",
                nrow = 4,
                ncol = 5)
myTab[1,] <- c("",
               "All pairs",
               "Mouse probes mapping to multiple human probes",
               "Human probes mapping to multiple mouse probes",
               "Single matches only")
myTab[,1] <- c("",
               "# rows",
               "# unique mouse probes",
               "# unique human probes")
fillTab <- lapply(list(chrLocMatchesHumanMouse,
                       chrLocMatchesHumanMouse %>%
                         filter(probeName.mouse %in% mouseToManyHuman$probeName.mouse),
                       chrLocMatchesHumanMouse %>%
                         filter(probeName.human %in% humanToManyMouse$probeName.human),
                       singleMatches),
                  function(x) {
                    res <- c(nrow(x),
                             length(unique(x$probeName.mouse)),
                             length(unique(x$probeName.human)))
                    return(res)
                  })
for (k in (1:length(fillTab)))
{
  myTab[(2:nrow(myTab)), k+1] <- fillTab[[k]]
}
rm(k)
myTab <- myTab %>%
  as.data.frame() %>%
  column_to_rownames(var = "V1") %>%
  janitor::row_to_names(1)
png(file.path(projectPath,
              outputPath,
              "orthology",
              "mouse_human_orthology_summary.png"),
    width = 1000,
    height = 320)
plot(1,
     type="n",
     xlab="",
     ylab="",
     main = "Summary - mouse to human orthology",
     sub = "",
     xaxt = "n",
     yaxt = "n",
     bty = "n")
gridExtra::grid.table(myTab,
                      theme = gridExtra::ttheme_minimal())
dev.off()
rm(myTab, fillTab)

# Save content of match
save.image(file.path(projectPath,
                     outputPath,
                     "orthology",
                     "mouse_human_orthology_output.RData"))
write.csv(chrLocMatchesHumanMouse,
          file = file.path(projectPath,
                           outputPath,
                           "orthology",
                           "mouse_human_orthology_with_one_to_many_matches.csv"),
          quote = F,
          row.names = F)
write.csv(singleMatches,
          file = file.path(projectPath,
                           outputPath,
                           "orthology",
                           "mouse_human_orthology_only_unique_pairs.csv"),
          quote = F,
          row.names = F)

# Other quick checks on gene symbol
x <- chrLocMatchesHumanMouse %>%
  mutate(geneSymbolMatching = if_else(tolower(geneSymbol.mouse) == tolower(geneSymbol.human),
                                      T, F))
plyr::count(x$geneSymbolMatching)
# FALSE 2510
# TRUE 7357
# => in most cases, gene symbols match
x %>% filter(geneSymbolMatching == F) %>% dplyr::select(geneSymbol.mouse, geneSymbol.human, everything()) %>% View()
rm(x)
# Same checks but on single matches
x <- singleMatches %>%
  mutate(geneSymbolMatching = if_else(tolower(geneSymbol.mouse) == tolower(geneSymbol.human),
                                      T, F))
plyr::count(x$geneSymbolMatching)
# FALSE 723
# TRUE 4038
# => in most cases, gene symbols match
x %>% filter(geneSymbolMatching == F) %>% dplyr::select(geneSymbol.mouse, geneSymbol.human, everything()) %>% View()
# => often gene symbols do not strictly match but are similar (e.g. Zfp212 - ZFP212)
rm(x)

# Cleaning
rm(list = ls())
gc()


