######################
### LOAD LIBRARIES ###
######################


library(ggplot2)
library(preprocessCore)
library(matrixStats)
library(sva)
library(limma)
library(mltools)
library(parallel)
library(caret)
library(dplyr)
library(glmnetcr)
library(rsample)
library(ggh4x)
library(reshape2)
library(openxlsx)
library(patchwork)
library(ggnewscale)
library(tidyr)
library(ggpubr)
library(clusterProfiler)
library(msigdbr)
library(cluster)
library(mclust)


#######################################
### LOAD DATA AND GENERATE METADATA ###
#######################################


# Load orthology data
orthology = read.table("Analyses/orthology/mouse_human_orthology_only_unique_pairs.csv",sep=",",h=T)
human2mouse = orthology$probeName.mouse
names(human2mouse) = orthology$probeName.human

# Load mouse array annotations
arrayMap = read.table("Data/GPL10787-9758.txt",h=T,sep="\t",stringsAsFactors=F, fill = TRUE,  comment.char = "#")
probe2geneName = arrayMap[,"GENE_SYMBOL"]
names(probe2geneName) = arrayMap[,"ID"]

# Load msigdb hallmark pathways
msigdb = msigdbr(species = "Mus musculus", category = "H")

# Load expression and metadata data
load("Data/metadata_data_preprocessed.RData")

# for(i in 1:length(DSname)){
#     data[[DSname$name[i]]] = NULL
#     try(data[[DSname$name[i]]] <- read.table(
#         file=paste0("Analyses/biovacsafe_",DSname$species[i],"/QC/biovacsafe_",DSname$species[i],"_count_data_filtering_low_count_genes_",DSname$compartment[i],".csv"),
#         sep=",",h=T
#     ))
#     if(!is.null(data[[DSname$name[i]]])){

#         # Correct TP=-24h in sample name
#         colnames(data[[DSname$name[i]]]) = sub("\\.24h","-24h",colnames(data[[DSname$name[i]]]))

#         # Get rownames from column X
#         rownames(data[[DSname$name[i]]]) = data[[DSname$name[i]]][,"X"]
#         data[[DSname$name[i]]] = as.matrix(data[[DSname$name[i]]][,-1])
#         data[[DSname$name[i]]] = data[[DSname$name[i]]][,!colnames(data[[DSname$name[i]]])%in%outliers]

#         # Generate metadata from colnames
#         metadata[[DSname$name[i]]] = as.data.frame(t(sapply(colnames(data[[DSname$name[i]]]), function(x) strsplit(x, '_')[[1]][1:3])))
#         colnames(metadata[[DSname$name[i]]]) = c("vaccine","compartment","timePoint")
#         metadata[[DSname$name[i]]]$reactoClass = vaccine2ReactoClasses[metadata[[DSname$name[i]]]$vaccine]
#         metadata[[DSname$name[i]]]$vaccineTP = paste(metadata[[DSname$name[i]]]$vaccine,metadata[[DSname$name[i]]]$timePoint,sep="_")
#     }
# }


#########################
### TP=0h CORRECTION  ###
#########################


# for(name in DSname$name){
#     data[[paste0(name,"_D0")]] = matrix(NA, nrow=nrow(data[[name]]), ncol=ncol(data[[name]]))
#     colnames(data[[paste0(name,"_D0")]]) = colnames(data[[name]])
#     rownames(data[[paste0(name,"_D0")]]) = rownames(data[[name]])
#     for(vaccine in unique(metadata[[name]]$vaccine)){
#         D0 = rowMedians(as.matrix(data[[name]][,grep(paste0(vaccine,ifelse(name=="humanBlood",".*_-24h_",".*_0h_")),colnames(data[[name]]),value=T)]))
#         iVaccine = grep(vaccine,colnames(data[[name]]))
#         data[[paste0(name,"_D0")]][,iVaccine] = as.matrix(data[[name]][,iVaccine]) - D0
#     }
#     metadata[[paste0(name,"_D0")]] = metadata[[name]]
# }
#
# save(metadata, data, file="Data/metadata_data_preprocessed.RData")


#########################################################
### GET SHARED PROBES ACROSS MOUSE AND HUMAN DATASETS ###
#########################################################


probeMap = intersect(rownames(data[["humanBlood"]]),names(human2mouse))
names(probeMap) = human2mouse[probeMap]
probeMap = probeMap[intersect(rownames(data[["mouseMuscle"]]),names(probeMap))]


##############################
### REMOVE DATANAME != _D0 ###
##############################


data = data[grep("_D0",names(data),value=T)]
metadata = metadata[grep("_D0",names(data),value=T)]


##################################################
### REPLACE VACCINE ABBREVIATION BY THEIR NAME ###
##################################################


for(name in names(data)){
    colnames(data[[name]]) = str_replace_all(colnames(data[[name]]), abbreviation2Name)
    metadata[[name]]$vaccine = str_replace_all(metadata[[name]]$vaccine, abbreviation2Name)
}


##################################################################
### GENERATE sampleTP AND reactoClasses FOR ALL DATASET AND TP ###
##################################################################


sampleTP = reactoClasses = NULL
for(name in names(data)){
    species = sub("(mouse|human).*","\\1",name)
    for(tp in TPs){
        sampleTP[[name]][[tp]][[3]] = grep(paste0("_",tp),colnames(data[[name]]), value=T)
        reactoClasses[[name]][[tp]][[3]] = factor(vaccine2ReactoClasses[sub("^(.*)_(MU|HU|BL)_.*$","\\1",sampleTP[[name]][[tp]][[3]])],levels=levels[[species]][[3]], ordered=T)
        
        if(grepl("mouse",name)){
            i2Classes = which(reactoClasses[[name]][[tp]][[3]]!="medium")
        }else{i2Classes = 1:length(reactoClasses[[name]][[tp]][[3]])}
        sampleTP[[name]][[tp]][[2]] = sampleTP[[name]][[tp]][[3]][i2Classes]
        reactoClasses[[name]][[tp]][[2]] = factor(reactoClasses[[name]][[tp]][[3]][i2Classes],levels=levels[[species]][[2]], ordered=T)
    }
}

