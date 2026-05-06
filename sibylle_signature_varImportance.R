####################################################
### LOAD FUNCTIONS, PARAMETERS DATA AND METADATA ###
####################################################


setwd("/sps/bioaster/Projects/SIBYLLE_BA039/BA039-WP4-XDATA-02/SIBYLLE_2025_clean_pipeline_for_paper/")
source("src/sibylle_parameters.R")
source("src/sibylle_functions.R")
source("src/sibylle_data_loader.R")


##############################################################################
### FIT AND EVALUATE GLMNET IN MOUSE MUSCLE (100 NESTED CROSS VALIDATIONS) ###
##############################################################################


probeSubset = "mouse"
nClasses = 3
lambdaFromMoBl = F
reweight = T

#for(probeSubset in c("mouse","human")){
    #for(nClasses in 2:3){
        #for(reweight in c(F,T)){
            print(paste(nClasses, reweight, probeSubset))
            p_metrics = fit_evaluate_glmnet_in_mouseMuscle(nClasses, reweight, probeSubset, lambdaFromMoBl)
        #}
    #}
#}

stabilityRes = run_stability_selection_in_mouseMuscle(nClasses, reweight, probeSubset="mouse")
figure2(p_metrics, stabilityRes)
Supplementary_Figure_3(probeSubset, nClasses, reweight, stabilityRes[["freqSelection"]], degFilename["mouseMuscle"])


############################################################################
### TRAIN glmnetcr IN MOUSE MUSCLE AND EVALUATE IN MOUSE AND HUMAN BLOOD ###
############################################################################


lambdaOptim = T
#nClasses = 2
probeSubset = "mouse"
pairwiseTtest_LatentVar = list(list())
# for(probeSubset in c("mouse","human")){
    reweight = ifelse(probeSubset=="mouse",T,F)
    for(nClasses in 2:3){
        try(pairwiseTtest_LatentVar[[as.character(nClasses)]][[probeSubset]] <- train_glmnet_in_mouseMuscle_evaluate_in_blood.wrap(nClasses, reweight, probeSubset, lambdaOptim))            
    }
# }

Figure_3(pairwiseTtest_LatentVar[["3"]], stabilityRes[["topGenes"]])
Supplementary_Figure_4(degFilename, stabilityRes[["freqSelection"]], pairwiseTtest_LatentVar)


###########################################################
### BOXPLOT ACCURACY ACROSS MODELS (= MODEL COMPARISON) ###
###########################################################


boxplot_accuracy_across_models()

