###########################
### FUNCTION DEFINITION ###
###########################


weighted_f1 <- function(data, reference){

    confMat = confusionMatrix(data = data, reference = reference)
    byClass = confMat$byClass

    if (is.null(dim(byClass))) {
        # Binary classification: byClass is a named vector
        precision = byClass["Pos Pred Value"]
        recall = byClass["Sensitivity"]
        f1 = 2 * precision * recall / (precision + recall)
        f1[is.na(f1)] = 0
        return(list(weighted_f1=as.numeric(f1), confMat=confMat))

    } else {
        
        # Multi-class: byClass is a matrix
        precision = byClass[, "Pos Pred Value"]
        recall = byClass[, "Sensitivity"]
        f1_per_class = 2 * precision * recall / (precision + recall)
        f1_per_class[is.na(f1_per_class)] = 0

        # Align support with the order of classes in byClass
        support = table(reference)
        common_levels = sub("Class: ","",rownames(byClass))
        support = support[common_levels]  # reorders correctly

        weighted_f1 = sum(f1_per_class * support / sum(support))
        return(list(weighted_f1=weighted_f1, confMat=confMat$table))
    }
}




glmModel <- function(x, y, lambda, maxit, weights){

    # fit model : glmnetcr or glmnet
    if(nClasses > 2){
        glmnetcr(y = y, x = x, lambda = lambda, maxit = maxit, weights = weights)
    }else{glmnet(y = y, x = x, lambda = lambda, maxit = maxit, weights = weights, family='binomial')}

}


glmModel.wrap <- function(y, x, lambdaSeq, nSelectedVar, returnLambda){
    
    # Define weights
    if(reweight){
        class_weights = 1/table(y)
        weights = as.numeric(class_weights[y])
    }else{weights = rep(1,length(y))/length(y)}

    # Fit model to get lambda range
    if(all(is.na(lambdaSeq))){
        mod = glmModel(x, y, NULL, maxit, weights)
        lambda_max = max(mod$lambda)
        lambda_min = min(mod$lambda)
        lambdaSeq = exp(seq(log(lambda_max), log(lambda_min), length.out = nLambda))
        if(returnLambda){return(lambdaSeq)}
    }

    # Fit model on the full lambda range
    mod = glmModel(x, y, lambdaSeq, maxit, weights)
    nselectedvar = computeSelectedVar(mod)
    
    if(!is.na(nSelectedVar)){
        iBestLambda = which.min(abs(nSelectedVar - nselectedvar[["n"]]))
        mod = glmModel(x, y, lambdaSeq[iBestLambda], maxit, weights)
        nselectedvar = computeSelectedVar(mod)
    }

    return(list(mod = mod, nSelectedVar = nselectedvar[["n"]], selected = nselectedvar[["selected"]]))
}


computeSelectedVar <- function(mod){

    beta_names = rownames(mod$beta)
    valid_idx = grep("^cp[12]", beta_names, invert = TRUE)
    res = logical(length(valid_idx))
    
    if(ncol(mod$beta)==1){
        selected = which(abs(mod$beta[valid_idx,])>sqrt(.Machine$double.eps))
        res[selected] = T
        n = length(selected)
    }else{
        n = apply(mod$beta[valid_idx,] ,2, function(x) {
            length(which(abs(x)>sqrt(.Machine$double.eps)))
        })
        selected = NULL
    }
    list(selected=res, n=n)
}


compute_phat_glmnet <- function(x, y, B = 50, nSelectedVar, lambdaSeq = NA) {

    p = ncol(x)
    n = nrow(x)
    phat_array = array(FALSE, dim = c(p, length(lambdaSeq), B),
                    dimnames = list(colnames(x), paste0("s",seq_len(length(lambdaSeq))-1), paste0("B", 1:B)))

    selection_list  = mclapply(1:B, function(b) {

        # Subsampling
        idx = sample(1:n, size = floor(0.5 * n), replace = FALSE)

        # Run model wrapper
        mod = glmModel.wrap(y = y[idx], x = x[idx, ], lambdaSeq = lambdaSeq, nSelectedVar = nSelectedVar, returnLambda = FALSE)

        # Build full selection matrix with all p variables
        beta_path = mod$mod$beta[grep("^cp[12]", rownames(mod$mod$beta), invert = TRUE), , drop = FALSE]
        return((abs(beta_path) > sqrt(.Machine$double.eps)))
    }, mc.cores=nbcores)

    # Compute element-wise mean across matrices
    phat = as.matrix(Reduce("+", selection_list) / length(selection_list))

    return(phat)
}


cv.glmnetcr <- function(name, tp, probenames, lambdaSeq){

    # Get reactoclasses and sampletp
    reactoclasses = reactoClasses[[name]][[tp]][[nClasses]]
    sampletp = sampleTP[[name]][[tp]][[nClasses]]


    # Get data
    dat = data[[name]][probenames,]

    # Outer CV: estimate generalization performance
    outerFolds = createFolds(reactoclasses, k = nFolds, returnTrain=T)
    latent = yOutterPred = numeric(length(reactoclasses))

    # Fit model to get lambda range
    if(is.na(lambdaSeq)){
        lambdaSeq = glmModel.wrap(
            y = reactoclasses,
            x = t(dat[,sampletp]),
            lambdaSeq = NA,
            nSelectedVar = NA,
            returnLambda = T
        )
    }

    # Initialization
    outterBestNSelectedVar = NULL

    for (i in seq_along(outerFolds)) {

        # Inner CV: tune lambda
        innerFolds = createFolds(reactoclasses[outerFolds[[i]]], k = nFolds, returnTrain = TRUE)
        innerAccuracy = innerBestNSelectedVar = NULL

        for (j in seq_along(innerFolds)) {
            
            # Fit inner model
            innerModel = glmModel.wrap(
                y = reactoclasses[outerFolds[[i]]][innerFolds[[j]]],
                x = t(dat[,sampletp[outerFolds[[i]]][innerFolds[[j]]]]),
                lambdaSeq = lambdaSeq,
                nSelectedVar = NA,
                returnLambda = F)
            innerBestNSelectedVar = rbind(innerBestNSelectedVar, innerModel$nSelectedVar)

            # Predict
            yInnerPred = predict(
                innerModel$mod,
                t(dat[,sampletp[outerFolds[[i]]][-innerFolds[[j]]]]), 
                type = "class")
            if(nClasses>2){yInnerPred = yInnerPred$class}

            innerAccuracy = rbind(
                innerAccuracy,
                apply(yInnerPred, 2, function(y_pred) as.character(reactoclasses[outerFolds[[i]]][-innerFolds[[j]]]) == y_pred)
            )
        }

        # Select best lambda
        outterBestNSelectedVar[i] = max(innerBestNSelectedVar[,which.max(colMeans(innerAccuracy))])

        # Retrain model on full outer training set with best lambda
        outterModel = glmModel.wrap(y = reactoclasses[outerFolds[[i]]],
                            x = t(dat[,sampletp[outerFolds[[i]]]]),
                            lambdaSeq = lambdaSeq,
                            nSelectedVar = outterBestNSelectedVar[i],
                            returnLambda = F
                        )

        # Predict on outer test set
        youtterpred = predict(outterModel$mod, t(dat[,sampletp[-outerFolds[[i]]]]), type = "class")
        if(nClasses>2){
            yOutterPred[-outerFolds[[i]]] = youtterpred$class
        }else{yOutterPred[-outerFolds[[i]]] = youtterpred}
    }

    weighted_f1 = weighted_f1(data=factor(yOutterPred, levels=levels[[sub(".*(mouse|human).*","\\1",name)]][[nClasses]]), reference=reactoclasses)

    # Return metrics
    return(
        list(
            confMat = weighted_f1$confMat,
            metrics = data.frame(
                name = name, 
                tp = tp,
                nSelectedVar = median(outterBestNSelectedVar),
                weighted_f1 = weighted_f1$weighted_f1
            )
        )
    )
}

plotMetrics <- function(metrics){

    if("tp"%in%colnames(metrics)){
        metrics$tp = factor(metrics$tp,levels=TPs)
    }
    metrics$ylabel = "Mouse Blood"
    
    # plot prediction performances
    if("tp"%in%colnames(metrics)){
        p <- ggplot(metrics, aes(y=weighted_f1, x=tp, fill=name)) + theme_bw() + 
            labs(fill="Probe Set:") + xlab("Time Point") + scale_fill_grey(start = 0.7, end = 0.4)
    }else{
        # if(length(unique(metrics$name))>1){
        #     p <- ggplot(metrics, aes(y=weighted_f1, x=tpTest, fill=name)) 
        # }
        if(length(unique(metrics$nClasses))>1){
            p <- ggplot(metrics, aes(y=weighted_f1, x=tpTrain, fill=as.character(nClasses))) +
                scale_fill_manual(values=wesanderson::wes_palette(n=2, name="GrandBudapest2")[2:1])  + labs(fill="Number of Classes:")
        }
        p <- p + theme_bw() + facet_nested(. ~ ylabel + tpTest) + xlab("Mouse Muscle")
    }
    p <- p + geom_boxplot(outliers=F) + ylab("Weighted F1-score") + theme(legend.position="top")
    return(p)
}

fit_evaluate_glmnet_in_mouseMuscle <- function(nClasses, reweight, probeSubset, lambdaFromMoBl){

    suffix = paste0("_",probeSubset,"Probes_",nClasses,"classes",ifelse(reweight,"_weighted",""),ifelse(lambdaFromMoBl,"_lambdaFromMoBl",""))
    file = paste0("Data/ModelEvaluation/CV_metrics_nSelectedVar_MoMu",suffix,".RData")

    # # If lambdaFromMoBl, get lambda from MoBl lambda optim
    # lambdaSeq = rep(NA, length(TPs))
    # names(lambdaSeq) = TPs
    # if(lambdaFromMoBl){
    #     metricsLatentVar = read.table(paste0("Data/ModelEvaluation/LatentVariable_MoMuOnMoBl","_",nClasses,"classes",ifelse(reweight,"_weighted",""),"_lambdaOptim.txt"), h=T)
    #     tmp = aggregate(lambda~tpTrain, FUN=median, metricsLatentVar)
    #     lambdaSeq[tmp$tpTrain] = tmp$lambda
    # }

    # # Get mouse muscle dataname(s) on which models will be fitted and evaluated 
    # dataName <- grep("mouseMuscle",names(data),value=T) %>%
    #     grep(ifelse(correctMethodComparison,"",correctMethod),., value=T)

    # # Run cv.glmnetcr on multiple CV partitions
    # metricsConfMat = mclapply(1:(nCVRepetition*length(dataName)), function(iName){

    #     # Set seed and get data name
    #     set.seed(iName)
    #     name = dataName[1+(iName%%length(dataName))]
    #     metrics = confMat = NULL

    #     # Get probenames
    #     if(probeSubset == "human"){
    #         probenames = names(probeMap)    
    #     }else{probenames = rownames(data[[name]])}

    #     # Get dataset specific time points
    #     TPs = paste0(sort(as.numeric(sub("h","",unique(metadata[[name]]$timePoint)))),"h")

    #     for(tp in TPs[-(1)]){
    #         # Compute glmnetcr in CV and return metrics 
    #         cvglmnetcr = cv.glmnetcr(name, tp, probenames, lambdaSeq[tp])
    #         metrics = rbind(metrics, cvglmnetcr$metrics)
    #         confMat = rbind(confMat,cvglmnetcr$confMat)
    #     }
    #     return(list(metrics=metrics, confMat=confMat))
    # }, mc.cores=nbcores)

    # # Concat metrics and compute element-wise conMat
    # metrics = do.call(rbind, lapply(metricsConfMat, `[[`, "metrics"))
    # uniqTP = unique(metrics$tp)
    # confMatList = lapply(metricsConfMat, `[[`, "confMat")
    # confMatSum = Reduce(`+`, confMatList)
    # confMatAvg = melt(confMatSum / length(confMatList))
    # confMatAvg$tp = factor(rep(rep(uniqTP,each=nClasses),nClasses), levels=uniqTP)

    # # Save metrics and best lambda
    # nSelectedVar = aggregate(nSelectedVar~name+tp,FUN=median, data=metrics)
    # nSelectedVar = nSelectedVar[nSelectedVar$name=="mouseMuscle_D0",]
    # rownames(nSelectedVar) = nSelectedVar$tp
    # save(metrics, nSelectedVar, confMatAvg, file=file)

    # Plot accuracy and kappa scores
    if(all(file.exists(c(sub("human","mouse",file), sub("mouse","human",file))))){
        metricsAll = NULL
        for(probeSubset in c("human","mouse")){
            suffix = paste0("_",probeSubset,"Probes_",nClasses,"classes",ifelse(reweight,"_weighted",""))
            load(file=paste0("Data/ModelEvaluation/CV_metrics_nSelectedVar_MoMu",suffix,".RData"))
            metrics$name = ifelse(probeSubset=="human","Human Orthologous","Full Array")
            metricsAll = rbind(metricsAll, metrics)
        }

        return(plotMetrics(metricsAll)+ggtitle("a"))
        
    }
}


plot_confusion_matrices <- function(probeSubset, nClasses, reweight){

    suffix = paste0("_",probeSubset,"Probes_",nClasses,"classes",ifelse(reweight,"_weighted",""))
    load(file=paste0("Data/ModelEvaluation/CV_metrics_nSelectedVar_MoMu",suffix,".RData"))
    p <- ggplot(confMatAvg,aes(x=Var1, y=Var2)) +
        geom_tile(aes(fill=value)) +
        geom_text(aes(label = round(value,1)), color = "black") +
        scale_fill_gradient(low = "white", high = "steelblue") +
        facet_wrap(~tp, ) +
        labs(x = "Predicted", y = "Actual", fill = "Average Count") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return(p)
}

pathway_enrichment <- function(freqSelection){

    gsea = NULL
    freqSelection$geneName = probe2geneName[rownames(freqSelection)]
    for(tp in TPTest){
        ranked_gene_list = freqSelection[,tp]
        names(ranked_gene_list) = freqSelection$geneName
        ranked_gene_list = ranked_gene_list[!is.na(names(ranked_gene_list))]
        ranked_gene_list = ranked_gene_list[ranked_gene_list>0]
        tmp = GSEA(
            geneList = sort(ranked_gene_list, decreasing = TRUE),
            TERM2GENE = msigdb[, c("gs_name", "gene_symbol")],
            pvalueCutoff = 0.7,
            pAdjustMethod = "none"
        )
        if(!is.null(tmp)){
            gsea = rbind(gsea, data.frame(TP=tp,tmp@result))
        }
    }
    gsea$log10_pvalue = -log10(gsea$pvalue)
    p <- ggplot(gsea, aes(y= Description, x=log10_pvalue, size=setSize)) + geom_point(aes(col=log10_pvalue)) + labs(size="Gene Set Size", col="-log10(p-value)") +
        ylab("") + xlab("-log10(p-value)") + scale_color_gradient(low = "blue", high ="red") + facet_grid(.~TP) + scale_size(range = c(3, 5)) + theme_bw()
    return(p)
}

Supplementary_Figure_3 <- function(probeSubset, nClasses, reweight, freqSelection, degfilename){

    p_confusion_matrices <- plot_confusion_matrices(probeSubset, nClasses, reweight) + ggtitle("a")
    p_logFC_vs_freqSelection <- plot_logFC_vs_freqSelection(freqSelection, degfilename) + ggtitle("b")
    p_pathway_enrichment <- pathway_enrichment(freqSelection) + ggtitle("c")

    combined_plot <-  patchwork::free(p_confusion_matrices) / patchwork::free(p_logFC_vs_freqSelection) / patchwork::free(p_pathway_enrichment) + plot_layout(heights = c(3.5,3,2.5))
    ggsave("Analyses/ModelEvaluation/Supplementary_Figure_3.pdf", combined_plot, height=10, width=8, device = cairo_pdf)
}


# Optimization function: tries thresholds that split latent into 3 predicted classes
accuracy_fn <- function(latentVar, thresholds, classes) {

    # Sort thresholds
    thresholds = sort(thresholds)

    # penalize near-equal thresholds
    if (diff(thresholds) < 1e-3) return(0)
    # Discretize latent scores using thresholds
    preds = cut(latentVar,
                breaks = c(-Inf, thresholds, Inf),
                labels = 1:3,
                right = TRUE)
    
    # Compute accuracy
    mean(as.numeric(preds) == as.numeric(classes))
}


twoClassesPredictionFromLatentVar <- function(latentvar, classes, threshold, returnThres){

    if(!is.numeric(threshold)){
        latentSeq = seq(min(latentvar), max(latentvar), length.out = nLambda)
        kappa = sapply(1:nLambda,function(i){
            classPred = factor(ifelse(latentvar<latentSeq[i], levels[[probeSubset]][[nClasses]][1], levels[[probeSubset]][[nClasses]][2]), levels=levels[[probeSubset]][[nClasses]], ordered=T)
            confusionMatrix(data=classPred, reference=classes)$overall["Kappa"]
        })
        threshold = latentSeq[which.max(kappa)]
    }

    if(returnThres){
        return(threshold)
    }else{
        return(ifelse(latentvar<threshold+sqrt(.Machine$double.eps), levels[[probeSubset]][[nClasses]][1], levels[[probeSubset]][[nClasses]][2]))
    }
}

threeClassesPredictionFromLatentVar <- function(latentvar, classes, thresholds, returnThres){

    if(!is.numeric(thresholds)){
        if(var(latentvar)>sqrt(.Machine$double.eps)){
            # Initial guess using quantiles
            init_thresh = quantile(latentvar, probs = c(1/3, 2/3))
            # Use optim to find best thresholds
            opt = optim(par = init_thresh,
                        fn = function(par) -accuracy_fn(latentvar, par, classes),  # negative accuracy to minimize
                        method = "L-BFGS-B",
                        lower = min(latentvar),
                        upper = max(latentvar))

            # Extract optimized thresholds
            thresholds = sort(opt$par)
        }else{thresholds = rep(NA, 2)}
    }

    if(returnThres){
        return(thresholds)
    }else{
        if(any(is.na(thresholds))){
            return(rep(NA,length(latentvar)))
        }else{
            return(cut(latentvar,
                breaks = c(-Inf, thresholds, Inf),
                labels = levels(classes),
                right = TRUE)
            )
        }
    }
}

classPredictionFromLatentVar <- function(latentVar, classes, threshold, returnThres){

    apply(latentVar, 2, function(latentvar){
        if(length(unique(classes))==3){
            threeClassesPredictionFromLatentVar(latentvar, classes, threshold, returnThres)
        }else{
            twoClassesPredictionFromLatentVar(latentvar, classes, threshold, returnThres)
        }
    })
}


optimizeLambdaByCV <- function(latentvar, classes){

    iFolds = createFolds(classes, k = nFolds, returnTrain=T)
    classPred = matrix(0, nrow=length(classes), ncol=nLambda)
    
    for(i in seq_along(iFolds)){
        # Predict class from latent variable
        threshold = classPredictionFromLatentVar(latentvar[iFolds[[i]],], classes[iFolds[[i]]], threshold=NULL, returnThres=T)
        if(length(unique(classes)) == 2){threshold = t(threshold)}
        classPred[-iFolds[[i]],] = sapply(1:nLambda, function(iLambda){
            classPredictionFromLatentVar(as.matrix(latentvar[-iFolds[[i]],iLambda]), classes[-iFolds[[i]]], threshold=threshold[,iLambda], returnThres=F)
        })
    }
    # Compute performance metrics
    kappa = apply(classPred, 2, function(classpred){
        confusionMatrix(data=factor(classpred, levels=levels[[probeSubset]][[nClasses]]), reference=classes)$overall["Kappa"]
    })
    ibestlambda = which.max(kappa)
    return(list(latentvar = as.matrix(latentvar[,ibestlambda]), iBestLambda = ibestlambda))
}


train_glmnet_in_mouseMuscle_evaluate_in_blood <- function(nClasses, reweight, probeSubset, lambdaOptim, name, suffix){


    load(paste0("Data/ModelEvaluation/CV_metrics_nSelectedVar_MoMu_",probeSubset,"Probes_",nClasses,"classes",ifelse(reweight,"_weighted",""),".RData"))

    # Get probenames
    if(probeSubset == "human"){
        trainProbes = names(probeMap)
        testProbes = probeMap
    }else{trainProbes = testProbes = rownames(data[[name]])}

    # Initialization
    metricsLatentVar = NULL

    for(tpTrain in TPTest){

        # If lambdaOptim : optimize the number of variables, otherwise use the number estimated in mouseMuscle
        if(lambdaOptim){
            nselectedvar = NA
        }else{
            lambda = NA
            nselectedvar = ifelse(tpTrain=="0h" , median(nSelectedVar[,"nSelectedVar"]) , nSelectedVar[tpTrain,"nSelectedVar"])
        }

        # Fit glmnetcr on mouse muscle
        mouseMuscleMod = glmModel.wrap(y = reactoClasses[["mouseMuscle_D0"]][[tpTrain]][[nClasses]],
                            x = t(data[["mouseMuscle_D0"]][trainProbes,sampleTP[["mouseMuscle_D0"]][[tpTrain]][[nClasses]]]),
                            lambdaSeq = NA,
                            nSelectedVar = nselectedvar,
                            returnLambda = F)

        for(tpTest in TPTest){

            # Compute latent variable manually: X %*% beta + intercept, nClasses is set to 3 to have predictions on all samples
            latentVar = sweep(t(data[[name]][testProbes,sampleTP[[name]][[tpTest]][[nClasses]]]) %*% mouseMuscleMod$mod$beta[trainProbes,], 2, mouseMuscleMod$mod$a0, "+")
            latentVarPos = as.matrix(sweep(latentVar, 2, abs(colMins(as.matrix(latentVar))), "+"))

            metricsLatentVar = rbind(metricsLatentVar, do.call("rbind", mclapply(1:ifelse(grepl("mouse",name),nBootstrapRepetition,1), function(iSeed){

                # Set seed
                set.seed(iSeed)

                # Generate bootstrap sample if iSeed > 1
                bootstrap = data.frame(samples=sampleTP[[name]][[tpTest]][[nClasses]], 
                                classes = reactoClasses[[name]][[tpTest]][[nClasses]])
                if(iSeed > 1){bootstrap <- bootstrap %>%  group_by(classes) %>% group_modify(~ slice_sample(.x, n = nrow(.x), replace = TRUE)) %>% ungroup()}
                latentVarBoost = latentVarPos[bootstrap$samples,]

                # Determine best lambda by cross validation
                if(lambdaOptim){
                    res = optimizeLambdaByCV(latentVarBoost, bootstrap$classes)
                    latentVarBoost = res$latentvar
                    latentVarPos = latentVarPos[,res$iBestLambda]
                    lambda = mouseMuscleMod$mod$lambda[res$iBestLambda]
                }
                # Predict class from latentVarBoost and compute accuracy+Kappa
                if(probeSubset == "mouse"){
                    classPred = classPredictionFromLatentVar(latentVarBoost, bootstrap$classes, threshold=NULL, returnThres=F)
                    weighted_f1 = weighted_f1(data=factor(classPred, levels=levels[[probeSubset]][[nClasses]]), reference=bootstrap$classes)$weighted_f1
                }else{weighted_f1 = NA}

                # Do not keep latent variables in bootstrap samples
                if(iSeed>1){latentVarPos=NA}

                return(data.frame(
                    sample = sampleTP[[name]][[tpTest]][[nClasses]],
                    name = name,
                    tpTrain = tpTrain,
                    tpTest = tpTest,
                    vaccine = sub("^([^_]+)_.*","\\1",sampleTP[[name]][[tpTest]][[nClasses]]),
                    weighted_f1 = weighted_f1,
                    lambda = lambda,
                    latent = latentVarPos
                ))
            }, mc.cores = nbcores)))
        }
    }
    return(metricsLatentVar)
}


train_glmnet_in_mouseMuscle_evaluate_in_blood.wrap <- function(nClass, reweight, probeSubset, lambdaOptim){

    # Get data name, suffix and load nSelectedVar for prediction
    name = paste0(probeSubset,"Blood",correctMethod)
    suffix = paste0(ifelse(grepl("mouse",name),"MoBl","HuBl"),"_",nClass,"classes",ifelse(reweight,"_weighted",""),ifelse(lambdaOptim,"_lambdaOptim",""))
    file = paste0("Data/ModelEvaluation/LatentVariable_MoMuOn",suffix,".txt")
    # Compute latentVar (w/o lambdaOptim) and predict classes
    # metricsLatentVar = train_glmnet_in_mouseMuscle_evaluate_in_blood(nClass, reweight, probeSubset, lambdaOptim, name, suffix)

    # Save latent variables
    # write.table(metricsLatentVar, quote=F, row.names=F, file=file)
    metricsLatentVar = NULL
    if((probeSubset == "mouse") && (all(file.exists(c(sub("2","3",file), sub("3","2",file)))))){
        for(nClasses in 2:3){
            metricsLatentVar = rbind(metricsLatentVar, data.frame(nClasses=nClasses, read.table(file=sub("[23]",nClasses,file), h=T)))
        }
        p_metrics <- plotMetrics(metricsLatentVar)+ggtitle("a")
    }
    if(probeSubset == "human"){
        metricsLatentVar = read.table(file=file, h=T)
        p_metrics = NULL
    }

    # Compute pairwise t-test and plot latent variables
    p_ttest_latentVar = compute_ttest_plot_latentVariables(metricsLatentVar %>% filter(nClasses == nClass), probeSubset, suffix)

    return(list(p_pairwiseTtest=p_ttest_latentVar[["p_pairwiseTtest"]], p_latentVar=p_ttest_latentVar[["p_LatentVar"]], p_metrics=p_metrics))
}


compute_ttest_plot_latentVariables <- function(metricsLatentVar, probeSubset, suffix){


    # Format metricsLatentVar
    nClass = ifelse(grepl("2classes", suffix), 2, 3)
    metricsLatentVar$vaccine = str_replace_all(metricsLatentVar$vaccine, abbreviation2Name)
    metricsLatentVar$vaccine = factor(metricsLatentVar$vaccine,levels=names(vaccine2ReactoClasses))
    metricsLatentVar$xlabel = "Mouse Muscle"
    metricsLatentVar$ylabel = paste0(str_to_title(probeSubset)," blood(", nClass, " classes)")

    # Compute pairwise t-test
    t_test_result = NULL
    for(tpTrain in TPTest){
        for(tpTest in TPTest){
            iLines = which((metricsLatentVar$tpTrain==tpTrain)&(metricsLatentVar$tpTest==tpTest))
            t_test_result = rbind(t_test_result, data.frame(tpTrain=tpTrain, tpTest=tpTest,
                melt(pairwise.t.test(metricsLatentVar[iLines,"latent"], metricsLatentVar[iLines,"vaccine"], na.rm = TRUE, p.adjust.method = "BH")$p.value)
            ))
        }
    }

    # Format t_test_result
    t_test_result = t_test_result[!is.na(t_test_result$value),]
    t_test_result$discreteValue = cut(t_test_result$value, breaks=pvalueBreaks[[1]])
    t_test_result$xlabel = "Mouse Muscle"
    t_test_result$ylabel = paste(str_to_title(probeSubset),"blood")
    t_test_result$ylabelWithClass = paste0(str_to_title(probeSubset)," blood (", nClass, " classes)")
    t_test_fill = t_test_result[1:length(palette[["pvalue1"]]),]
    t_test_fill$discreteValue = levels(t_test_result$discreteValue)
    t_test_fill[,"value"] = NA
    
    # Plot pairwise t-test p-values
	p_pairwiseTtest <- ggplot(rbind(t_test_fill,t_test_result), aes(Var1, Var2, fill=discreteValue)) +
	xlab("") + ylab("") +  labs(fill = "p-value") + scale_fill_manual(values=palette[["pvalue1"]], drop=F) + geom_tile(na.rm = T) + 
    facet_nested(xlabel + tpTrain ~ ylabelWithClass + tpTest) + labs(fill = "p-value") + theme_bw()  + ggtitle(ifelse(probeSubset=="mouse","c","d"))
    if(probeSubset=="mouse"){
        p_pairwiseTtest <- p_pairwiseTtest + theme(axis.text.x = element_text(angle=45, hjust=1, size=7), axis.text.y = element_text(size=7), legend.position="none")
    }else{p_pairwiseTtest <- p_pairwiseTtest + theme(axis.text.x = element_text(angle=45, hjust=1, size=7), axis.text.y = element_text(size=7))}
    
    # Plot latent variable / vaccine / TP
    p_LatentVar <- ggplot(metricsLatentVar, aes(y=latent, x=vaccine, fill=vaccine)) + 
        geom_boxplot(outliers=F) + theme_bw() +
        scale_fill_manual(values=palette[["vaccine"]]) +
        facet_nested(xlabel + tpTrain ~ ylabel + tpTest, scales = "free_y",  independent = "y") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7), axis.text.y = element_text(size=7), legend.position="none") + 
        xlab("") + ylab("Latent variable") + ggtitle(ifelse(probeSubset=="mouse","a","b"))
    if(probeSubset=="human"){p_LatentVar <- p_LatentVar + scale_y_continuous(labels = scales::label_scientific(digits = 2))}

    return(list(p_pairwiseTtest = p_pairwiseTtest, p_LatentVar = p_LatentVar))
}


boxplot_accuracy_across_models <- function(){

    Accuracy = NULL
    for(probeSubset in c("mouse","human")){
        for(nClasses in 2:3){
            for(reweight in c(F,T)){
                for(lambdaOptim in c(F,T)){
                    name = paste0(probeSubset,"Blood",correctMethod)
                    suffix = paste0(ifelse(grepl("mouse",name),"MoBl","HuBl"),"_",nClasses,"classes",ifelse(reweight,"_weighted",""),ifelse(lambdaOptim,"_lambdaOptim",""))
                    metricsLatentVar = read.table(paste0("Data/ModelEvaluation/LatentVariable_MoMuOn",suffix,".txt"), h=T)
                    Accuracy = rbind(Accuracy, 
                        data.frame(species=probeSubset, nClasses=nClasses, reweight=reweight, lambdaOptim=lambdaOptim,
                                    aggregate(Accuracy~tpTrain+tpTest,FUN=mean, data=metricsLatentVar),
                                    aggregate(Kappa~tpTrain+tpTest,FUN=mean, data=metricsLatentVar)
                        )
                    )
                }
            }
        }
    }

    Accuracy$lambdaOptim = ifelse(Accuracy$lambdaOptim, "Lambda Optimization","#Variables estimated in MoMu")
    Accuracy$reweight = ifelse(Accuracy$reweight, "Weighted Observation","Equal weights")
    Accuracy$nClasses = ifelse(Accuracy$nClasses==2, "nClasses=2","nClasses=3")

    pdf("Analyses/ModelEvaluation/Boxplot_accuracy_across_models.pdf", width=10, height=5)
    for(metric in c("Accuracy","Kappa")){
        for(species in c("mouse","human")){
            print(ggplot(Accuracy[Accuracy$species==species,], aes(y=get(metric),x=tpTrain, fill=lambdaOptim)) + 
            geom_boxplot() + facet_wrap(~reweight+nClasses, nrow=1) + theme_bw() + ggtitle(species) + labs(y=metric, fill=NULL))
        }
    }
    dev.off()
}


run_stability_selection_in_mouseMuscle <- function(nClasses, reweight, probeSubset){

    suffix = paste0("_",probeSubset,"Probes_",nClasses,"classes",ifelse(reweight,"_weighted",""))

    # # Get nSelectedVar
    # load(file=paste0("Data/ModelEvaluation/CV_metrics_nSelectedVar_MoMu",suffix,".RData"))
    
    # # Get probenames
    # if(probeSubset == "human"){
    #     probenames = names(probeMap)    
    # }else{probenames = rownames(data[["mouseMuscle_D0"]])}
    
    # # Initialization
    # freqSelection = matrix(0,nrow=length(probenames), ncol=length(TPTest))
    # rownames(freqSelection) = probenames
    # colnames(freqSelection) = TPTest

    # for(tp in TPTest){

    #     # Get lambda sequence
    #     lambdaSeq = glmModel.wrap(y = reactoClasses[["mouseMuscle_D0"]][[tp]][[nClasses]],
    #                         x = t(data[["mouseMuscle_D0"]][probenames,sampleTP[["mouseMuscle_D0"]][[tp]][[nClasses]]]),
    #                         lambdaSeq = NA,
    #                         nSelectedVar = NA,
    #                         returnLambda = T)

    #     # Compute selection frequency
    #     freqselection = compute_phat_glmnet(y = reactoClasses[["mouseMuscle_D0"]][[tp]][[nClasses]],
    #                         x = t(data[["mouseMuscle_D0"]][probenames,sampleTP[["mouseMuscle_D0"]][[tp]][[nClasses]]]),
    #                         B = 1000,
    #                         nSelectedVar = nSelectedVar[tp,"nSelectedVar"],
    #                         lambdaSeq = lambdaSeq)

    #     freqSelection[rownames(freqSelection),tp] = freqselection[,1]
    # }
    # write.table(freqSelection, quote=F, sep="\t", file=paste0("Data/ModelEvaluation/Variable_selection_frequency_MoMu",suffix,".txt"))
    freqSelection = read.table(file=paste0("Data/ModelEvaluation/Variable_selection_frequency_MoMu",suffix,".txt"), h=T, stringsAsFactors=F, sep="\t", check.names = FALSE)

    iGenes2show = unique(as.numeric(apply(freqSelection, 2, function(x) order(x,decreasing=T)[1:nTopGeneStabs])))
    topGenes = as.data.frame(freqSelection[iGenes2show,])
    
    topGenes$geneName = sapply(rownames(topGenes), function(x){
        y = probe2geneName[x]
        if(is.na(y)){x}else{y}
    })

    p_freqSelec = plot_variable_selection_frequency(topGenes, suffix)
    p_trajectories = plot_trajectories(topGenes, nClasses, name="mouseMuscle_D0", title="c")
    return(list(p_freqSelec=p_freqSelec, topGenes=topGenes, p_trajectories=p_trajectories, freqSelection=freqSelection))
}


plot_variable_selection_frequency <- function(topGenes, suffix){

    topGenes.df = melt(topGenes)
    colnames(topGenes.df) = c("Name","TimePoint","Frequency")

    # Step 1: Reorder Name factor by total Frequency
    topGenes.df <- topGenes.df %>%
    group_by(Name) %>%
    mutate(TotalFreq = sum(Frequency)) %>%
    ungroup() %>%
    mutate(Name = reorder(Name, TotalFreq))

    p <- ggplot(topGenes.df, aes(x = Name)) +
        geom_bar(aes(y = Frequency, fill = TimePoint), position = position_stack(reverse = T), stat="identity", width = .7) + 
        coord_flip() + theme_bw() + theme(legend.position="top") + xlab("") + ylab("Cumulative Frequency") + ggtitle("b") + labs(fill="Time Point:")
    return(p)
}


plot_trajectories <- function(topGenes, nClasses, name, title){

    # Offset expression values so they are positives 
    offset = abs(min(data[[name]][rownames(topGenes),unique(unlist(sampleTP[[name]]))]))

    # Combine time points and melt
    df = do.call("rbind",lapply(TPs, function(tp){
        df = melt(offset + data[[name]][rownames(topGenes),sampleTP[[name]][[tp]][[nClasses]]])
        df$tp = sub("^.*_(MU|BL)_([0-9]+)h_[1-5]","\\2",df$Var2)
        df$class = vaccine2ReactoClasses[sub("^(.*)_(MU|BL)_.*","\\1",df$Var2)]
        return(df)
    }))
    df$tp = factor(df$tp, levels=sub("h","",TPs))
    df$geneName = topGenes[df$Var1,"geneName"]
    df$class = factor(df$class, c("high","medium","low"))

    # Remove outliers
    df_no_outliers <- df %>%
    group_by(geneName, class) %>%
    filter(
        value > quantile(value, 0.25) - 1.5 * IQR(value),
        value < quantile(value, 0.75) + 1.5 * IQR(value),
        !(tp %in% c("-24","0"))
    ) %>% ungroup()

    # Plot trajectories
    p <- ggplot(df_no_outliers, aes(x = tp, y = value, color = class, group = class)) + 
        geom_jitter(alpha=0.25, width = 0.3, height = 0) + theme_bw() + 
        geom_smooth(method = "loess") + facet_wrap(~geneName, scales="free", ncol=6) +
        scale_color_manual(values = palette[["reactoClasses"]]) +
        theme(axis.text.x = element_text(angle=45, hjust = 1, size=7), axis.text.y = element_text(size=7), panel.spacing.x=unit(0.1, "lines"),
            legend.background=element_blank(), legend.key=element_blank(), legend.position="top", strip.text = element_text(size=8)) +
        guides(color=guide_legend(override.aes=list(fill=NA))) +
        labs(x = "Time Point (hours)", y = "log2(Gene Expression)", color = "Reactogenicity Class:", fill = "") + ggtitle(title)
    return(p)
}

figure2 <- function(p_metrics, stabilityRes){

    layout <- "
    ACC
    ACC
    BCC
    BCC
    BCC
    "
    combined_plot <-  patchwork::free(p_metrics) + patchwork::free(stabilityRes[["p_freqSelec"]]) + patchwork::free(stabilityRes[["p_trajectories"]]) + plot_layout(design = layout)
    ggsave("Analyses/ModelEvaluation/Figure_2.pdf", combined_plot, height=8, width=11, device = cairo_pdf)
}

plot_logFC_vs_freqSelection <- function(freqSelection, degFilename){
    
    dfAll = NULL
    DiffAn = read.xlsx(degFilename)
    
    for(tp in unique(colnames(freqSelection))){

        selectedGenes = rownames(freqSelection[freqSelection[,tp]>0,])
        df = data.frame(tp=tp,  
            DiffAn[(DiffAn$timepoint==tp)&(DiffAn$geneID%in%selectedGenes),c("adj.P.Val","logFC","reactoClass","geneID")]
        )
        df$freqSelection = freqSelection[df$geneID,tp]
        dfAll = rbind(dfAll, df)
    }

    dfAll$logAdjPVal = -log10(dfAll$adj.P.Val)
    dfAll$absLogFC = abs(dfAll$logFC)
    dfAll$discLogAdjPVal = cut(dfAll$adj.P.Val, breaks=pvalueBreaks[[2]])

    p <- ggplot(dfAll, aes(x=freqSelection, y=absLogFC, col=discLogAdjPVal)) +
        scale_color_manual(values=palette[["pvalue2"]]) +
        geom_point() +
        scale_y_log10() + 
        theme_bw() +
        ylab("|log2(Fold-Change)|") + xlab("Selection Frequency") + labs(col="p-value") + 
        facet_wrap(.~tp,scales="free")
    return(p)
}


# Old_Figure_3 <- function(pairwiseTtest_LatentVar, topGenes){

#     trajectories <- plot_trajectories(topGenes, nClasses=nClasses, name="mouseBlood_D0", title="c")

#     layout <- "
#     AACCC
#     BBCCC
#     "
#     combined_plot <-  patchwork::free(pairwiseTtest_LatentVar[["mouse"]][["p_latentVar"]]) + patchwork::free(pairwiseTtest_LatentVar[["human"]][["p_latentVar"]]) +
#         patchwork::free(trajectories) + plot_layout(design = layout)
#     ggsave("Analyses/ModelEvaluation/Figure_3.pdf", combined_plot, height=7, width=11, device = cairo_pdf)
# }

Figure_3 <- function(pairwiseTtest_LatentVar, topGenes){

    trajectories <- plot_trajectories(topGenes, nClasses=nClasses, name="mouseBlood_D0", title="b")
    layout <- '
    #B
    AB
    AB
    AB
    #B
    '
    combined_plot <-  free(pairwiseTtest_LatentVar[["mouse"]][["p_latentVar"]]) + free(trajectories) + 
    plot_layout(design=layout, widths=c(2,3))
    ggsave("Analyses/ModelEvaluation/Figure_3.pdf", combined_plot, height=8, width=11, device = cairo_pdf)
}

scatterPlot_LFC <- function(degFilename, freqSelection){

    freqSelection$geneID = rownames(freqSelection)
    freqSelection.df <- reshape2::melt(freqSelection,variable.name="timepoint",value.name="Frequency") %>%
        filter(Frequency > 0.1)

    DiffAn = list()
    for(degfilename in degFilename){
        DiffAn[[sub(".*(blood|muscle).*","\\1",degfilename)]] = read.xlsx(degfilename)
    }
    df <- full_join(DiffAn[["muscle"]],DiffAn[["blood"]], by=c("timepoint","geneID"), suffix = c(".mu", ".bl")) %>%
        filter(timepoint %in% TPTest) %>%
        right_join(freqSelection.df, by=c("timepoint","geneID"))

    p_scatterPlot_LFC <- ggplot(df, aes(x=logFC.mu, y=logFC.bl)) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method="pearson") +
        labs(col="Frequency") +
        xlab("log2(Fold-Change) in Mouse Muscle") + ylab("log2(Fold-Change) in Mouse Blood") +
        theme_bw() + facet_wrap(.~timepoint, scales="free") + ggtitle("b")

    return(p_scatterPlot_LFC)
}

Supplementary_Figure_4 <- function(degFilename, freqSelection, pairwiseTtest_LatentVar){


    p_scatterPlot_LFC <- scatterPlot_LFC(degFilename, freqSelection)

    layout <- "
    AAAAAAAAA
    BBBBBBBBB
    CCCCDDDDD
    "
    figd <- pairwiseTtest_LatentVar[["3"]][["mouse"]][["p_pairwiseTtest"]] +
        ggtitle("d") +
        theme(legend.position = "right")


    combined_plot <-  patchwork::free(pairwiseTtest_LatentVar[["2"]][["mouse"]][["p_metrics"]]) + patchwork::free(p_scatterPlot_LFC) + 
        patchwork::free(pairwiseTtest_LatentVar[["2"]][["mouse"]][["p_pairwiseTtest"]]) + patchwork::free(figd) +
        plot_layout(design = layout)
    ggsave("Analyses/ModelEvaluation/Supplementary_Figure_4.pdf", combined_plot, width = 10, height = 10, device = cairo_pdf)

}
