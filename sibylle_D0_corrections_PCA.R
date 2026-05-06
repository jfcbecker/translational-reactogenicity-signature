######################
### SVA CORRECTION ###
######################


data[paste0(DSname$name,"_sva")] = mclapply(DSname$name, function(name){
    mod <- model.matrix(~ vaccineTP, data = metadata[[name]])
    mod0 <- model.matrix(~ 1, data = metadata[[name]])

    n.sv = num.sv(data[[name]],mod,method="leek")
    svobj = sva(data[[name]],mod,mod0,n.sv=max(n.sv,2))
    return(removeBatchEffect(data[[name]], covariate=svobj$sv))
})

metadata[paste0(DSname$name,"_sva")] = metadata[DSname$name]


###############################################
### MERGE AND RE-NORMALIZE THE TWO DATASETS ###
###############################################


### MERGE THE TWO DATASETS AND METADATA INTO A SINGLE MATRIX
sharedProbes = intersect(rownames(data[["mouseMuscle_D0"]]),rownames(data[["mouseBlood_D0"]]))
data[["mouse_D0"]] = cbind(data[["mouseMuscle_D0"]][sharedProbes,],data[["mouseBlood_D0"]][sharedProbes,])

### RE-NORMALIZE DATA ACROSS COMPARTMENTS
data[["mouse_Combined"]] = normalize.quantiles(as.matrix(data[["mouse_D0"]]),copy=TRUE)
rownames(data[["mouse_Combined"]]) = rownames(data[["mouse_D0"]])
colnames(data[["mouse_Combined"]]) = colnames(data[["mouse_D0"]])
metadata[["mouse_Combined"]] = as.data.frame(do.call(rbind, metadata[1:2]))

### REMOVE mouse_D0 FROM DATA
data = data[names(data)!="mouse_D0"]


#########################
### PCA VISUALIZATION ###
#########################


percentExplained = dfAll = NULL

for(name in names(data)){
    igeneMostVar = order(rowSds(data[[name]]),decreasing=T)[1:min(c(nTopVarGenes,nrow(data[[name]])))]
    res = prcomp(data[[name]][igeneMostVar,])
    percentExplained = rbind(percentExplained, round(100*res$sdev[1:10]/sum(res$sdev)))
    df = data.frame(res$rotation[,1:2],metadata[[name]],dataset=name)
    dfAll = rbind(dfAll,df)
}

dfAll = dfAll[!dfAll$timePoint%in%c("-24h","0h"),]
dfAll$dataset = factor(str_to_title(sub("([MB])","_\\1",sub("_D0","",dfAll$dataset))), levels=c("Mouse_muscle","Mouse_blood","Mouse_combined","Human_blood"))
dfAll$timePoint = factor(dfAll$timePoint, TPs)
dfAll$Class = factor(vaccine2ReactoClasses[dfAll$vaccine], levels=c("low","medium","high"))
dfAll$Vaccine = factor(dfAll$vaccine, levels=names(palette[["vaccine"]]))


p <- ggplot(dfAll,aes(x=PC1, y=PC2)) + 
    geom_point(aes(col = Vaccine, shape = compartment), size=1) + 
    scale_color_manual(values = palette[["vaccine"]]) + 
    ggnewscale::new_scale_color() +
    stat_ellipse(aes(col=Class)) + 
    scale_color_manual(values = c("low" = "black", "medium" = "darkgray", "high" = "lightgray")) +
    theme_bw() + 
    theme(axis.text.x= element_blank(), 
    axis.text.y= element_blank(), 
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom", 
    legend.box = "vertical",        
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing = unit(-5, "pt"),
    legend.spacing.y = unit(10, "pt")) +
    guides(col = guide_legend(title.position="left")) +
    ggh4x::facet_grid2(timePoint~dataset, scales = "free", independent = "all")
ggsave("Analyses/correctionMethods/Supplementary_Figure_2.pdf", width=5.5, height=8)


percentExplained.df = melt(percentExplained)
colnames(percentExplained.df) = c("dataset","PC","percentage")
percentExplained.df[,"dataset"] = factor(str_to_title(sub("([MB])","_\\1",sub("_D0","",names(data)[percentExplained.df[,"dataset"]]))), levels=c("Mouse_muscle","Mouse_blood","Mouse_combined","Human_blood"))


p = list()
for(var in c("Vaccine","timePoint")){

    p[[var]] <- ggplot(dfAll,aes(x=PC1, y=PC2)) + 
    geom_point(aes(col = .data[[var]], shape = compartment)) + 
    stat_ellipse(aes(col=.data[[var]])) + 
    scale_color_manual(values = palette[["all"]]) + 
    theme_bw() + 
    ggh4x::facet_grid2(.~dataset, scales = "free", independent = "all")
    if(var=="Vaccine"){p[[var]] <- p[[var]] +  guides(shape = "none") + xlab("") + ggtitle("b")}

}
p[["varExplained"]] <- ggplot(percentExplained.df, aes(y=percentage,x=PC)) + 
    geom_bar(stat="identity", col="gray") + ggtitle("a") + 
    ggh4x::facet_grid2(.~dataset) + theme_bw() + 
    xlab("Principal Component") + ylab("Variance Explained (%)")


combined_plot <-  p[["varExplained"]] / p[["Vaccine"]] / p[["timePoint"]] + plot_layout(axis_titles = "collect", guides = "collect")
ggsave("Analyses/correctionMethods/Supplementary_Figure_1.pdf",combined_plot, width=12, height=8)

