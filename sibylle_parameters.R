#############################
### SET GLOBAL PARAMETERS ###
#############################


# LOAD MINIMAL LIBRARIES FOR PARAMETER SETTING 
library(RColorBrewer)
library(stringr)

# SET SEED
set.seed(123456)

# INITIALIZATION
map = data = metadata = levels = palette = pvalueBreaks = list()

### NUMERIC/BOOLEAN/CHARACTER PARAMETERS
correctMethodComparison = F
nbcores = 80
nFolds = 5
nTopVarGenes = 1000
nCVRepetition = 100
nBootstrapRepetition = 100
nLambda = 40
maxit = 10000
nTopGeneStabs = 10
adjPValThres = 0.05
LFCThres = 1
bgrdNoiseThres = 5.5 # Median log2(expression) value in mouse muscle and blood 
correctMethod = "_D0"
degFilename = c(mouseMuscle = "Analyses/biovacsafe_mouse/DEG/biovacsafe_mouse_differential_analysis_reactoclass_muscle_results_all_genes.xlsx",
                mouseBlood = "Analyses/biovacsafe_mouse/DEG/biovacsafe_mouse_differential_analysis_reactoclass_blood_results_all_genes.xlsx")

### PVALUE DISCRETIZATION
pvalueBreaks[[1]] = c(0,1e-3,1e-2,0.05,1e-1,1)
pvalueBreaks[[2]] = c(0,1e-10,1e-5,1e-2,1)

### REACTOGENICITY CLASSES, TPs, OUTLIERS
levels[["mouse"]][2:3] = list(c("low","high"), c("low","medium","high"))
levels[["human"]][2:3] = list(c("low","medium"), c("low","medium"))
TPTest = c("24h","48h","72h")
TPs = c("-24h","0h","4h","8h","24h","48h","72h","96h","120h","168h")
outliers = c("TriFLU_BL_168h_CRC305C9072343","SAL_BL_120h_CRC305B8042210","YFV_BL_168h_CRC305A7069148","YFV_BL_168h_CRC305A7042134","VZV_BL_.24h_CRC305A7039128","YFV_BL_96h_CRC305A7016107","VZV_BL_.24h_CRC305A7073144","SAL_BL_72h_CRC305A7038126","VZV_BL_168h_CRC305A7072127","VZV_BL_72h_CRC305A7048125" )

### DSname, abbreviation2Name, vaccine2ReactoClasses
DSname = list(
    name = c("mouseMuscle","mouseBlood","humanBlood"),
    species = c("mouse","mouse","human"),
    compartment = c("muscle","blood","blood")
)

abbreviation2Name = c(
    SAL="Saline",
    Placebo="Saline",
    Engerix="Engerix_B",
    ENG="Engerix_B",
    Pentavac="Pentavac_SD",
    PERT="Pentavac_SD",
    TriFLU.MF59="Fluad",
    TriFLU="Agrippal",
    VZV="Varilrix",
    Varilix="Varilrix",
    YFV="Stamaril"
)
# names(abbreviation2Name) = paste0(names(abbreviation2Name),"$")
vaccine2ReactoClasses = c(SAL="low", TriFLU="low", ENG="low", PolyIC="low", VZV="low", YFV="low", TriFLU.MF59="medium", IFA="medium", LPS="high", PERT="high")
names(vaccine2ReactoClasses) = str_replace_all(names(vaccine2ReactoClasses), abbreviation2Name)

### DEFINE COLOR PALETTE
palette[["vaccine"]] = c("#fc8d62","brown4","#00ABFF","#FFA500","cornflowerblue","navajowhite2","purple4","magenta","#a6d854","darkgrey")
names(palette[["vaccine"]]) = c("Engerix_B","IFA","LPS","Pentavac_SD","Stamaril","Varilrix","PolyIC","Agrippal","Fluad","Saline")
palette[["vaccine"]] = palette[["vaccine"]][names(vaccine2ReactoClasses)]
palette[["vaccine"]] = palette[["vaccine"]][length(palette[["vaccine"]]):1]
palette[["time"]] = brewer.pal(length(TPs)-2, "Set1")
names(palette[["time"]]) = TPs[3:(length(TPs))]
palette[["all"]] = unlist(palette[c("vaccine","time")])
names(palette[["all"]]) = sub("vaccine\\.|time\\.","",names(palette[["all"]]))
palette[["pvalue1"]] = brewer.pal(length(pvalueBreaks[[1]])-1, "RdYlBu")
palette[["pvalue2"]] = brewer.pal(length(pvalueBreaks[[2]])-1, "RdYlBu")
palette[["reactoClasses"]] = brewer.pal(3, "YlOrRd")

for(i in 1:2){
    palette[[paste0("pvalue",i)]] = palette[[paste0("pvalue",i)]][length(palette[[paste0("pvalue",i)]]):1]
}

