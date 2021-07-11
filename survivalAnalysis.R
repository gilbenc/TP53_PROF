#' ---
#' title: "TP53 Survival Analysis"
#' output: -
#' ---

library(RVAideMemoire)
library(survminer)
library(survival)
library(caret)
library(pROC)
library(dplyr)
library(ggplot2)
library(tidyr)
library(glue)
#library(ModelMetrics)
library(OpenMPController) # for Kaggle backend
library(readr)
library(vtreat)
library(xgboost)
library(ComplexHeatmap)
library(reshape)
library(ggfortify)
library(ramify)
library(e1071)
library(Matrix)
library(DiagrammeR)
library(mboost)
library(randomForestSRC)
library("RColorBrewer")

###############
## Load Data ##
###############

# path to survival data
survival.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\survival_data\\"

### Clinical data ###
cancer_type <- "TCGA_full"
clinical.path = file.path(survival.path, paste0(cancer_type, "_clinical_data.csv"), fsep = .Platform$file.sep)
clinical.dat<-read.table(clinical.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)


### mutation data ###
mutation.path = file.path(survival.path, paste0(cancer_type, "_mutations_data.csv"), fsep = .Platform$file.sep)
mutation.dat<-read.table(mutation.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)


# change to comp/all_features
# input.path2<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\GBM- Test 2 models (all features, comp only)//"
# all_data.path = file.path(input.path2, "all.mutations.predictions_with_features.test.all_features.csv", fsep = .Platform$file.sep)

### all_data ###
input.path2<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\GBM- Test final model\\"
all_data.path = file.path(input.path2, "all.mutations.predictions_with_features.test.functional.csv", fsep = .Platform$file.sep)
all_data<-read.table(all_data.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
# fix row names
row.names(all_data) <- all_data[,1]
all_data <- all_data[,-1]
# combine predictions and labels into one column. predictions > label.
all_data$predictions[all_data$predictions == "null"] <- all_data$label[all_data$predictions =="null"]
# or should label > predictions?
# all_data$label[all_data$label == "null"] <- all_data$predictions[all_data$label =="null"]


### orig.dat ###
input.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53"
alldata.path = file.path(input.path, "all_data.csv", fsep = .Platform$file.sep)
orig.dat<-read.table(alldata.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F, na.strings = c("", "-", "No_Data", "No data (synonymous)", "#N/A", "No_Data", " No_Data", "No_Data (synonymous)", "*", "?"),
                     row.names='cDNA_Variant')
# include only missense data
orig.dat<-orig.dat[orig.dat$Variant_Classification=='Missense',]
# only unique protein change values.
orig.dat <- orig.dat[!duplicated(orig.dat$Protein_p1_TP53),]

###############
## Functions ##
###############

#Decides which features to use (comp, act, vivo)
# and whether to include labels in data (includeLabel).
choose.features<-function(dat, comp, act, vivo, includeLabel = F) {
        

        # do we want activity features? (unique for T53, therefore better without, so model can  
        # indicate for other genes as well.)
        if (comp==T) {
                if(act == T) {
                        if(vivo == T)
                                cols <- features
                        else
                                cols <- features.no.vivo
                }
                else {
                        if(vivo == T)
                                cols <- features.no.act
                        else
                                cols <- features.comp
                }
        } 
        else {
                if(act == T) {
                        if(vivo == T)
                                cols <- features.no.comp
                        else
                                cols <- features.act
                }
                else {
                        if(vivo == T)
                                cols <- features.vivo
                        else
                                print("no features selected. redefine comp, act, vivo.")
                }
        }
        
        # do we want labels to be included or concluded from data? this depends on the learning  
        #algorithm's properties (some require labels seperately)
        if (includeLabel==T) {
                cols<-c(cols, "Overall Survival (Months)", "Overall Survival Status")
        }
        dat <- dat[,cols]
        
        return(dat)
}

# add the column with changed protein (from orig.dat) to the features dataset (all_data) in column protein_clean
add_protein_clean_to_features <- function(all_data, orig_dat) {
        all_data$mut_name <- row.names(all_data)
        orig.dat$mut_name <- row.names(orig.dat)
        mut_name <- orig.dat$mut_name
        
        #protein change data.
        orig.dat$protein_clean <- substr(orig.dat$Protein_p1_TP53, 3, 7)
        protein_clean <- orig.dat$protein_clean
        
        #all mutation names and protein changes from orig.dat
        mut_prot <- cbind(mut_name, protein_clean)
        
        #add protein data to all_data.
        all_data <- merge(x = all_data, y=mut_prot, by = "mut_name")
        all_data <- all_data[!duplicated(x=all_data$protein_clean),]
        return(all_data)
}


# input: in one matrix all features (dataset4 without label), than the matching months and status for             # those mutations. 
survival_ML <- function(data_survival) {
        set.seed(123456)
        #divide dataset. keep same amounts of y=dataset$labels on every division
        index <- createDataPartition(y=data_survival$label, p=0.8, list=FALSE) 
        trainSet <- data_survival[ index,]
        testSet <- data_survival[-index,]
        
        x <- as.matrix(choose.features(dat = trainSet, comp = T, vivo = T, act = T, includeLabel = F))
        y_months <- trainSet$`Overall Survival (Months)`
        y_status <- as.integer(trainSet$`Overall Survival Status` == "DECEASED")
        data <- data.frame(x, y_months, y_status)
        model_rfsrc <- rfsrc(Surv(y_months, y_status) ~ ., data = data, ntree = 500, nodesize = 5, na.action = "na.impute")
        print(model_rfsrc)
        
        testSet <- choose.features(dat = testSet, comp = T, vivo = T, act = T, includeLabel = T)
        x <- as.matrix(testSet[,1:42])
        y_months <- testSet$`Overall Survival (Months)`
        y_status <- testSet$`Overall Survival Status`
        y_status <- as.integer(y_status == "DECEASED")
        data_test <- data.frame(x, y_months, y_status)
        
        #create predictions
        pred <- predict.rfsrc(model_rfsrc, data_test)
        print(pred)
        return (model_rfsrc)
}




# print coxph values for comparison of D-ND, D-NotP53, ND-NotP53.
compare_coxph <- function(dataframe) {
        compare_D_ND <- dataframe[dataframe$label != "Not_P53",]
        compare_NotP53_D <- dataframe[dataframe$label != "ND",]
        compare_ND_NotP53 <- dataframe[dataframe$label != "D",]
        # head(compare_ND_NotP53$label)
        print("compare D to ND")
        use_coxph(compare_D_ND)
        print("compare Not P53 to ND")
        use_coxph(compare_ND_NotP53)
        print("compare D to Not P53")
        use_coxph(compare_NotP53_D)
}
# print coxph values for dataframe.
use_coxph <- function(dataframe_compare) {
        dataframe_compare$label <- as.factor(as.character(dataframe_compare$label))
        fit_compare <- coxph(Surv(dataframe_compare$months, dataframe_compare$status)~dataframe_compare$label, data = dataframe_compare)
        print(cox.zph(fit_compare))
        cat('\n')
}


############
### MAIN ###
############


input.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53"
analysis.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis"


#read feature selection file and create subcategories of features.
FEATURES_DF_PATH = file.path(input.path, "data_list.csv", fsep = .Platform$file.sep)
df_features = read.table(FEATURES_DF_PATH, sep=",", header=T, stringsAsFactors = F, check.names = F)
features = df_features[grepl("Yes",df_features$`as a feature`),]$'Field Name'
features.act<-df_features["Yes (activity)" == df_features$`as a feature`,]$'Field Name'
features.vivo <- df_features["Yes (vivo)" == df_features$'as a feature',]$'Field Name'
features.no.act <-setdiff(features, features.act)
features.no.vivo <- setdiff(features, features.vivo)
features.comp <- setdiff(features.no.act, features.vivo)
features.no.comp <- c(features.act, features.vivo)

### prepare data for survival analysis ###
# cancer_type <- "TCGA_full"
print(paste0("survival analysis for ", cancer_type))

### Clinical data ###
# Remove Duplicates (samples that appear more than once- 
# mostly same sample in 2 different tumor defenitions - glioma+glioblastoma)
clinical.dat <- clinical.dat[!(duplicated(clinical.dat$`Sample ID`)),]
# change column names for survival features.
colnames(clinical.dat)[grepl(pattern = "Months", x = colnames(clinical.dat))] <- "OS_MONTHS" 
colnames(clinical.dat)[grepl(pattern = "Status", x = colnames(clinical.dat))] <- "OS_STATUS"          
# remove NAs in survival data (months/status) from clinical table 
clinical.dat <- clinical.dat[!is.na(clinical.dat$OS_MONTHS) & !is.na(clinical.dat$OS_STATUS),]

### mutation data ###
# head(mutation.dat)
# Remove Duplicates in mutations: if same Sample ID and same Protein Change- remove.
check.mut.dups <- data.frame(mutation.dat$`Sample ID`, mutation.dat$`Protein Change`)
mutation.dat <- mutation.dat[!duplicated(check.mut.dups),]
rm(check.mut.dups)

### Data cleaning:
# remove from entire analysis (both mutation and clinical):
# samples that appear more than once in mutation.dat (1 sample, more than 1 mutation).
# samples are removed COMPLETELY, not just the duplicates.
mutation.dup <- mutation.dat[duplicated(x = mutation.dat$`Sample ID`),]
mutation.dat <- mutation.dat[!(mutation.dat$`Sample ID` %in% mutation.dup$`Sample ID`),]
clinical.dat <- clinical.dat[!(clinical.dat$`Sample ID` %in% mutation.dup$`Sample ID`),]
rm(mutation.dup)
# remove from entire analysis (both mutation and clinical):
# patients that had more than one sample. again remove completely.
patient.dup <- clinical.dat[duplicated(clinical.dat$`Patient ID`),]
mutation.dat <- mutation.dat[!(mutation.dat$`Sample ID` %in% patient.dup$`Sample ID`),]
clinical.dat <- clinical.dat[!(clinical.dat$`Patient ID` %in% patient.dup$`Patient ID`),]
rm(patient.dup)

# create Not P53 clinical dataset. (samples not in mutation.dat)
# must be done before changes are made in mutation.dat, that could cause removal of 
# protein change values.
not_p53.dat <- clinical.dat[!(clinical.dat$`Sample ID` %in% mutation.dat$`Sample ID`),]
not_p53.dat$label <- "Not_P53"

# divide missense from non-missense mutations.
mutation.dat.nonMissense <- mutation.dat[mutation.dat$'Mutation Type' != "Missense_Mutation",]
mutation.dat.Missense <- mutation.dat[mutation.dat$`Mutation Type` == "Missense_Mutation",]

# merge survival and mutation (the residue changed) data.
clinical.dat.Missense <- clinical.dat[(clinical.dat$`Sample ID` %in% mutation.dat.Missense$`Sample ID`),]
survival.dat.Missense <- merge(x = clinical.dat.Missense, y = mutation.dat.Missense, by= "Sample ID")
rm(clinical.dat.Missense, mutation.dat.Missense)

# create dataset for Truncating mutations label them 'D'.
## is Truncating really what we think?:
# "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Region", "Splice_Site"
truncate_mutations <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Region", "Splice_Site")
mutation.dat.truncate <- mutation.dat.nonMissense[(mutation.dat.nonMissense$`Mutation Type` %in% truncate_mutations),]
clinical.dat.truncate <- clinical.dat[(clinical.dat$`Sample ID` %in% mutation.dat.truncate$`Sample ID`),]
survival.dat.truncate <- merge(x = clinical.dat.truncate, y = mutation.dat.truncate, by= "Sample ID")
survival.dat.truncate$label <- 'T'
rm(clinical.dat.truncate, truncate_mutations, mutation.dat.truncate, mutation.dat.nonMissense)


### merge survival.dat.missense + all_data (TCGA + UMD)
# change 'protein Change' into 'protein_clean').
ind = which(colnames(survival.dat.Missense)== "Protein Change")
colnames(survival.dat.Missense)[ind] <- "protein_clean"
# add the column with changed protein (from orig.dat) to the features dataset (all_data)
# in column protein_clean.
all_data <- add_protein_clean_to_features(all_data, orig_dat)
# merge features dataset with survival dataset for Missense mutations of p53.
data_survival <- merge(survival.dat.Missense, all_data, by = "protein_clean")
data_survival$cancer_type <- cancer_type
# write.csv(data_survival, file = paste0("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\survival_data\\tables_by_type\\", cancer_type, "_final_table_13_02_2020.csv"))
rm(survival.dat.Missense, clinical.dat, mutation.dat)

### run survival_analysis. ###
months <- data_survival$OS_MONTHS
months_truncate <- survival.dat.truncate$OS_MONTHS
months_notp53 <- not_p53.dat$OS_MONTHS
months <- combine(months, months_truncate, months_notp53)

status <- data_survival$OS_STATUS
status_truncate <- survival.dat.truncate$OS_STATUS
status_notp53 <- not_p53.dat$OS_STATUS
status <- combine(status, status_truncate, status_notp53)
status <- as.integer(status == "DECEASED")

label <- data_survival$predictions
label_truncate <- survival.dat.truncate$label
label_notp53 <- not_p53.dat$label
label <- combine(label, label_truncate, label_notp53)
dataframe <- data.frame(months, status, label)



# compare_coxph(dataframe)

fit <- survfit(Surv(months, status) ~ label, data = dataframe)
fit_coxph <- coxph(Surv(months, status) ~ label, data = dataframe)
compare_D_ND <- dataframe[dataframe$label == "D" | dataframe$label == "ND",]
compare_NotP53_D <- dataframe[dataframe$label == "D" | dataframe$label == "Not_P53",]
compare_ND_NotP53 <- dataframe[dataframe$label == "Not_P53" | dataframe$label == "ND",]
compare_T_NotP53 <- dataframe[dataframe$label == "Not_P53" | dataframe$label == "T",]
compare_ND_T <- dataframe[dataframe$label == "T" | dataframe$label == "ND",]
compare_D_T <- dataframe[dataframe$label == "D" | dataframe$label == "T",]

print(coxph(Surv(compare_D_ND$months, compare_D_ND$status)~compare_D_ND$label))
print(coxph(Surv(compare_NotP53_D$months, compare_NotP53_D$status)~compare_NotP53_D$label))
print(coxph(Surv(compare_ND_NotP53$months, compare_ND_NotP53$status)~compare_ND_NotP53$label))
print(coxph(Surv(compare_ND_T$months, compare_ND_T$status)~compare_ND_T$label))
print(coxph(Surv(compare_D_T$months, compare_D_T$status)~compare_D_T$label))
print(coxph(Surv(compare_T_NotP53$months, compare_T_NotP53$status)~compare_T_NotP53$label))
print(fit_coxph)
print(fit)

# display.brewer.pal(n = 8, name = 'Set1')
# brewer.pal(n = 8, name = "Set1")
# jpeg(filename = paste0("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\survival_data\\plots\\", data_survival$cancer_type[1], "`_notp53.jpg"), width = 533, height = 437, units = "px") +
# ggplot2::autoplot(fit, col = c(2,4), conf.int = FALSE) + ggtitle(data_survival$cancer_type[1])
#dev.off()
# display.brewer.pal(n = 12, name = 'Set3')
# brewer.pal(n = 12, name = "Set3")
myCol <- brewer.pal(4, "Set3")

ggsurv <- ggsurvplot(
        fit,                     # survfit object with calculated statistics.
        size = 5,
        data = dataframe,             # data used to fit survival curves.
        # risk.table = TRUE,       # show risk table.
        # pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,         # show confidence intervals for 
        # point estimates of survival curves.
        palette = c("#80B1D3", "#FFED6F", "#B3DE69", "#BC80BD"),
        xlim = c(0,400),         # present narrower X axis, but not affect
        # survival estimates.
        # xlab = FALSE, #"Time in months",   # customize X axis label.
        # ylab = FALSE,
        break.time.by = 200,     # break X axis in time intervals by 500.
        ggtheme = theme_light(), # customize plot and risk table with a theme.
        # risk.table.y.text.col = T,# colour risk table text annotations.
        # risk.table.height = 0.25, # the height of the risk table
        # risk.table.y.text = FALSE,# show bars instead of names in text annotations
        # in legend of risk table.
        # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
        # ncensor.plot.height = 0.25,
        # conf.int.style = "step",  # customize style of confidence intervals
        
        # surv.median.line = "hv",  # add the median survival pointer.
        legend.labs = 
                c(1:4)    # change legend labels.
)
ggsurv <- ggpar(
        ggsurv,
        # font.title    = c(16, "bold", "darkblue"),         
        # font.subtitle = c(15, "bold.italic", "purple"), 
        font.caption  = c(14, "#E41A1C", "#4DAF4A", "#377EB8"),
        # font.x        = c(14, "bold.italic", "red"),          
        # font.y        = c(14, "bold.italic", "darkred"),      
        font.xtickslab = c(20, "plain", "#9C9C9C"),
        font.ytickslab = c(20, "plain", "#9C9C9C"),
        legend = "top"
)
setwd("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\plots_for_paper\\tiff plots from R\\Figure_6_Survival")
# tiff("func_model_survival_curve_10_22_2020.tiff", units="in", width=15, height=8, res=300)
ggsurv
# dev.off()
