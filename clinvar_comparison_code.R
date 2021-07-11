library(easyGgplot2)
library("RColorBrewer")
library(pROC)
library(ggplot2)
library(scales)
###############
## Load Data ##
###############

### load clinvar missense mutations
clinvar <- read.csv("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\clinvar_oncokb_analysis\\clinvar_TP53_Missense_table.csv")
# clinvar_old <- read.csv("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\clinvar_oncokb_analysis\\clinvar_UMD_merge.csv")

### load data with predictions (no protein change duplications)
input.path2<- "C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\GBM- Test final model"
all_data.path = file.path(input.path2, "\\all.mutations.predictions_with_features.test.functional.csv", fsep = .Platform$file.sep)
all_data<-read.table(all_data.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
row.names(all_data) <- all_data[,1]
colnames(all_data)[1] <- "cDNA"



###############
## Functions ##
###############

choose.features<-function(dat, comp, act, vivo, includeLabel=F) {
        
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
        # algorithm's properties (some require labels seperately)
        if (includeLabel==T) {
                cols<-c(cols, "label")
        }
        dat <- dat[,cols]
        
        return(dat)
}

create_mutational_event <- function(clinvar_row) {
        str <- toString(clinvar_row[1])
        ind <- which(strsplit(str, "")[[1]]==">")
        clinvar_row["Mutational_event"] <- substr(str, (ind-1), (ind+1))
}

# remove part of string from "(" to end of string.
# e.g. for 'Pathogenic(Last reviewed: Apr 2, 2018)', keep 'pathogenic' 
remove_part_of_string <- function(row) {
        # row <- orig.dat_clinvar[1,]
        str <- toString(row["clinvar"])
        # split str every character (strsplit), find the index of "(" (which).
        ind <- which(strsplit(str, "")[[1]]=="(")
        if(sum(strsplit(str, "")[[1]]=="(")==0)
                ind <- nchar(str)+1
        row["clinvar"] <- substr(str, 1, (ind-1))
}

##############
#### MAIN ####
##############


### OLD ###
### clinvar
# make sure all variants are in chromosome 17. (indeed they are)
# sum(clinvar$GRCh38Chromosome == 17)
# create mutational event for clinvar
# clinvar$Mutational_event <- apply(clinvar, 1, create_mutational_event)
# # create a merging column
# clinvar$HG38wvar <- paste0(clinvar$GRCh38Location, clinvar$Mutational_event)
# ### orig.dat
# # modify rownames of orig.dat for after the merge.
# orig.dat$cDNA <- rownames(orig.dat)
# # create a merging column
# orig.dat$HG38wvar <- paste0(orig.dat$HG38_Start, orig.dat$Mutational_event)
# ### orig.dat_clinvar: merge clinvar + orig.dat
# orig.dat_clinvar <- merge(orig.dat, clinvar, "HG38wvar")


### NEW ###
# clinvar: 778
## merge clinvar+orig.dat

# include only unique protein change
orig.dat<-all_data
# clinvar's duplicates: 26*2 = 52
clinvar_dups <- clinvar[clinvar$Protein.change.1 %in% clinvar$Protein.change.1[duplicated(clinvar$Protein.change.1)],]
### Test duplicates for conflict interpretations
clinvar_dups <- clinvar_dups[,c("Protein.change.1", "New.Clinical.significance..Last.reviewed.", "Change", "Data.change")]
# remove 5 variants due to conflict in label
remove_from_clinvar <- c("p.C176R", "p.K132N", "p.M246L", "p.L35F", "p.M237I") # change from old analysis: p.V272L, p.V216L removed. p.C176R, p.K132N added. Gil, 10.18.20
clinvar <- clinvar[!clinvar$Protein.change.1 %in% remove_from_clinvar,]
# clinvar: 768
clinvar <- clinvar[!duplicated(clinvar$Protein.change.1),]
# 21 duplicates removed. clinvar: 747
# 2 conflicting labels: p.R248W defined LP, p.R158L defined P/LP. no adjustments needed.

clinvar$Protein.change.1 <- as.character(clinvar$Protein.change.1)
# change p.R72 to p.P72.
clinvar$Protein.change.1[grepl("p.R72", clinvar$Protein.change.1, fixed = T)] <- c("p.P72C", "p.P72S", "p.P72T", "p.P72A", "p.P72H")

# adjust column name
colnames(clinvar)[which(grepl("Protein.change.1", colnames(clinvar), fixed = T))] <- "Protein_p1_TP53" 

# in clinvar not in UMD:
# 3 variants: p.E271L, p.A159F, p.P72C
clinvar_not_in_UMD <- clinvar[!clinvar$Protein_p1_TP53 %in% orig.dat$Protein_p1_TP53,]


### Load chosen model (change 'functional' to 'all_features' or to 'comp', change directory)
current_model <- "functional"

input.path2<- "C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\GBM- Test final model"
all_data.path = file.path(input.path2, "\\all.mutations.predictions_with_features.test.functional.csv", fsep = .Platform$file.sep) # change this
all_data<-read.table(all_data.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
row.names(all_data) <- all_data[,1]
colnames(all_data)[1] <- "cDNA"


# merge clinvar with orig.dat
orig.dat_clinvar <- merge(all_data, clinvar, by = "Protein_p1_TP53")

# take only relevant parts: cDNA and clinical significance
orig.dat_clinvar <- data.frame(orig.dat_clinvar$cDNA, orig.dat_clinvar$Protein_p1_TP53, orig.dat_clinvar$New.Clinical.significance..Last.reviewed.)
colnames(orig.dat_clinvar) <- c("cDNA", "Protein change", "clinvar")


# adjust clinvar.
orig.dat_clinvar$`clinvar` <- apply(orig.dat_clinvar, 1, remove_part_of_string)
orig.dat_clinvar$`clinvar`[orig.dat_clinvar$`clinvar` == "Uncertain significance" | orig.dat_clinvar$`clinvar` == "Conflicting interpretations of pathogenicity"] <- "VUS"
orig.dat_clinvar$`clinvar`[orig.dat_clinvar$`clinvar` == "not provided"] <- "No data"
orig.dat_clinvar$`clinvar`[orig.dat_clinvar$`clinvar` == "Pathogenic/Likely pathogenic"] <- "Pathogenic"
orig.dat_clinvar$`clinvar`[orig.dat_clinvar$`clinvar` == "drug response"] <- "No data"
orig.dat_clinvar$`clinvar`[orig.dat_clinvar$`clinvar` == "Likely pathogenic, drug response"] <- "Likely pathogenic"
orig.dat_clinvar$`clinvar`[orig.dat_clinvar$`clinvar` == "Benign/Likely benign"] <- "Likely benign"
# create factor with desired order of levels.
# levels(as.factor(orig.dat_clinvar$`clinvar`))
orig.dat_clinvar$`clinvar` <- factor(orig.dat_clinvar$`clinvar`, levels = c("Pathogenic", "Likely pathogenic", "Benign", "Likely benign", "VUS", "No data"))


### all data: model's predictions
# remove unnecessary columns from all_data
all_data <- all_data[,c("Protein_p1_TP53", "Set", "label", "predictions", "score")]
colnames(all_data)[1] <- "Protein change"
# Validation is considered training set.
all_data$Set[all_data$Set == "Validation"] <- "train"

### all_data_clinvar: merge orig.dat_clinvar with all_data, for labels and predictions
all_data_clinvar <- merge(all_data, orig.dat_clinvar, "Protein change")
row.names(all_data_clinvar) <- all_data_clinvar$cDNA
# seperate all_data_clinvar to train (D/ND) and predictions
clinvar_train_D <- all_data_clinvar[all_data_clinvar$label == "D",]
clinvar_train_ND <-all_data_clinvar[all_data_clinvar$label == "ND",]
clinvar_pred <- all_data_clinvar[all_data_clinvar$Set == "prediction",]

#take only mutations that DON'T appear in clinvar
all_data_not_clinvar <- all_data[!(all_data$`Protein change` %in% orig.dat_clinvar$`Protein change`),]
# add columns to all_data so that rows can be attached to all_data_clinvar
all_data_not_clinvar$cDNA <- row.names(all_data_not_clinvar)
all_data_not_clinvar$'clinvar' <- "No data"
# all_data_not_clinvar$HG38wvar <- NA
# seperate all_data_not_clinvar to D/ND labels. 
# this is only Train (+ Val) Set! (since predictions don't have label)
D_all_data <- all_data_not_clinvar[all_data_not_clinvar$label == "D",] # 31 mutations.
ND_all_data <- all_data_not_clinvar[all_data_not_clinvar$label == "ND",] # 882 mutations.
pred_all_data <- all_data_not_clinvar[all_data_not_clinvar$Set == "prediction",] # 657 mutations.

###################################
####### 2020_12_08_UPDATE: ########
###################################
# Thierry's modifications to figure 7
# apply next rows (203-219) to create supp.figure 4 for final paper.
# clinvar_train_D <- clinvar_train_D[clinvar_train_D$clinvar != "No data",]
# clinvar_train_ND <- clinvar_train_ND[clinvar_train_ND$clinvar != "No data",]
# clinvar_train_D$`clinvar` <- factor(clinvar_train_D$`clinvar`, levels = c("Benign", "Likely benign", "VUS", "Pathogenic", "Likely pathogenic"))
# clinvar_train_ND$clinvar <- factor(clinvar_train_ND$`clinvar`, levels = c("Benign", "Likely benign", "VUS", "Pathogenic", "Likely pathogenic"))
# # plot D:
# setwd("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\plots_for_paper\\tiff plots from R\\Figure_7_Clinvar")
# temp1 <- data.frame(table(clinvar_train_D$clinvar))
# # save plot
# tiff("Clinvar_D_barplot_all_features_2020_12_08.tiff", units="in", width=10, height=5, res=300)
# ggplot(temp1, aes(x=Var1, y=Freq, fill=Var1)) + scale_fill_manual(values = c("#FFED6F", "#FFED6F", "#BC80BD", "#80B1D3", "#BEBADA")) + geom_bar(position="dodge", stat="identity")
# dev.off()
# # plot ND:
# temp2 <- data.frame(table(clinvar_train_ND$clinvar))
# # Save plot
# tiff("Clinvar_ND_barplot_2020_12_08.tiff", units="in", width=10, height=5, res=300)
# ggplot(temp2, aes(x=Var1, y=Freq, fill=Var1)) + scale_fill_manual(values = c("#FFED6F", "#FFED6F", "#BC80BD", "#80B1D3", "#BEBADA")) + geom_bar(position="dodge", stat="identity")
# dev.off()



#attach them to all_data_clinvar
clinvar_train_D <- rbind(clinvar_train_D, D_all_data)
# unique(clinvar_train_D$clinvar)
clinvar_train_D$`clinvar` <- factor(clinvar_train_D$`clinvar`, levels = c("Pathogenic", "Likely pathogenic", "VUS"))

clinvar_train_ND <- rbind(clinvar_train_ND, ND_all_data)
# unique(clinvar_train_ND$clinvar)
clinvar_train_ND$clinvar <- factor(clinvar_train_ND$`clinvar`, levels = c("Likely benign", "VUS", "No data"))

clinvar_pred <- rbind(clinvar_pred, pred_all_data)
# unique(clinvar_pred$clinvar)
clinvar_pred$clinvar <- factor(clinvar_pred$clinvar, levels = c("Pathogenic", "Likely pathogenic", "Benign", "Likely benign", "VUS", "No data"))


### Plot training data (D/ND)### 

# see colors. change name to: Set1/Set2/Set3
# look here for more: http://www.sthda.com/english/wiki/colors-in-r#using-rcolorbrewer-palettes
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n = 12, name = "Set3")

# plot D:
setwd("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\plots_for_paper\\tiff plots from R\\Figure_7_Clinvar")
clinvar_D_barplot <- ggplot2.barplot(data=clinvar_train_D, xName='clinvar',groupName="clinvar", groupColors =c("#80B1D3", "#BEBADA", "#FFFFB3","#D9D9D9")) + ggtitle("Training set deleterious mutations compared to ClinVar")
tiff("Clinvar_D_barplot_all_features_10_23_2020.tiff", units="in", width=10, height=5, res=300)
ggplot2.setAxis(clinvar_D_barplot, xtitle = "Clinical Significance (ClinVar)", ytitle = "Count")
dev.off()

# plot ND:
clinvar_ND_barplot <- ggplot2.barplot(data=clinvar_train_ND, xName='clinvar', groupName='clinvar', groupColors =c("#4DAF4A",  "#FFFFB3","#D9D9D9")) + ggtitle("Training set non-deleterious mutations compared to ClinVar")

tiff("Clinvar_ND_barplot_10_18_2020.tiff", units="in", width=10, height=5, res=300)
ggplot2.setAxis(clinvar_ND_barplot, xtitle = "Clinical Significance (ClinVar)", ytitle = "Counts")
dev.off()

# plots in numbers:
sum(clinvar_train_D$clinvar == "Pathogenic")
sum(clinvar_train_D$clinvar == "Likely pathogenic")
sum(clinvar_train_D$clinvar == "VUS")
sum(clinvar_train_D$clinvar == "No data")

sum(clinvar_train_ND$clinvar == "Benign")
sum(clinvar_train_ND$clinvar == "Likely benign")
sum(clinvar_train_ND$clinvar == "VUS")
sum(clinvar_train_ND$clinvar == "No data")

sum(clinvar_pred$clinvar == "Pathogenic")
sum(clinvar_pred$clinvar == "Likely pathogenic")
sum(clinvar_pred$clinvar == "Likely benign")
sum(clinvar_pred$clinvar == "Benign")
sum(clinvar_pred$clinvar == "VUS")
sum(clinvar_pred$clinvar == "No data")

### predictions set
# make adjustments in clinvar_pred
# colnames(clinvar_pred)[7] <- "clinvar"
clinvar_pred$clinvar <- as.character(clinvar_pred$clinvar)
clinvar_pred$clinvar[clinvar_pred$clinvar == "Pathogenic" | clinvar_pred$clinvar=="Likely pathogenic"] <- "D"
clinvar_pred$clinvar[clinvar_pred$clinvar == "Likely benign" | clinvar_pred$clinvar == "Benign"] <- "ND"

# calculate confusion matrix for predictions:
sum(clinvar_pred$predictions == "D" & clinvar_pred$clinvar == "D")
sum(clinvar_pred$predictions == "D" & clinvar_pred$clinvar == "ND")
sum(clinvar_pred$predictions == "D" & clinvar_pred$clinvar == "VUS")
sum(clinvar_pred$predictions == "D" & clinvar_pred$clinvar == "No data")

sum(clinvar_pred$predictions == "ND" & clinvar_pred$clinvar == "D")
sum(clinvar_pred$predictions == "ND" & clinvar_pred$clinvar == "ND")
sum(clinvar_pred$predictions == "ND" & clinvar_pred$clinvar == "VUS")
sum(clinvar_pred$predictions == "ND" & clinvar_pred$clinvar == "No data")

# create table of false negatives
clinvar_model_false_negatives <- clinvar_pred[clinvar_pred$predictions == "ND" & clinvar_pred$clinvar == "D",]
# write.csv(clinvar_model_false_negatives_comp, "C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\clinvar_oncokb_analysis\\clinvar_model_comp_false_negatives_10_18_2020.csv")

# calculate ROC for clinvar mutations with defined clinical effect (D/ND).
predict_clinvar_global <- rbind(clinvar_pred, clinvar_train_D, clinvar_train_ND)
predict_clinvar <- predict_clinvar_global[predict_clinvar_global$clinvar == "D" | predict_clinvar_global$clinvar == "ND",]
predict_clinvar_global$clinvar <- as.character(predict_clinvar_global$clinvar)
predict_clinvar_global$clinvar[predict_clinvar_global$clinvar == "Pathogenic" | predict_clinvar_global$clinvar=="Likely pathogenic"] <- "D"
predict_clinvar_global$clinvar[predict_clinvar_global$clinvar == "Benign" | predict_clinvar_global$clinvar == "Likely benign"] <- "ND"
predict_clinvar_global$predictions[predict_clinvar_global$predictions == "null"] <- predict_clinvar_global$label[predict_clinvar_global$predictions == "null"]
predict_clinvar_global_all <- predict_clinvar_global 

# confusion matrix for clinvar pred+train (clinvar_global)
sum(predict_clinvar_global_all$predictions == "D" & predict_clinvar_global_all$clinvar == "D")
sum(predict_clinvar_global_all$predictions == "D" & predict_clinvar_global_all$clinvar == "ND")
sum(predict_clinvar_global_all$predictions == "D" & predict_clinvar_global_all$clinvar == "VUS")
sum(predict_clinvar_global_all$predictions == "D" & predict_clinvar_global_all$clinvar == "No data")

sum(predict_clinvar_global_all$predictions == "ND" & predict_clinvar_global_all$clinvar == "D")
sum(predict_clinvar_global_all$predictions == "ND" & predict_clinvar_global_all$clinvar == "ND")
sum(predict_clinvar_global_all$predictions == "ND" & predict_clinvar_global_all$clinvar == "VUS")
sum(predict_clinvar_global_all$predictions == "ND" & predict_clinvar_global_all$clinvar == "No data")

# take only what's relevant for AUROC.
predict_clinvar_global <- predict_clinvar_global[predict_clinvar_global$clinvar == "D" | predict_clinvar_global$clinvar == "ND",c("Protein change", "cDNA", "clinvar", "score", "predictions")]


# compute ROC for clinvar-model 
# include training
ROC_global <- roc(predictor=predict_clinvar_global$score,
           response=predict_clinvar_global$clinvar,
           levels=rev(levels(as.factor(as.character(predict_clinvar_global$clinvar)))))
print(paste0("AUC ", ROC_global$auc))
# predictions only.

ROC_pred <- roc(predictor=predict_clinvar$score,
           response=predict_clinvar$clinvar,
           levels=rev(levels(as.factor(as.character(predict_clinvar$clinvar)))))
print(paste0("AUC ", ROC_pred$auc))


# occurence <- orig.dat[,c("cDNA", "TCGA_Freq", "ICGC_Freq", "Sanger_freq", "MSKCC_Freq", "GNOMAD_Freq")]
# occurence <- merge(clinvar_model_false_negatives, splice, "cDNA")

# false negatives analysis
input.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53"
FEATURES_DF_PATH = file.path(input.path, "data_list.csv", fsep = .Platform$file.sep)
df_features = read.table(FEATURES_DF_PATH, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
features = df_features[grepl("Yes",df_features$`as a feature`),]$'Field Name'
features.act<-df_features["Yes (activity)" == df_features$`as a feature`,]$'Field Name'
features.vivo <- df_features["Yes (vivo)" == df_features$'as a feature',]$'Field Name'
features.no.act <-setdiff(features, features.act)
features.no.vivo <- setdiff(features, features.vivo)
features.comp <- setdiff(features.no.act, features.vivo)
features.no.comp <- c(features.act, features.vivo)

row.names(orig.dat) <- orig.dat$Protein_p1_TP53
orig.dat_features <- choose.features(orig.dat, T, T, T, F)
compare_AUROC <- orig.dat_features[predict_clinvar$`Protein change`, c("Envision_predictions", "CADD_raw_Rankscore" ,"Polyphen2_HDIV_score_cor", "Polyphen2_HVAR_score_cor", "SIFT_score_cor")]
orig.dat_features$"Protein change" <- orig.dat$Protein_p1_TP53
compare_AUROC <- merge(orig.dat_features, predict_clinvar, by = "Protein change")

# PolyPhen2 HDIV
ROC_polyHDIV <- roc(predictor=compare_AUROC$Polyphen2_HDIV_score_cor,
                response=compare_AUROC$clinvar,
                levels=rev(levels(as.factor(as.character(compare_AUROC$clinvar)))))
print(paste0("AUC ", ROC_polyHDIV$auc))
# PolyPhen2 HVAR
ROC_polyHVAR <- roc(predictor=compare_AUROC$Polyphen2_HVAR_score_cor,
                    response=compare_AUROC$clinvar,
                    levels=rev(levels(as.factor(as.character(compare_AUROC$clinvar)))))
print(paste0("AUC ", ROC_polyHVAR$auc))
# CADD
ROC_CADD <- roc(predictor=compare_AUROC$CADD_raw_Rankscore,
                    response=compare_AUROC$clinvar,
                    levels=rev(levels(as.factor(as.character(compare_AUROC$clinvar)))))
print(paste0("AUC ", ROC_CADD$auc))
# SIFT
ROC_SIFT <- roc(predictor=compare_AUROC$SIFT_score_cor,
                response=compare_AUROC$clinvar,
                levels=rev(levels(as.factor(as.character(compare_AUROC$clinvar)))))
print(paste0("AUC ", ROC_SIFT$auc))
# Envision
ROC_Env <- roc(predictor=as.numeric(compare_AUROC$Envision_predictions),
                response=compare_AUROC$clinvar,
                levels=rev(levels(as.factor(as.character(compare_AUROC$clinvar)))))
print(paste0("AUC ", ROC_Env$auc))

# compare scores for global
compare_AUROC_global <- orig.dat_features[predict_clinvar_global$`Protein change`, c("Envision_predictions", "CADD_raw_Rankscore" ,"Polyphen2_HDIV_score_cor", "Polyphen2_HVAR_score_cor", "SIFT_score_cor")]
compare_AUROC_global <- merge(orig.dat_features, predict_clinvar_global, "Protein change")
compare_AUROC_global <- cbind(predict_clinvar_global, compare_AUROC_global)


# PolyPhen2 HDIV
ROC_polyHDIV_global <- roc(predictor=compare_AUROC_global$Polyphen2_HDIV_score_cor,
                    response=compare_AUROC_global$clinvar,
                    levels=levels(as.factor(as.character(compare_AUROC_global$clinvar))))
print(paste0("AUC ", ROC_polyHDIV_global$auc))
# PolyPhen2 HVAR
ROC_polyHVAR_global <- roc(predictor=compare_AUROC_global$Polyphen2_HVAR_score_cor,
                    response=compare_AUROC_global$clinvar,
                    levels=levels(as.factor(as.character(compare_AUROC_global$clinvar))))
print(paste0("AUC ", ROC_polyHVAR_global$auc))
# CADD
ROC_CADD_global <- roc(predictor=compare_AUROC_global$CADD_raw_Rankscore,
                response=compare_AUROC_global$clinvar,
                levels=rev(levels(as.factor(as.character(compare_AUROC_global$clinvar)))))
print(paste0("AUC ", ROC_CADD_global$auc))
# SIFT
ROC_SIFT_global <- roc(predictor=compare_AUROC_global$SIFT_score_cor,
                response=compare_AUROC_global$clinvar,
                levels=rev(levels(as.factor(as.character(compare_AUROC_global$clinvar)))))
print(paste0("AUC ", ROC_SIFT_global$auc))
# Envision
ROC_Env_global <- roc(predictor=as.numeric(compare_AUROC_global$Envision_predictions),
               response=compare_AUROC_global$clinvar,
               levels=rev(levels(as.factor(as.character(compare_AUROC_global$clinvar)))))
print(paste0("AUC ", ROC_Env_global$auc))

### SIFT compared to Model ###
sum(orig.dat$Sift_Prediction == "Damaging" & orig.dat$predictions == "D")
sum(orig.dat$Sift_Prediction == "Damaging" & orig.dat$predictions == "ND")
sum(orig.dat$Sift_Prediction == "Tolerated" & orig.dat$predictions == "D")
sum(orig.dat$Sift_Prediction == "Tolerated" & orig.dat$predictions == "ND")

### SIFT compared to ClinVar ###
orig.dat$"Protein change" <- orig.dat$Protein_p1_TP53
SIFT_clinvar_comparison <- merge(predict_clinvar, orig.dat, "Protein change")
sum(SIFT_clinvar_comparison$clinvar == "D" & SIFT_clinvar_comparison$Sift_Prediction == "Damaging")
sum(SIFT_clinvar_comparison$clinvar == "ND" & SIFT_clinvar_comparison$Sift_Prediction == "Damaging")
sum(SIFT_clinvar_comparison$clinvar == "D" & SIFT_clinvar_comparison$Sift_Prediction == "Tolerated")
sum(SIFT_clinvar_comparison$clinvar == "ND" & SIFT_clinvar_comparison$Sift_Prediction == "Tolerated")

SIFT_clinvar_comparison_global <- merge(predict_clinvar_global, orig.dat, "Protein change")
sum(SIFT_clinvar_comparison_global$clinvar == "D" & SIFT_clinvar_comparison_global$Sift_Prediction == "Damaging")
sum(SIFT_clinvar_comparison_global$clinvar == "D" & SIFT_clinvar_comparison_global$Sift_Prediction == "Tolerated")
sum(SIFT_clinvar_comparison_global$clinvar == "ND" & SIFT_clinvar_comparison_global$Sift_Prediction == "Damaging")
sum(SIFT_clinvar_comparison_global$clinvar == "ND" & SIFT_clinvar_comparison_global$Sift_Prediction == "Tolerated")

### create AUROC values for models 5.2-5.5
# load models and environment.
load("C:/Users/USER/Desktop/gil/lab/TP53/workspaces/Model_Results_Neg_Set_corrected_30_03_2020/model_results_neg_set_corrected_validation_models_tuned_02_04_2020.RData")

get_scores_from_model <- function(model, to_pred, comp, act, vivo) {
        predsSet <- orig.dat
        predsSet <- predsSet[!duplicated(predsSet$Protein_p1_TP53),]
        predsSet <- choose.features(dat = predsSet, comp = comp, act = act, vivo = vivo, includeLabel = F)
        predsSet[1:(ncol(predsSet)-1)]=data.frame(apply(predsSet[1:(ncol(predsSet)-1)],2,orderColumns))
        scores <-  predict(model,newdata = predsSet, na.action = na.pass, type="prob")
        scores <- scores$D
        cDNA = row.names(predsSet)
        protein_change <- orig.dat[cDNA,]$Protein_p1_TP53
        data <- data.frame(cDNA, protein_change, scores)
        row.names(data) <- protein_change
        return(data[to_pred,]$scores)
}
predict_clinvar$gbm_only_score <- get_scores_from_model(gbm_test_model5.5, predict_clinvar$`Protein change`, T, F, F)
predict_clinvar$score5.1 <- get_scores_from_model(gbm_test_model5.1, predict_clinvar$`Protein change`, T, T, T)
predict_clinvar$score5.2 <- get_scores_from_model(gbm_model5.2, predict_clinvar$`Protein change`, T, T, F)
predict_clinvar$score5.3 <- get_scores_from_model(gbm_model5.3, predict_clinvar$`Protein change`, T, F, T)
predict_clinvar$score5.4 <- get_scores_from_model(gbm_model5.4, predict_clinvar$`Protein change`, F, T, T)

# comp
ROC_comp_only <- roc(predictor=as.numeric(predict_clinvar$gbm_only_score),
               response=predict_clinvar$clinvar,
               levels=rev(levels(as.factor(as.character(predict_clinvar$clinvar)))))
print(paste0("AUC ", ROC_comp_only$auc))
# comp act
ROC_comp_act <- roc(predictor=as.numeric(predict_clinvar$score5.2),
                     response=predict_clinvar$clinvar,
                     levels=rev(levels(as.factor(as.character(predict_clinvar$clinvar)))))
print(paste0("AUC ", ROC_comp_act$auc))
# comp vivo
ROC_comp_vivo <- roc(predictor=as.numeric(predict_clinvar$score5.3),
                    response=predict_clinvar$clinvar,
                    levels=rev(levels(as.factor(as.character(predict_clinvar$clinvar)))))
print(paste0("AUC ", ROC_comp_vivo$auc))
ROC_act_vivo <- roc(predictor=as.numeric(predict_clinvar$score5.4),
                     response=predict_clinvar$clinvar,
                     levels=rev(levels(as.factor(as.character(predict_clinvar$clinvar)))))
print(paste0("AUC ", ROC_act_vivo$auc))

# predict clinvar_global
predict_clinvar_global$gbm_only_score <- get_scores_from_model(gbm_test_model5.5, predict_clinvar_global$`Protein change`, T, F, F)
predict_clinvar_global$score5.1 <- get_scores_from_model(gbm_test_model5.1, predict_clinvar_global$`Protein change`, T, T, T)
predict_clinvar_global$score5.2 <- get_scores_from_model(gbm_model5.2, predict_clinvar_global$`Protein change`, T, T, F)
predict_clinvar_global$score5.3 <- get_scores_from_model(gbm_model5.3, predict_clinvar_global$`Protein change`, T, F, T)
predict_clinvar_global$score5.4 <- get_scores_from_model(gbm_model5.4, predict_clinvar_global$`Protein change`, F, T, T)

# comp
ROC_comp_only <- roc(predictor=as.numeric(predict_clinvar_global$gbm_only_score),
                     response=predict_clinvar_global$clinvar,
                     levels=rev(levels(as.factor(as.character(predict_clinvar_global$clinvar)))))
print(paste0("AUC ", ROC_comp_only$auc))
# comp act
ROC_comp_act <- roc(predictor=as.numeric(predict_clinvar_global$score5.2),
                    response=predict_clinvar_global$clinvar,
                    levels=rev(levels(as.factor(as.character(predict_clinvar_global$clinvar)))))
print(paste0("AUC ", ROC_comp_act$auc))
# comp vivo
ROC_comp_vivo <- roc(predictor=as.numeric(predict_clinvar_global$score5.3),
                     response=predict_clinvar_global$clinvar,
                     levels=rev(levels(as.factor(as.character(predict_clinvar_global$clinvar)))))
print(paste0("AUC ", ROC_comp_vivo$auc))
ROC_act_vivo <- roc(predictor=as.numeric(predict_clinvar_global$score5.4),
                    response=predict_clinvar_global$clinvar,
                    levels=rev(levels(as.factor(as.character(predict_clinvar_global$clinvar)))))
print(paste0("AUC ", ROC_act_vivo$auc))


# PIE Charts for model comparison with clinvar
# functional/all_features model (they are the same, gil 18_10_2020)
df <- data.frame(
        group = c("True Positive", "True Negative", "False Negative"),
        Variants = c(157, 26, 7)
)

pie<- ggplot(df, aes(x="", y=Variants, fill=group))+
        geom_bar(width = 1, stat = "identity")  + 
        ggtitle("All features model predictions compared with ClinVar annotations") +
        coord_polar("y", start=0) +
        scale_fill_manual(values=c("#E41A1C", "#B3DE69",  "#4DAF4A")) +
        theme(axis.text.x=element_blank(), plot.title = element_text(size=15))
pie
setwd("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\plots_for_paper\\tiff plots from R")
# tiff("Pie_Chart_compare_functional_model_to_clinvar_10_18_2020.tiff", units="in", width=10, height=5, res=300)
pie
# dev.off()


# computational model
df <- data.frame(
        group = c("True Positive", "True Negative", "False Negative"),
        Variants = c(183, 20, 4)
)

pie <- ggplot(df, aes(x="", y=Variants, fill=group))+
        geom_bar(width = 1, stat = "identity")  + 
        ggtitle("Model predictions compared with ClinVar annotations") +
        coord_polar("y", start=0) +
        scale_fill_manual(values=c("#E41A1C", "#B3DE69",  "#4DAF4A")) +
        theme(axis.text.x=element_blank(), plot.title = element_text(size=15))
pie
setwd("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\plots_for_paper\\tiff plots from R")
# tiff("Pie_Chart_compare_comp_model_to_clinvar_19_05_2020.tiff", units="in", width=10, height=5, res=300)
pie
# dev.off()





### heatmaps ###
orig.dat_features <- apply(orig.dat_features, 2, as.numeric)
orig.dat_features <- as.matrix(orig.dat_features)
range01 <- function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}
orig.dat_features<-apply(orig.dat_features,2, range01 )
row.names(orig.dat_features) <- row.names(orig.dat)

false_neg_features <- orig.dat_features[orig.dat$HG38wvar %in% clinvar_model_false_negatives$HG38wvar,]
false_neg_features <- as.matrix(apply(false_neg_features, 2, as.numeric))
for_heatmap <- rbind(orig.dat_features, false_neg_features)
# tiff("heatmap.tiff", units="in", width=15, height=8, res=600)
heatmap(false_neg_features, cexCol=0.3, cexRow = 1)
# dev.off()

median_orig.dat <- apply(orig.dat_features, 2, function(x){median(x, na.rm=T)})
train_D_features <- orig.dat_features[row.names(orig.dat) %in% clinvar_train_D$cDNA,]
train_ND_features <- orig.dat_features[row.names(orig.dat) %in% clinvar_train_ND$cDNA,]

median_D <- apply(train_D_features, 2, function(x){median(x, na.rm=T)})
median_ND <- apply(train_ND_features, 2, function(x){median(x, na.rm=T)})
for_heatmap <- rbind(false_neg_features, median_orig.dat, median_D, median_ND)



orig.dat_features <- choose.features(orig.dat, T, T, T, F)
predict_clinvar_features <- orig.dat_features[row.names(orig.dat_features) %in% predict_clinvar$cDNA,]
predict_clinvar_features[1:(ncol(predict_clinvar_features)-1)]=data.frame(apply(predict_clinvar_features[1:(ncol(predict_clinvar_features)-1)],2,orderColumns))
clinvar_preds <- predict(gbm_test_model5.1, newdata = predict_clinvar_features, na.action = na.pass, type="prob")
clinvar_preds$cDNA <- row.names(predict_clinvar_features)
predict_clinvar <- merge(clinvar_preds, predict_clinvar, "cDNA")
predict_clinvar$clinvar <- as.factor(predict_clinvar$clinvar)
temp_model <- gbm_test_model5.1
temp_model$cutoff<-calc.cutoff(predict_clinvar$D.x, labels=predict_clinvar$clinvar, sens.levepredict_clinvar_features)
                               
ROC <- roc(predictor=predict_clinvar$score,
           response=predict_clinvar$clinvar,
           levels=rev(levels(as.factor(as.character(predict_clinvar$clinvar)))))
print(paste0("AUC ", ROC$auc))

predicted <- predict.cutoff(gbm_test_model5.1, newdata = predict_clinvar_features, gbm_test_model5.1$cutoff)
predicted <- factor(predicted, c("D", "ND"))
predicted <- data.frame(row.names(predict_clinvar_features), predicted)
colnames(predicted) <- c("cDNA", "predict")
predict_clinvar <- merge(predict_clinvar, predicted, "cDNA")
print(caret::confusionMatrix(predict_clinvar$predict, predict_clinvar$clinvar))

### Test 41 variants comparing functional model to SIFT.
# load all_data, fix row names.
input.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\GBM- Test final model"
all_data.path = file.path(input.path, "all.mutations.predictions_with_features.test.functional.csv", fsep = .Platform$file.sep)
all_data<-read.table(all_data.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
row.names(all_data) <- all_data[,1]
all_data <- all_data[,-1]

scores_for_comparison <- all_data[row.names(all_data) %in% variants_41$mutation,]
scores_for_comparison <- scores_for_comparison[,c("Sift_Prediction", "SIFT_score_cor", "Polyphen_2_HumVar", "Polyphen2_HVAR_score_cor","Polyphen_2_HumDiv", "Polyphen2_HDIV_score_cor", "Provean_Prediction", "PROVEAN_score_cor", "Mutassessor_Prediction", "Mutassessor_Score","Condel", "Condel_Score")]
scores_for_comparison$mutation <- row.names(scores_for_comparison)

# load 41 variants with model scores and predictions.
variants_41 <- read.csv("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\41 set\\41_variants_3_models_predictions.csv")
### merge with scores_for_comparison ###
variants_41 <- merge(variants_41, scores_for_comparison, by = "mutation")
# create label for the true variants position.
variants_41$true_label <- variants_41$predictions_func
variants_41$true_label[variants_41$protein_change == "p.G154S" | variants_41$protein_change == "p.R156H"] <- "ND"

apply(variants_41, 2, unique)

# calculate ACC, AUROC
# Sift
variants_41$Sift_Prediction[variants_41$Sift_Prediction == "Damaging"] <- "D"
variants_41$Sift_Prediction[variants_41$Sift_Prediction == "Tolerated"] <- "ND"
# Mutassessor
variants_41$Mutassessor_Prediction[variants_41$Mutassessor_Prediction == "Neutral"] <- "ND"
variants_41$Mutassessor_Prediction[variants_41$Mutassessor_Prediction == "neutral"] <- "ND"
variants_41$Mutassessor_Prediction[variants_41$Mutassessor_Prediction == "Medium"] <- "D"
variants_41$Mutassessor_Prediction[variants_41$Mutassessor_Prediction == "Low"] <- "ND"
# Polyphen
variants_41$Polyphen_2_HumVar[variants_41$Polyphen_2_HumVar == "Benign"] <- "ND"
variants_41$Polyphen_2_HumVar[variants_41$Polyphen_2_HumVar == "Probably damaging"] <- "D"
variants_41$Polyphen_2_HumVar[variants_41$Polyphen_2_HumVar == "Possibly damaging"] <- "D"

variants_41$Polyphen_2_HumDiv[variants_41$Polyphen_2_HumDiv == "Benign"] <- "ND"
variants_41$Polyphen_2_HumDiv[variants_41$Polyphen_2_HumDiv == "Probably damaging"] <- "D"
variants_41$Polyphen_2_HumDiv[variants_41$Polyphen_2_HumDiv == "Possibly damaging"] <- "D"
# Provean
variants_41$Provean_Prediction[variants_41$Provean_Prediction == "Deleterious"] <- "D"
variants_41$Provean_Prediction[variants_41$Provean_Prediction == "Neutral"] <- "ND"
# Condel
variants_41$Condel[variants_41$Condel == "Deleterious"] <- "D"
variants_41$Condel[variants_41$Condel == "Neutral"] <- "ND"

# Acc:
sum(variants_41$Sift_Prediction == variants_41$true_label)/41 
sum(variants_41$Mutassessor_Prediction == variants_41$true_label)/41 
sum(variants_41$Provean_Prediction == variants_41$true_label)/41 
# remove NAs from polyphen, condel (3 NAs)
sum(variants_41$Polyphen_2_HumVar == variants_41$true_label, na.rm = T)/38 
sum(variants_41$Polyphen_2_HumDiv == variants_41$true_label, na.rm = T)/38 
sum(variants_41$Condel == variants_41$true_label, na.rm = T)/38 
# models
sum(variants_41$predictions_func == variants_41$true_label)/41 
sum(variants_41$predictions_all_features == variants_41$true_label)/41 
sum(variants_41$predictions_comp == variants_41$true_label)/41 

#AUC
# SIFT
ROC <- roc(predictor=variants_41$SIFT_score_cor,
           response=variants_41$true_label,
           levels=rev(levels(as.factor(as.character(variants_41$true_label)))))
print(paste0("AUC ", ROC$auc))
# Mutassesor
ROC <- roc(predictor=variants_41$Mutassessor_Score,
           response=variants_41$true_label,
           levels=rev(levels(as.factor(as.character(variants_41$true_label)))))
print(paste0("AUC ", ROC$auc))
# polyphen HDIV
ROC <- roc(predictor=variants_41$Polyphen2_HDIV_score_cor,
           response=variants_41$true_label, na.rm = T,
           levels=rev(levels(as.factor(as.character(variants_41$true_label)))))
print(paste0("AUC ", ROC$auc))
# polyphen HVAR
ROC <- roc(predictor=variants_41$Polyphen2_HVAR_score_cor,
           response=variants_41$true_label, na.rm = T,
           levels=rev(levels(as.factor(as.character(variants_41$true_label)))))
print(paste0("AUC ", ROC$auc))
# Provean
ROC <- roc(predictor=variants_41$PROVEAN_score_cor,
           response=variants_41$true_label,
           levels=rev(levels(as.factor(as.character(variants_41$true_label)))))
print(paste0("AUC ", ROC$auc))
# Condel
ROC <- roc(predictor=variants_41$Condel_Score,
           response=variants_41$true_label, na.rm = T,
           levels=rev(levels(as.factor(as.character(variants_41$true_label)))))
print(paste0("AUC ", ROC$auc))
# Functional
ROC <- roc(predictor=variants_41$score_func,
           response=variants_41$true_label,
           levels=rev(levels(as.factor(as.character(variants_41$true_label)))))
print(paste0("AUC ", ROC$auc))
# All features
ROC <- roc(predictor=variants_41$score_all_features,
           response=variants_41$true_label,
           levels=rev(levels(as.factor(as.character(variants_41$true_label)))))
print(paste0("AUC ", ROC$auc))
# comp
ROC <- roc(predictor=variants_41$score_comp,
           response=variants_41$true_label,
           levels=rev(levels(as.factor(as.character(variants_41$true_label)))))
print(paste0("AUC ", ROC$auc))


