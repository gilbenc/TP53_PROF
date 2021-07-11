### LFS analysis ###

# 80 variants
LFS_denovo <- read.csv("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\LFS analysis\\LFS_denovo.csv")
# 84 variants
LFS_familial <- read.csv("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\LFS analysis\\LFS_familial.csv")

LFS_familial$status <- "Familial"
LFS_denovo$status <- "Denovo"
colnames(LFS_denovo)[1] <- "Protein_p1_TP53"
colnames(LFS_familial)[1] <- "Protein_p1_TP53"
# 164 variants
LFS <- rbind(LFS_familial, LFS_denovo)


# check for duplications: 22 duplicates appear 3 times each, twice as FM, once denovo.
dups <- LFS[LFS$Protein_p1_TP53 %in% LFS$Protein_p1_TP53[duplicated(LFS$Protein_p1_TP53)],]
# change their status to include both, and remove duplicates. 
LFS$status[LFS$Protein_p1_TP53 %in% dups$Protein_p1_TP53] <- "Denovo/familial"
# dups analysis:
# how many times each variant in dups appear?
sapply(unique(dups$Protein_p1_TP53), function(x){sum(dups$Protein_p1_TP53 == x)})
# how many duplicates of 2 are there? (use same command with ==3).
sum(sapply(unique(dups$Protein_p1_TP53), function(x){sum(dups$Protein_p1_TP53 == x)})==2)
# remove duplicates. 22 variants appear 3 times, 40 appear twice. overall: 84 duplicates removed.
LFS <- LFS[!duplicated(LFS$Protein_p1_TP53),]
# after this change, should have 80 variants.

### load data with predictions (no protein change duplications)
input.path2<- "C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\GBM- Test final model"
all_data.path = file.path(input.path2, "\\all.mutations.predictions_with_features.test.functional.csv", fsep = .Platform$file.sep)
all_data<-read.table(all_data.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
row.names(all_data) <- all_data[,1]
colnames(all_data)[1] <- "cDNA"
all_data_for_merge <- all_data[,c("cDNA", "Protein_p1_TP53","GNOMAD_Rec","Set","label","predictions","score")]

# put aside variants that are not in UMD (2)
LFS_not_in_UMD <- LFS[!LFS$Protein_p1_TP53 %in% all_data_for_merge$Protein_p1_TP53,]
# after merge: 78 variants.
LFS_with_model_predictions <- merge(LFS, all_data_for_merge, "Protein_p1_TP53")

# analysis
# 48 variants in train set
sum(LFS_with_model_predictions$Set == "train")
# 16 in validation set
sum(LFS_with_model_predictions$Set == "Validation")
# 14 in prediction set
sum(LFS_with_model_predictions$Set == "prediction")
# 8 variants predicted ND from prediction set.
sum(LFS_with_model_predictions$Set == "prediction" & LFS_with_model_predictions$predictions == "D")
# 6 variants predicted ND from prediction set.
sum(LFS_with_model_predictions$Set == "prediction" & LFS_with_model_predictions$predictions == "ND")
# get variants predicted as ND. 
LFS_with_model_predictions$Protein_p1_TP53[LFS_with_model_predictions$Set == "prediction" & LFS_with_model_predictions$predictions == "ND"]


# keep table used for analysis.
write.csv(LFS_with_model_predictions, "C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\LFS analysis\\LFS_with_model_predictions_10_18_2020.csv")
