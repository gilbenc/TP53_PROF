library(plyr)
library(ggplot2)
library(tidyr)
library(easyGgplot2)
library(dplyr)
library("RColorBrewer")

# load all_data, fix row names.
input.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\GBM- Test final model"
all_data.path = file.path(input.path, "all.mutations.predictions_with_features.test.functional.csv", fsep = .Platform$file.sep)
all_data<-read.table(all_data.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
row.names(all_data) <- all_data[,1]
all_data <- all_data[,-1]

# remove 41 variants from this analysis.
variants_41 <- read.csv("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\41 set\\41_variants_3_models_predictions.csv")
all_data <- all_data[!(row.names(all_data) %in% variants_41$mutation),]

# put predictions in label. labels > predictions.
all_data$predictions[all_data$label != "null"] <- all_data$label[all_data$label != "null"]  


# order data by prevalence.
all_data_order <- all_data[order(all_data$Records_Number_prot),]

# create df from prevalence, label and set
df <- all_data_order[,c("Records_Number_prot", "predictions", "Set")]
row.names(df) <- row.names(all_data_order)
colnames(df)[1] <- "count"

df$Mutations <- 1:2273
df$count <- log2(df$count+2)
df$Set[df$Set == "Validation"] <- "train"

df$combined_category <- paste(df$predictions, df$Set)
df_train <- df
df_pred <- df 


# send prediction set to negative axis in df.
df$count[df$Set == "prediction"] <- df$count[df$Set == "prediction"]*(-1)

setwd("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\plots_for_paper\\tiff plots from R\\Figure_2_UMD_mds_scatter_freq")

### 1. plot by set.
# plot (high res)
barplot_train_pred_set <- ggplot2.barplot(data=df, xName='Mutations', yName="count", width = 2,
                groupName='Set', groupColors =  c("#8DD3C7", "#FB8072"), aes(x="Mutations")) +  theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
tiff("Barplot_mutations_log2(count+2)_train_pred_flipped_10_22_2020.tiff", units="in", width=15, height=8, res=300)
ggpar(barplot_train_pred_set, ylim = c(-12, 12))
dev.off()

### 2. plot train by D/ND
df_train$count[df_train$combined_category == "ND train"] <- df_train$count[df_train$combined_category == "ND train"]*(-1)

# train
barplot_train_D_ND <- ggplot2.barplot(data=df_train, xName='Mutations', yName= 'count', width = 2, groupName='combined_category', groupColors = c("white", "#80B1D3", "white", "#FFED6F")) +  theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
tiff("Barplot_mutations_train_only_D_flipped_ND_10_22_2020.tiff", units="in", width=15, height=8, res=300)
ggpar(barplot_train_D_ND, ylim = c(0, 12))
dev.off()

### 3. plot pred by D/ND

df_pred$count[df_pred$predictions == "ND"] <- df_pred$count[df_pred$predictions == "ND"]*(-1)

barplot_pred_D_ND <- ggplot2.barplot(data=df_pred, xName='Mutations', yName= 'count', width = 2,
                groupName='combined_category', groupColors = c("#80B1D3", "white", "#FFED6F", "white")) +  theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
tiff("Barplot_mutations_Preds_only_D_flipped_ND_10_22_2020.tiff", units="in", width=15, height=8, res=300)
ggpar(barplot_pred_D_ND, ylim = c(-10, 10))
dev.off()



# see colors. change name to: Set1/Set2/Set3
# look here for more: http://www.sthda.com/english/wiki/colors-in-r#using-rcolorbrewer-palettes
display.brewer.pal(n = 8, name = 'Set1')
brewer.pal(n = 8, name = "Set1")


# t test for pred set: D vs. ND
pred <- df[df$Set == "prediction",]

group_by(pred, predictions) %>%
        summarise(
                count = n(),
                mean = mean(count, na.rm = TRUE),
                sd = sd(count, na.rm = TRUE)
        )
