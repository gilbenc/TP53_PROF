

library(caret)
library(pROC)
library(dplyr)
library(plyr)
library(tidyr)
library(ComplexHeatmap)
library(reshape)
library(ggfortify)
library(ramify)
library(dendextend)
library(factoextra)
library(NbClust)
library(magrittr)
library(ggpubr)
library("RColorBrewer")

### Setup data ###
input.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53"
analysis.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis"

FEATURES_DF_PATH = file.path(input.path, "data_list.csv", fsep = .Platform$file.sep)
df_features = read.table(FEATURES_DF_PATH, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
features = df_features[grepl("Yes",df_features$`as a feature`),]$'Field Name'
features.act<-df_features["Yes (activity)" == df_features$`as a feature`,]$'Field Name'
features.vivo <- df_features["Yes (vivo)" == df_features$'as a feature',]$'Field Name'
features.no.act <-setdiff(features, features.act)
features.no.vivo <- setdiff(features, features.vivo)
features.comp <- setdiff(features.no.act, features.vivo)
features.no.comp <- c(features.act, features.vivo)

#Decides which features to use (comp, act, vivo)
# and whether to include labels in data (includeLabel).
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
        #algorithm's properties (some require labels seperately)
        if (includeLabel==T) {
                cols<-c(cols, "label")
        }
        dat <- dat[,cols]
        
        return(dat)
}

#  creates dat, which contains both pos&neg with their labels.
addLabel<-function(pos,neg) {
        #  pos<-prepare.data(pos)
        #  neg<-prepare.data(neg)
        
        #  make sure that no mutation appears in both neg and pos.
        stopifnot(nrow(merge(pos,neg,by='ID'))==0)
        
        #  create col named label, add functional variable and merge data set
        pos$label<-rep("D", nrow(pos))
        neg$label<-rep("ND", nrow(neg))
        dat<-rbind(pos,neg)
        #  make R treat labels ("D", "ND") as factor levels, not as strings.
        dat$label<-factor(dat$label, levels=c("D","ND"), labels=c("D","ND"))
        return(dat)
}

#create a small table, michal created this for convenience in showing the data.
data.to.desc<-function(name, c_pos, c_neg, removedActivity) {
        return(c(name,sprintf("Dataset %s. Pos: %s, Neg: %s, ActivityRemoved: %s", name, paste(c_pos, collapse = ","), paste(c_neg, collapse = ","), removedActivity)))
}

#input: dataset&path, what column in dataset we choose to be positive. which is negative.
# do we want to remove activity?
#output: dataset divided into pos&neg with labels and chosen features.
prepare.dataset <- function(orig.dat, dataset.path, c_pos, c_neg, remove.Act) {
        if(file.exists(dataset.path)){
                dataset<-read.table(dataset.path, sep=",", header=T, stringsAsFactors = F, check.names = F, fill=F)
        } else {
                
                set.pos<-setNames(data.frame(matrix(ncol = length(names(orig.dat)), nrow = 0)), names(orig.dat)) #create a matrix set.pos with the col names of orig.dat. and no rows.
                for (pos in c_pos)
                {
                        set.pos<-rbind(set.pos,orig.dat[!is.na(orig.dat[,pos]),])
                }
                
                set.neg<-setNames(data.frame(matrix(ncol = length(names(orig.dat)), nrow = 0)), names(orig.dat))
                for (neg in c_neg)
                {
                        set.neg<-rbind(set.neg,orig.dat[!is.na(orig.dat[,neg]),])
                }
                
                dataset <- addLabel(pos=set.pos,neg=set.neg)
                
                #print(Heatmap(log(dataset[,features.act]+1), name = "hm", split = dataset$label,show_row_names = FALSE,                  #show_row_dend = FALSE, show_column_dend = FALSE))
                
                #rownames(dataset1)<-dataset1$cDNA_variant
                dataset <- choose.features(dataset, remove.Act = remove.Act, includeLabel=T)
                
                #write.table(dataset1, file=dataset1.path, quote = F, row.names = F, sep = ',')
        }
        return(dataset)
}
orderColumns=function(x){
        x<-as.numeric(as.character(x)) #first convert each column into numeric if it is from factor
        x[is.na(x)] =median(x, na.rm=TRUE) #convert the item with NA to median value from the column
        x #display the column
}

k_clust <- function(mds, clust_num){
        kmeans(mds, clust_num)$cluster %>%
                as.factor() %>% return()
        
}

visualize_mds <- function(mds, label, label2=NULL) {
        colnames(mds) <- c("Dim.1", "Dim.2")
        mds$label <- label
        mds$label2 <- label2
        if(!is.null(label2)) {
                ggscatter(mds, x = "Dim.1", y = "Dim.2",
                          palette = c("#FFFF33", "#E41A1C", "#377EB8"),
                          color = "label", 
                          shape = "label2",
                          size = 1.5,
                          repel = TRUE) + ggtitle(label = "Multi-Dimension Scaling of TP53 Mutations in UMD")
        }
        else {
                ggscatter(mds, x = "Dim.1", y = "Dim.2",
                          palette = c("#E41A1C", "#4DAF4A", "#377EB8"),
                          color = "label",
                          size = 1.5,
                          repel = TRUE) + ggtitle(label = "Multi-Dimension Scaling of TP53 Mutations in UMD")
                
        }
        
}

#################
##   Main
#################

# Setup Data
alldata.path = file.path(input.path, "all_data.csv", fsep = .Platform$file.sep)
dataset1.path = file.path(input.path, "dataset1.csv", fsep = .Platform$file.sep)
dataset2.path = file.path(input.path, "dataset2.csv", fsep = .Platform$file.sep)

# read orig.dat
orig.dat<-read.table(alldata.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F, na.strings = c("", "-", "No_Data", "No data (synonymous)", "#N/A", "No_Data", " No_Data", "No_Data (synonymous)", "*", "?"),
                     row.names='cDNA_Variant')

# include only missense data
orig.dat<-orig.dat[orig.dat$Variant_Classification=='Missense',]
        
### read all mutations: Train + Prediction ###
all_mutations_PATH = file.path(input.path, "\\data_presentation\\GBM- Test final model\\all.mutations.predictions_with_features.test.functional.csv", fsep = .Platform$file.sep)
all_data <- read.table(all_mutations_PATH, sep=",", header=T, stringsAsFactors = F, check.names = F)
row.names(all_data) <- all_data[,1]
all_data <- all_data[,-1]

# $predictions: predictions > lavel. choose third row over these two for label > predictions.
all_data$predictions[all_data$predictions == "null"] <- all_data$label[all_data$predictions == "null"]  
all_data$label <- all_data$predictions
#all_data$label[all_data$label == "null"] <-all_data$predictions[all_data$label == "null"]



# choose features for all mutations
all_mutations_features <- choose.features(all_data, TRUE, TRUE, TRUE, includeLabel=F)

# create a table with D mutations only. to cluster.
D_mutations <- all_data[all_data$label == "D",]

D_mutations_features <- choose.features(D_mutations, TRUE, TRUE, TRUE, includeLabel=F)

all_mutations_features=data.frame(apply(all_mutations_features,2,orderColumns))
row.names(all_mutations_features) <- row.names(all_data)

# Compute MDS
mds <- all_mutations_features %>%
        dist() %>% # dist's default is "euclidean". 
        cmdscale() %>%
        as_tibble()

#set
all_data$Set[all_data$Set == "Validation"] <- "train"

mds$label <- all_data$label
mds$Set <- all_data$Set

display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n = 12, name = "Set3")

setwd("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\plots_for_paper\\tiff plots from R\\Figure_3_UMD_mds_scatter_freq\\")

tiff("MDS_all_TP53_mutations_by_set_2020_12_09.tiff", units="in", width=15, height=8, res=300)
# plot by set
ggscatter(mds, x = "V1", y = "V2", 
          xlab = FALSE, ylab = FALSE,
          size = 6, repel = TRUE) + font("xy.text", size = 30) + geom_point(aes(fill = factor(Set)), shape = 21, size = 6,colour = "black") + scale_fill_manual(values=c("#8DD3C7", "#FB8072"))

dev.off()

p <- ggplot(mds %>% arrange(label == "D"), aes(V1, V2))
tiff("MDS_all_TP53_mutations_by_label_2020_12_09.tiff", units="in", width=15, height=8, res=300)
# plot by label
ggscatter(mds %>% arrange(label == "D"), x = "V1", y = "V2", 
          xlab = FALSE, ylab = FALSE,
          size = 6, repel = TRUE) + font("xy.text", size = 30) + geom_point(aes(fill = factor(label)), shape = 21, size = 6,colour = "black") + scale_fill_manual(values=c("#80B1D3", "#FFED6F"))
dev.off()
