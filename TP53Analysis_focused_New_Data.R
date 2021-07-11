
#' ---
#' title: "TP53 Analysis"
#' output: csv files
#' ---


library(caret)
library(pROC)
library(dplyr)
library(ggplot2)
library(tidyr)
library(glue)
library(OpenMPController)
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

input.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53"
analysis.path<-"C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis"


# Decides which features to use (comp, act, vivo)
# and whether to include labels in data (includeLabel).
choose.features<-function(dat, comp, act, vivo, includeLabel=F) {
        
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


#  creates dat, which contains both pos&neg with their labels.
addLabel<-function(pos,neg) {

    pos$cDNA_Variant <- row.names(pos)
    neg$cDNA_Variant <- row.names(neg)
    stopifnot(nrow(merge(pos,neg,by='cDNA_Variant'))==0)
   
#  create col named label, add functional variable and merge data set
  pos$label<-rep("D", nrow(pos))
  neg$label<-rep("ND", nrow(neg))
  dat<-rbind(pos,neg)
#  make R treat labels ("D", "ND") as factor levels, not as strings.
  dat$label<-factor(dat$label, levels=c("D","ND"), labels=c("D","ND"))
  return(dat)
}


#input: dataset&path, what column in dataset we choose to be positive. which is negative.
# do we want to remove activity?
#output: dataset divided into pos&neg with labels and chosen features.
prepare.dataset <- function(orig.dat, c_pos, c_neg) {
    
        #create two sub-tables of orig.dat- one with all pos mutations, one with all neg. 
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
    
    #combine pos with neg and add labels.
    dataset <- addLabel(pos=set.pos,neg=set.neg)
    
    #print(Heatmap(log(dataset[,features.act]+1), name = "hm", split = dataset$label,show_row_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE))
    
    #rownames(dataset1)<-dataset1$cDNA_variant
    dataset <- choose.features(dataset, comp = T, act = T, vivo = T, includeLabel=T)
    
    #write.table(dataset1, file=dataset1.path, quote = F, row.names = F, sep = ',')
    dataset[1:(ncol(dataset)-1)]=data.frame(apply(dataset[1:(ncol(dataset)-1)],2,orderColumns))
  return(dataset)
}


run.test<-function(model, dataset, testSet, to.predict, desc, feature_selection, comp, act, vivo) {

  probs <- predict(model,newdata = testSet, na.action = na.pass, type="prob")
  #predicted <- predict(model, newdata= testSet, na.action = na.pass)
  #head(probs)
  ROC <- roc(predictor=probs$D,
                 response=testSet$label,
                 levels=rev(levels(testSet$label)))
  print(paste0("AUC ", ROC$auc))
  

  predicted <- predict.cutoff(model, newdata = testSet, model$cutoff)
  #print(predicted, testSet$label)
  predicted <- factor(predicted, c("D", "ND"))
  print(caret::confusionMatrix(predicted, testSet$label))
  
  
  testSet<-cbind(testSet,predicted)
  testSet<-cbind(testSet,probs$D)
  model$auc <- ROC$auc
  #save predictions for all data and create tables accordingly.
  dataPresentation(dataset, model, testSet, to.predict, feature_selection, comp, act, vivo)
  return(model)

  # write.table(testSet, file=sprintf(file.path(analysis.path, "%s_test_gbm_model.rds", fsep = .Platform$file.sep), desc), quote = F, row.names = F, sep = ',')
  ### To examine later. gbc. ###  
  # visualize results
  #visualize(testSet, file=sprintf(file.path(analysis.path, "%s_heatmap.png", fsep = .Platform$file.sep), desc))
  ##############################
  
  #predicted.test<-predict(model, newdata = to.predict, na.action = na.pass)
  #to.predict$predicted<-predicted.test
  #write.table(to.predict, file="C:\\Hadassah_Shai\\Projects\\TP53\\analysis\\2.1\\pred.dataset1loocv_test.rf.loocv.rds", quote = F, row.names = F, sep = ',')
  
  
  
}

visualize<-function(dat, file) {
  
  png(filename=file)
  print(Heatmap(log(dat[,features.no.act]+1), name = "hm", split = dat$predicted,show_row_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE))
  dev.off()
  
}

#run Gradient Boosting Machine using caret
#input:training set. output: trained model.
run.gbm <- function(trainSet, met, desc) {
        trainMatrix <- as.matrix(select(trainSet, -label)) # requires input as matrix
        trainLabel <- trainSet$label # label is needed outside the matrix.
        
        control <- caret::trainControl(method = "repeatedcv", classProbs = TRUE, summaryFunction = twoClassSummary, verboseIter = FALSE, allowParallel = TRUE, savePredictions = "final")
        
        # Tune Grid:
        # grid <- tune.gbm(trainMatrix, trainLabel, control)
        # print(grid)
        
        ##############
        # Test Grids #
        ##############
        
        ### all_features test
        if(desc[[1]] == 5.1) {
          grid <- expand.grid(nrounds = 26, eta = 0.09495805, max_depth = 5, gamma = 0.2551619, colsample_bytree = 0.7369411, min_child_weight = 1, subsample = 0.8628348)
        }
        ### functional (act_vivo) test
        if(desc[[1]] == 5.4){
          grid <- expand.grid(nrounds = 30, eta = 0.08622785, max_depth = 3, gamma = 0.2804512, colsample_bytree = 0.7805748, min_child_weight = 1, subsample = 0.729023)
        }
        ### comp only test
        if(desc[[1]] == 5.5) {
          grid <- expand.grid(nrounds = 30, eta = 0.09014091, max_depth = 3, gamma = 0.2257494, colsample_bytree = 0.4595327, min_child_weight = 1, subsample = 0.8270071)
        }
        
        ##############
        # Val Grids: #
        ##############
        # 
        # ### all_features val
        if(desc[[1]] == 5.1) {
          grid <- expand.grid(nrounds = 28, eta = 0.09188703, max_depth = 4, gamma = 0.272785, colsample_bytree = 0.4386464, min_child_weight = 1, subsample = 0.8605561)
        }
        # ### functional (act_vivo) val
        # if(desc[[1]] == 5.4) {
        #   grid <- expand.grid(nrounds = 29, eta = 0.08300878, max_depth = 5, gamma = 0.4094955, colsample_bytree = 0.9161638, min_child_weight = 1, subsample = 0.7730019)
        # }
        # ### comp only val
        # if(desc[[1]] == 5.5) {
        #   grid <- expand.grid(nrounds = 26, eta = 0.09014091, max_depth = 1, gamma = 0.2257494, colsample_bytree = 0.4595327, min_child_weight = 1, subsample = 0.8628348)
        # }
        
        set.seed(1234567)
        model <- caret::train(
                x = trainMatrix,
                y = as.factor(trainLabel),
                trControl = control,
                tuneGrid = grid,
                method = "xgbTree",
                verbose = TRUE,
                metric = met
        )
        return(model)
}

# randomly select parameters for tunning, select best ones.
tune.gbm<- function(trainMatrix, trainLabel, control) {
        
        ### 1. tune: eta###
        # nrounds not under 200.
        tune_grid <- expand.grid(
                nrounds = c(sample(200:360, 1), sample(360:520, 1), sample(520:680, 1), 
                            sample(680:840, 1), sample(840:1000, 1)),
                eta = c(runif(1, 0, 0.025), runif(1, 0.025, 0.05), runif(1, 0.05, 0.1), 
                        runif(1, 0.1, 0.3)),
                max_depth = c(2, 3, 4, 5),
                gamma = 0,
                colsample_bytree = 1,
                min_child_weight = 1,
                subsample = 1
        )
        
        
        xgb_tune <- caret::train(
                x = trainMatrix,
                y = trainLabel,
                trControl = control,
                tuneGrid = tune_grid,
                method = "xgbTree",
                verbose = TRUE,
                metric = met
        )
        
        # plot during tunning. 
        # tuneplot <- function(x, probs = .90) {
        #         ggplot(x) +
        #                 coord_cartesian(ylim = c(quantile(x$results$RMSE, probs = probs), min(x$results$RMSE))) +
        #                 theme_bw()
        # }
        
        # tuneplot(xgb_tune)
        
        ### 2. tune: max_depth, min_child_weight ###
        tune_grid2 <- expand.grid(
                nrounds = c(sample(200:360, 1), sample(360:520, 1), sample(520:680, 1), 
                            sample(680:840, 1), sample(840:1000, 1)),
                eta = xgb_tune$bestTune$eta,
                max_depth = c(3, 4, 5),
                gamma = 0,
                colsample_bytree = 1,
                min_child_weight = c(1, 2, 3),
                subsample = 1
        )
        
        xgb_tune2 <- caret::train(
                x = trainMatrix,
                y = trainLabel,
                trControl = control,
                tuneGrid = tune_grid2,
                method = "xgbTree",
                verbose = TRUE,
                metric = met
        )
        
        #tuneplot(xgb_tune2)
        
        ### 3. tune: colsample_bytree, subsample ###
        tune_grid3 <- expand.grid(
                nrounds = c(sample(200:360, 1), sample(360:520, 1), sample(520:680, 1), 
                            sample(680:840, 1), sample(840:1000, 1)),
                eta = xgb_tune$bestTune$eta,
                max_depth = xgb_tune2$bestTune$max_depth,
                gamma = 0,
                colsample_bytree = runif(4, 0.4, 1),
                min_child_weight = xgb_tune2$bestTune$min_child_weight,
                subsample = runif(3, 0.5, 1.0)
        )
        
        xgb_tune3 <- caret::train(
                x = trainMatrix,
                y = trainLabel,
                trControl = control,
                tuneGrid = tune_grid3,
                method = "xgbTree",
                verbose = TRUE,
                metric = met
        )
        
        # tuneplot(xgb_tune3, probs = .95)
        
        ### 4. tune: gamma
        tune_grid4 <- expand.grid(
                nrounds = c(sample(200:360, 1), sample(360:520, 1), sample(520:680, 1), 
                            sample(680:840, 1), sample(840:1000, 1)),
                eta = xgb_tune$bestTune$eta,
                max_depth = xgb_tune2$bestTune$max_depth,
                gamma = runif(7, 0, 1.0),
                colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
                min_child_weight = xgb_tune2$bestTune$min_child_weight,
                subsample = xgb_tune3$bestTune$subsample
        )
        
        xgb_tune4 <- caret::train(
                x = trainMatrix,
                y = trainLabel,
                trControl = control,
                tuneGrid = tune_grid4,
                method = "xgbTree",
                verbose = TRUE,
                metric = met
        )
        
        # tuneplot(xgb_tune4)
        
        #step 5: retune eta, tune nrounds
        tune_grid5 <- expand.grid(
                nrounds = sample(30, 100, 3000),
                eta = c(runif(1, 0, 0.015), runif(1, 0.015, 0.025), runif(1, 0.025, 0.05), 
                        runif(1, 0.05, 0.1)),
                max_depth = xgb_tune2$bestTune$max_depth,
                gamma = xgb_tune4$bestTune$gamma,
                colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
                min_child_weight = xgb_tune2$bestTune$min_child_weight,
                subsample = xgb_tune3$bestTune$subsample
        )
        
        xgb_tune5 <- caret::train(
                x = trainMatrix,
                y = trainLabel,
                trControl = control,
                tuneGrid = tune_grid5,
                method = "xgbTree",
                verbose = TRUE,
                metric = met
        )
        
        # tuneplot(xgb_tune5)
        
        final_grid <- expand.grid(
                nrounds = xgb_tune5$bestTune$nrounds,
                eta = xgb_tune5$bestTune$eta,
                max_depth = xgb_tune5$bestTune$max_depth,
                gamma = xgb_tune5$bestTune$gamma,
                colsample_bytree = xgb_tune5$bestTune$colsample_bytree,
                min_child_weight = xgb_tune5$bestTune$min_child_weight,
                subsample = xgb_tune5$bestTune$subsample
        )
        
        return(final_grid)
}

#run Randeom Forests. input: training set. output: trained model.
run.rf <- function(trainSet, met, tuneLength) {
        ### Random Forests ###
        # TODO: change method to LOOCV (leave one out. means for every iteration ALL data is training, one example is cv.) Remember: check for new parameters when changing method.
        # ctrl contains more 'hyperparamters' for the model. 
        #input: repeatedcv means divide into training&cv 'repeats' amount of times. classprobs means we want the model to return a linear value rather than a binary value. summart function- summary is shown in a 2X2 table. 
        ctrl <- trainControl(method = "repeatedcv",number = 5, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = "final")
        
        # Random Forest (rf) on train set.
        # use all features to predict label. when na: use median impute. 
        set.seed(1234567)
        model <- train(label~., 
                       data=trainSet,
                       method='rf',
                       trControl=ctrl,
                       tuneLength=tuneLength, 
                       na.action = na.pass, 
                       preProcess = "medianImpute",
                       metric=met)  
        return(model)
} 

#input: dataset, desc, tuneLength= do we want the model to tune hyperparameters? how many values should the model attempt before choosing? met=ROC- choose according to best ROC.
run.ml<-function(trainSet, CVSet, desc, to.predict, tuneLength=10, met='ROC', feature_selection, comp, act, vivo, rf = F) {
        
  ## temp code to run line by line ###
  # trainSet_temp <- trainSet
  # CVSet_temp <- CVSet
  # trainCVSet_temp <- trainCVSet
  # testSet_temp <- testSet
  # #
  # trainSet <- trainSet_temp
  # CVSet <- CVSet_temp
  # trainCVSet <- trainCVSet_temp
  # testSet <- testSet_temp
  # 

  # trainSet <- trainCVSet5
  #  CVSet <- testSet5
  # desc <- desc5.5[[1]]
  # to.predict <- orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),]
  # tuneLength = 10
  # met = 'ROC'
  # feature_selection = "comp"
  # comp = T
  # act = F
  # vivo = F

  ### Choose model type between GBM/RF ###
 
  if(rf) {
    model <- run.rf(trainSet, met, tuneLength)
  }
  else {
    model <- run.gbm(trainSet, met, desc)
  }
  print(model)      

   
  #calculate cutoff.
   
  # calculate and print some possible cutoffs  for model##
  # calcsomelevs = function(x) {calc.cutoff(model$pred$D, labels=model$pred$obs, x)}
  # sapply(c(0.94, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99), calcsomelevs)
  # sapply(c(0.95, 0.952, 0.954, 0.956, 0.958, 0.96), calcsomelevs)
  
  model$cutoff<-calc.cutoff(model$pred$D, labels=model$pred$obs, 0.954)
  
   
   filename = sprintf(file.path(input.path, "model_%s_%s.gbm.test.repeatedcv.rds", fsep = .Platform$file.sep), desc, format(Sys.time(), "%d-%b-%Y_%H.%M"))
   saveRDS(model, file=filename)
   print(varImp(model, scale=F))
   
   
   #run model for test set.
   model <- run.test(model, trainSet, CVSet, to.predict, desc, feature_selection, comp = comp, act = act, vivo = vivo)
  #########################################
  
  
  # # SVM, Radial
  # model <- train(label~.,
  #                   data=trainSet,
  #                   method = "svmRadial",   # Radial kernel
  #                   trControl=ctrl,
  #                   tuneLength = tuneLength,					# 10 values of the cost function
  #                   na.action = na.pass,
  #                   preProc = c("center","scale","medianImpute"),  # Center and scale data
  #                   metric=metric)
  # print(model)
  # run.test(model, testSet)
  # 
  # 
  #print(table(predicted, testSet$label))
  #print(confusionMatrix(predicted, testSet$label))
  return(model)
}

#create a small table, michal created this for convenience in showing the data.
data.to.desc<-function(name, c_pos, c_neg, removedActivity) {
  return(c(name,sprintf("Dataset %s. Pos: %s, Neg: %s, ActivityRemoved: %s", name, paste(c_pos, collapse = ","), paste(c_neg, collapse = ","), removedActivity)))
}

#predict.new: this function is for using the model to predict mutations that are in dispute, after having the final data.
#df- all dataset. to.predict- take only rows that are in dispute. model- the learned model.
predict.new<-function(name, df, to.predict, model) {
  predicted.col<-sprintf("predicted_%s", name)
  df[rownames(to.predict), predicted.col]<-factor(predict.cutoff(model, newdata = df[rownames(to.predict),], cutoff = model$cutoff), levels=model$levels)
  
  return(df)
}

checkAgreement<-function(df,pred_names) {
  for (pred_name in pred_names) {
    # SIFT
    col<-sprintf("agree_sift_%s", pred_name)
    df[,col]<-ifelse(((df[,pred_name]=="D") & (df[,"Sift_Prediction"]=="Damaging")) | 
                                                       ((df[,pred_name]=="ND") & (df[,"Sift_Prediction"]=="Tolerated")),T,F)
    print("SIFT agreement")
    print(table(df[,col]))
    
    # POLYPHEN
    col<-sprintf("agree_polyphen_HumVar_%s", pred_name)
    df[,col]<-ifelse(((df[,pred_name]=="D") & (grepl("damaging", df[,"PolyphenNo data (synonymous)2_HumVar"]))) | 
                                                       ((df[,pred_name]=="ND") & (df[,"PolyphenNo data (synonymous)2_HumVar"]=="Benign")),T,F)
    print("POLYPHEN agreement")
    print(table(df[,col]))
  }
  return(df)
  
}

analyze.diff<-function(dat, pred_col, name) {
  dat.scale<-as.data.frame(scale(dat[,names(dat)!=pred_col]))
  dat.scale[,pred_col]<-dat[pred_col]
  dat.melt<-melt(dat.scale, id=c(pred_col))
  filename<-file.path(analysis.path, sprintf("%s_boxplot.png", name), fsep = .Platform$file.sep)
  png(filename=filename)
  print(ggplot(dat.melt, aes_string(x="variable", y="value", fill=pred_col)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  dev.off()
  return(dat.scale)
}

topVars<-function(model, n=20) {
  vars<-rownames(head(varImp(model, scale=F)$importance,n))
  vars<-unlist(lapply(vars, function(x) gsub('`', "", x)))
  return(vars)
}

plot.mds<-function(dat, feat, col, name) {
  d <- dist(scale(dat[,feat])) # euclidean distances between the rows
  fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  filename<-file.path(analysis.path, sprintf("%s_mds.png", name), fsep = .Platform$file.sep)
  png(filename=filename)
  print(autoplot(fit, colour = dat[,col], alpha=0.4))
  dev.off()
  
}

#make sure sens.level works when new data arrives.
#calculate cutoff. we want sensitivity to be sens.level.
#algorithm checks for closest value and uses that value.
calc.cutoff<-function(predicted, labels, sens.level=0.95) {
        #predicted = model$pred$D 
        #labels= model$pred$obs
  ROC<-roc(predictor=predicted,
           response=labels,
           levels=c("ND", "D"))
  print(ROC)
  plot(ROC)
  # we want the sensitivity level to be one of the roc$sensitivities, otherwise threshold is NA
  #check distance between sensitivities and sens.level. choose minimum value's index. use that value as new.sens.
  new.sens<-ROC$sensitivities[[which.min(unlist(lapply(ROC$sensitivities, FUN = function(x) abs(x-sens.level))))]]
  #take new.sens's values (threshold, specificity) into res.
  res<-coords(ROC, new.sens, input= c("sensitivity"), ret=c("threshold", "specificity", "sensitivity"), transpose = TRUE)
  
  #choose res's threshold as cutoff.
  cutoff<-res[["threshold"]]
  cat("Original sensitivity level: ", sens.level, ", New sensitivity level: ", new.sens, ", specificity is ", res[["specificity"]], ", cutoff is ", cutoff, "\n")
  
  
  #pred<-predicted[labels=="D",]
  #pred<-pred[order(pred$D),]
  #cutoff<-pred[floor(nrow(pred)*(1-sens.level)),"D"]
  
  return(cutoff)
}


predict.cutoff<-function(model, newdata, cutoff) {
        pred<-predict(model, newdata = newdata, na.action = na.pass, type="prob")
  return(ifelse(pred$D>=cutoff, "D", "ND") )
}

# eliminate NAs, also make sure all columns are numeric and not factors.
orderColumns=function(x){
        x<-as.numeric(as.character(x)) #first convert each column into numeric if it is from factor
        x[is.na(x)] =median(x, na.rm=TRUE) #convert the item with NA to median value from the column
        x #display the column
}

# runs model for data called 'preds', which is all the data in orig.dat that is not in dataset.
# create tables of train, test and preds, with labels, predictions and features. 
# also creates a table of all false predictions; where label != prediction in Test set.
dataPresentation=function(dataset, model, testSet, to.predict, feature_selection, comp, act, vivo){
        
  # Use model on pred data (all mutations that are not included in dataset5's train/CV/test sets).
        predsSet<-choose.features(dat=to.predict, comp = comp, act = act, vivo = vivo)

        predsSet[1:ncol(predsSet)]=data.frame(apply(predsSet[1:ncol(predsSet)],2,orderColumns))

        #predicted <- predict(model, newdata= testSet, na.action = na.pass)
        #head(probs)
        preds <- predict(model,newdata = predsSet, na.action = na.pass, type="prob")
        predicted <- predict.cutoff(model, newdata = predsSet, model$cutoff)
        predicted <- factor(predicted, c("D", "ND"))
        predsSet<-cbind(predsSet,predicted, preds$D)

        # data presentation. mutation. test/train/pred, label (for train & test) ,predictions (for test & pred)
        #train.data
        model$pred$obs <- model$pred$obs[order(model$pred$rowIndex)]
        model$pred$D <- model$pred$D[order(model$pred$rowIndex)]
        train_mutations = row.names(model$trainingData)
        train_AA_change <- orig.dat[train_mutations,]$Protein_p1_TP53
        train.data <- data.frame(train_mutations, train_AA_change, "train", model$pred$obs, "null", model$pred$D)
        colnames(train.data) <- c("mutation", "protein_change", "Set", "label", "predictions", "score")
        #test.data
        test_mutations = row.names(testSet)
        test_AA_change <- orig.dat.protein[test_mutations,]$Protein_p1_TP53
        test.data <- data.frame(test_mutations, test_AA_change, "Validation", testSet$label, testSet$predicted, testSet$`probs$D`)
        colnames(test.data) <- c("mutation", "protein_change", "Set", "label", "predictions", "score")
        #pred.data
        pred_mutations <- row.names(predsSet)
        pred_AA_change <- orig.dat.protein[pred_mutations,]$Protein_p1_TP53
        pred.data <- data.frame(pred_mutations, pred_AA_change, "prediction", "null", predsSet$predicted, predsSet$`preds$D`)
        colnames(pred.data) <- c("mutation", "protein_change", "Set", "label", "predictions", "score")

        final.data <- rbind(train.data, test.data, pred.data)
        # write.csv(final.data, file = paste0("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\all.mutations.predictions.validation.", feature_selection, ".csv"))

        temp_orig.dat.protein <- orig.dat.protein
        temp_orig.dat.protein$mutation <- row.names(temp_orig.dat.protein)
        features_final.data <- merge(temp_orig.dat.protein, final.data, by = "mutation")
        row.names(features_final.data) <- features_final.data$mutation
        features_final.data <- features_final.data[,-1]
        # write.csv(features_final.data, file = paste0("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\all.mutations.predictions_with_features.validation.", feature_selection, ".csv"))

        false.predictions <- test.data
        false.predictions <- false.predictions[false.predictions$label!=false.predictions$predictions,]
        # write.csv(false.predictions, file = paste0("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\data_presentation\\false.predictions.validation.", feature_selection, ".csv"))
}

#################
##   Main
#################

FEATURES_DF_PATH = file.path(input.path, "data_list.csv", fsep = .Platform$file.sep)
df_features = read.table(FEATURES_DF_PATH, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
features = df_features[grepl("Yes",df_features$`as a feature`),]$'Field Name'
features.act<-df_features["Yes (activity)" == df_features$`as a feature`,]$'Field Name'
features.vivo <- df_features["Yes (vivo)" == df_features$'as a feature',]$'Field Name'
features.no.act <-setdiff(features, features.act)
features.no.vivo <- setdiff(features, features.vivo)
features.comp <- setdiff(features.no.act, features.vivo)
features.no.comp <- c(features.act, features.vivo)


# read data files
alldata.path = file.path(input.path, "all_data.csv", fsep = .Platform$file.sep)

#load data into orig.dat.
orig.dat<-read.table(alldata.path, sep=",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F, na.strings = c("", "-", "No_Data", "No data (synonymous)", "#N/A", "No_Data", " No_Data", "No_Data (synonymous)", "*", "?"),
                             row.names='cDNA_Variant')
# include only missense data
orig.dat<-orig.dat[orig.dat$Variant_Classification=='Missense',]


# remove duplications in protein change
orig.dat.protein <- orig.dat[!duplicated(orig.dat$Protein_p1_TP53),]

# create new negative set based on protein record number:
sum_rec_num <- function (prot_change) {
  return(sum(as.numeric(orig.dat$Records_Number[orig.dat$Protein_p1_TP53 == prot_change])))
}
orig.dat.protein$Records_Number_prot <- sapply(orig.dat.protein$Protein_p1_TP53, sum_rec_num)
orig.dat.protein$set_neg <- NA
orig.dat.protein$set_neg[orig.dat.protein$Records_Number_prot <= 1] <- row.names(orig.dat.protein)[orig.dat.protein$Records_Number_prot <= 1]

# read and remove 41 samples thierry analyzed
remove.path = file.path(input.path, "41_removed_variants_by_thierry_24_11_19.csv", fsep = .Platform$file.sep)
samples_to_remove <- read.table(remove.path, sep = ",", fileEncoding="UTF-8-BOM", header=T, stringsAsFactors = F, check.names = F)
orig.dat_removed <- orig.dat.protein[!(row.names(orig.dat.protein) %in% samples_to_remove$mutation),]


# positive: set6. negative: set14 + set15
c_pos<-c("set_6")
# c_neg<-c("set_9", "set_10")
c_neg<-c("set_neg")

# NOTE that if you change dataset5.1, the tables saved in data_presentation will be CHANGED.
desc5.1<-data.to.desc("5.1",c_pos, c_neg, "all_features")
dataset5.1<-prepare.dataset(orig.dat = orig.dat_removed, c_pos = c_pos, c_neg=c_neg)

desc5.2<-data.to.desc("5.2",c_pos, c_neg, "comp_act")
desc5.3<-data.to.desc("5.3",c_pos, c_neg, "comp_vivo")
desc5.4<-data.to.desc("5.4",c_pos, c_neg, "act_vivo")
desc5.5<- data.to.desc("5.5",c_pos, c_neg, "comp")
desc5.6<- data.to.desc("5.6",c_pos, c_neg, "act")
desc5.7<- data.to.desc("5.7",c_pos, c_neg, "vivo")

## create train, validation and test sets. ##
dataset <- dataset5.1
set.seed(123456)
index <- createDataPartition(y=dataset$label, p=0.8, list=FALSE) #divide dataset. keep same amounts of y=dataset$labels on every division.
trainCVSet <- dataset[ index,]
testSet <- dataset[-index,]

# #create Cross Validation, 25% of train set.
set.seed(123456)
index <- createDataPartition(y=trainCVSet$label, p=0.75, list=FALSE) #divide dataset. keep same amounts of y=dataset$labels on every division.
CVSet <- trainCVSet[-index,]
trainSet <- trainCVSet[ index,]


# prepare for 5.2: comp and act
trainSet2 <- choose.features(dat = trainSet, comp = TRUE, act = TRUE, vivo = F, includeLabel = TRUE)
CVSet2 <- choose.features(dat = CVSet, comp = TRUE, act = TRUE, vivo = F, includeLabel = TRUE)
trainCVSet2 <- choose.features(dat = trainCVSet, comp = TRUE, act = TRUE, vivo = F, includeLabel = TRUE)
testSet2 <- choose.features(dat = testSet, comp = TRUE, act = TRUE, vivo = F, includeLabel = TRUE)

# prepare for 5.3: comp and vivo
trainSet3 <- choose.features(dat = trainSet, comp = TRUE, act = F, vivo = T, includeLabel = TRUE)
CVSet3 <- choose.features(dat = CVSet, comp = TRUE, act = F, vivo = T, includeLabel = TRUE)
trainCVSet3 <- choose.features(dat = trainCVSet, comp = TRUE, act = F, vivo = T, includeLabel = TRUE)
testSet3 <- choose.features(dat = testSet, comp = TRUE, act = F, vivo = T, includeLabel = TRUE)

# prepare for 5.4: act and vivo
trainSet4 <- choose.features(dat = trainSet, comp = F, act = TRUE, vivo = T, includeLabel = TRUE)
CVSet4 <- choose.features(dat = CVSet, comp = F, act = TRUE, vivo = T, includeLabel = TRUE)
trainCVSet4 <- choose.features(dat = trainCVSet, comp = F, act = TRUE, vivo = T, includeLabel = TRUE)
testSet4 <- choose.features(dat = testSet, comp = F, act = TRUE, vivo = T, includeLabel = TRUE)

# prepare for 5.5: comp. only.
trainSet5 <- choose.features(dat = trainSet, comp = T, act = F, vivo = F, includeLabel = TRUE)
CVSet5 <- choose.features(dat = CVSet, comp = T, act = F, vivo = F, includeLabel = TRUE)
trainCVSet5 <- choose.features(dat = trainCVSet, comp = T, act = F, vivo = F, includeLabel = TRUE)
testSet5 <- choose.features(dat = testSet, comp = T, act = F, vivo = F, includeLabel = TRUE)

# prepare for 5.6: act only.
trainSet6 <- choose.features(dat = trainSet, comp = F, act = T, vivo = F, includeLabel = TRUE)
CVSet6 <- choose.features(dat = CVSet, comp = F, act = T, vivo = F, includeLabel = TRUE)
trainCVSet6 <- choose.features(dat = trainCVSet, comp = F, act = T, vivo = F, includeLabel = TRUE)
testSet6 <- choose.features(dat = testSet, comp = F, act = T, vivo = F, includeLabel = TRUE)

# prepare for 5.7: vivo only.
trainSet7 <- choose.features(dat = trainSet, comp = F, act = F, vivo = T, includeLabel = TRUE)
CVSet7 <- choose.features(dat = CVSet, comp = F, act = F, vivo = T, includeLabel = TRUE)
trainCVSet7 <- choose.features(dat = trainCVSet, comp = F, act = F, vivo = T, includeLabel = TRUE)
testSet7 <- choose.features(dat = testSet, comp = F, act = F, vivo = T, includeLabel = TRUE)

# run gbm
# run 5.1
cat(sprintf("----  %s  ------",desc5.1[[2]]))
gbm_model5.1 <- run.ml(trainSet, CVSet, desc5.1[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "all_features", comp = T, act = T, vivo = T)

# run 5.2
cat(sprintf("----  %s ------",desc5.2[[2]]))
gbm_model5.2 <- run.ml(trainSet2, CVSet2, desc5.2[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "comp_act", comp = T, act = T, vivo = F)

# run 5.3
cat(sprintf("----  %s ------",desc5.3[[2]]))
gbm_model5.3 <- run.ml(trainSet3, CVSet3, desc5.3[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "comp_vivo", comp = T, act = F, vivo = T)

# run 5.4
cat(sprintf("----  %s ------",desc5.4[[2]]))
gbm_model5.4 <- run.ml(trainSet4, CVSet4, desc5.4[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "act_vivo", comp = F, act = T, vivo = T)

# run 5.5
cat(sprintf("----  %s ------",desc5.5[[2]]))
gbm_model5.5 <- run.ml(trainSet5, CVSet5, desc5.5[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "comp", comp = T, act = F, vivo = F)

# run 5.6
cat(sprintf("----  %s ------",desc5.6[[2]]))
gbm_model5.6 <- run.ml(trainSet6, CVSet6, desc5.6[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "act", comp = F, act = T, vivo = F)

# run 5.7
cat(sprintf("----  %s ------",desc5.7[[2]]))
gbm_model5.7 <- run.ml(trainSet7, CVSet7, desc5.7[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "vivo", comp = F, act = F, vivo = T)


# run RF
# run 5.1
cat(sprintf("----  %s  ------",desc5.1[[2]]))
RF_model5.1 <- run.ml(trainSet, CVSet, desc5.1[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "all_features", comp = T, act = T, vivo = T, rf = T)

# run 5.2
cat(sprintf("----  %s ------",desc5.2[[2]]))
RF_model5.2 <- run.ml(trainSet2, CVSet2, desc5.2[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "comp_act", comp = T, act = T, vivo = F, rf = T)

# run 5.3
cat(sprintf("----  %s ------",desc5.3[[2]]))
RF_model5.3 <- run.ml(trainSet3, CVSet3, desc5.3[[1]], to.predict=orig.dat[!(rownames(orig.dat) %in% rownames(dataset5.1)),], feature_selection = "comp_vivo", comp = T, act = F, vivo = T, rf =T)

# run 5.4
cat(sprintf("----  %s ------",desc5.4[[2]]))
RF_model5.4 <- run.ml(trainSet4, CVSet4, desc5.4[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "act_vivo", comp = F, act = T, vivo = T, rf = T)

cat(sprintf("----  %s  ------",desc5.5[[2]]))
RF_model5.5 <- run.ml(trainSet5, CVSet5, desc5.5[[1]], to.predict = orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "comp", comp = T, act = F, vivo = F, rf = T)

# GBM On Test Set
# run 5.4
cat(sprintf("----  %s  ------",desc5.4[[2]]))
gbm_test_model5.4 <- run.ml(trainCVSet4, testSet4, desc5.4[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "act_vivo", comp = F, act = T, vivo = T)

## Re-calculate and Plot ROC
# calculate
probs <- predict(gbm_test_model5.4,newdata = testSet4, na.action = na.pass, type="prob")
ROC <- roc(predictor=probs$D,
           response=testSet4$label,
           levels=rev(levels(testSet4$label)))
print(paste0("AUC ", ROC$auc))
# plot
setwd("C:\\Users\\USER\\Desktop\\gil\\lab\\TP53\\analysis\\plots_for_paper\\tiff plots from R")
tiff("ROC_curve_Functional_features_test_model_05_04_2020.tiff", units="in", width=15, height=8, res=300)
plot(ROC,main=model$method)
dev.off()

# run 5.1
cat(sprintf("----  %s  ------",desc5.1[[2]]))
gbm_test_model5.1 <- run.ml(trainCVSet, testSet, desc5.1[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "all_features", comp = T, act = T, vivo = T)

# run 5.5
cat(sprintf("----  %s  ------",desc5.5[[2]]))
gbm_test_model5.5 <- run.ml(trainCVSet5, testSet5, desc5.5[[1]], to.predict=orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "comp", comp = T, act = F, vivo = F)


# 10 RUNS COMPARISON
# between GBM and RF, 4 feature options, on validation set.
### IMPORTANT: ### 
### REMOVE SET.SEED before running this chunk of code. 
### also make sure Tune.GBM is applied.
gbm1_AUC <- c()
for(i in 1:10) {
  gbm_model5.1 <- run.ml(trainSet, CVSet, desc5.1[[1]], to.predict=orig.dat[!(rownames(orig.dat) %in% rownames(dataset5.1)),], feature_selection = "all_features", comp = T, act = T, vivo = T)
  gbm1_AUC <- c(gbm1_AUC, gbm_model5.1$auc)
}
gbm2_AUC <- c()
for(i in 1:10) {
  gbm_model5.2 <- run.ml(trainSet2, CVSet2, desc5.2[[1]], to.predict=orig.dat[!(rownames(orig.dat) %in% rownames(dataset5.1)),], feature_selection = "comp_act", comp = T, act = T, vivo = F)
  gbm2_AUC <- c(gbm2_AUC, gbm_model5.2$auc)
}
gbm3_AUC <- c()
for(i in 1:10) {
  gbm_model5.3 <- run.ml(trainSet3, CVSet3, desc5.3[[1]], to.predict=orig.dat[!(rownames(orig.dat) %in% rownames(dataset5.1)),], feature_selection = "comp_vivo", comp = T, act = F, vivo = T)
  gbm3_AUC <- c(gbm3_AUC, gbm_model5.3$auc)
}
gbm4_AUC <- c()
for(i in 1:10) {
  gbm_model5.4 <- run.ml(trainSet4, CVSet4, desc5.4[[1]], to.predict=orig.dat[!(rownames(orig.dat) %in% rownames(dataset5.1)),], feature_selection = "act_vivo", comp = F, act = T, vivo = T)
  gbm4_AUC <- c(gbm4_AUC, gbm_model5.4$auc)
}

gbm5_AUC <- c()
for(i in 1:10) {
  gbm_model5.5 <- run.ml(trainSet5, CVSet5, desc5.5[[1]], to.predict=orig.dat[!(rownames(orig.dat) %in% rownames(dataset5.1)),], feature_selection = "comp", comp = T, act = F, vivo = F)
  gbm5_AUC <- c(gbm5_AUC, gbm_model5.5$auc)
}

RF1_AUC <- c()
for(i in 1:10) {
  RF_model5.1 <- run.ml(trainSet, CVSet, desc5.1[[1]], to.predict=orig.dat[!(rownames(orig.dat) %in% rownames(dataset5.1)),], feature_selection = "all_features", comp = T, act = T, vivo = T, rf = T)
  RF1_AUC <- c(RF1_AUC, RF_model5.1$auc)
}

RF4_AUC <- c()
for(i in 1:10) {
  RF_model5.4 <- run.ml(trainSet4, CVSet4, desc5.4[[1]], to.predict=orig.dat[!(rownames(orig.dat) %in% rownames(dataset5.1)),], feature_selection = "act_vivo", comp = F, act = T, vivo = T, rf = T)
  RF4_AUC <- c(RF4_AUC, RF_model5.4$auc)
}

RF5_AUC <- c()
for(i in 1:10) {
  RF_model5.5 <- run.ml(trainSet5, CVSet5, desc5.5[[1]], to.predict = orig.dat.protein[!(rownames(orig.dat.protein) %in% rownames(dataset5.1)),], feature_selection = "comp", comp = T, act = F, vivo = F, rf = T)
  RF5_AUC <- c(RF5_AUC, RF_model5.5$auc)
}

print("gbm AUC:")
print(gbm1_AUC)
print(paste0("mean: ", mean(gbm1_AUC)))

print("gbm2 AUC:")
print(gbm2_AUC)
print(paste0("mean: ", mean(gbm2_AUC)))

print("gbm3 AUC:")
print(gbm3_AUC)
print(paste0("mean: ", mean(gbm3_AUC)))

print("gbm4 AUC:")
print(gbm4_AUC)
print(paste0("mean: ", mean(gbm4_AUC)))

print("RF1 AUC:")
print(RF1_AUC)
print(paste0("mean: ", mean(RF1_AUC)))

print("RF2 AUC:")
print(RF2_AUC)
print(paste0("mean: ", mean(RF2_AUC)))

print("RF3 AUC:")
print(RF3_AUC)
print(paste0("mean: ", mean(RF3_AUC)))

print("RF4 AUC:")
print(RF4_AUC)
print(paste0("mean: ", mean(RF4_AUC)))

print("RF5 AUC:")
print(RF5_AUC)
print(paste0("mean: ", mean(RF5_AUC)))

save.image("C:/Users/USER/Desktop/gil/lab/TP53/analysis/10_runs_gbm_RF_8_models.RData")
