###Assignment 1###
#Quesion 4

rm(list=ls())
setwd("C:\\Users\\sakshi Mahajan\\Desktop\\hmds_assi1")
getwd()
###load required packages###
library(vegan)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(randomForest)
library(e1071)
library(caret)
####upload given data files### 
clr_species_abundance<-read.table("Assignment1_ClrTrans_Species.txt",sep="\t",header=TRUE)
metadata_species<-read.table("Assignment1_Metadata.txt",sep="\t",header=TRUE)

##########################################################################################################################################
#the classification performance on the three different groups.
#1. Random forest

metadata_species$WHO_severity <- as.numeric(factor(metadata_species$WHO_severity, levels=c("mild", "moderate", "critical_severe")))
##merge clr abundance with metadata
merge_data<-merge(clr_species_abundance,metadata_species,by="Sample_ID")
#convert WHO_servity as a factor
merge_data$WHO_severity <- as.factor(merge_data$WHO_severity)

###for to generate random number###
set.seed(123)

#This randomly samples 64% of the rows from the merged dataset to be used as the training dataset. The remaining 36% of the rows will be used as the test dataset.
train_idx <- sample(nrow(merge_data), 0.85* nrow(merge_data))
#creates the training dataset
train_data <- merge_data[train_idx, ]  
# remove Sample_ID from training data
train_data$Sample_ID<-NULL  
#creates the testing dataset
test_data <- merge_data[-train_idx, ]
# remove Sample_ID from training data
test_data$Sample_ID<-NULL

###builds a random forest model using the training dataset
model <- randomForest(WHO_severity ~ ., data = train_data, importance = TRUE, ntree = 1000,type="classification")

#generates predictions for the test dataset using the random forest model 
predictions <- predict(model, newdata = test_data)

# creates a confusion matrix to evaluate the performance of the model on the test dataset
confusion_matrix <- table(predictions, test_data$WHO_severity)
confusion_matrix

#calculates the overall accuracy of the Random forest model
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)  #calculate accuracy
print(paste0("Accuracy: ", round(accuracy, 3)))
###Create feature importance plot
varImpPlot(model)
importance <- importance(model,type=2)
importance

##plot top 50 features of model
varImpPlot(model, top=50) 

###plot Confusion Matrix of Random Forest model for all 3 grp using ggplot###
con<-confusionMatrix(predictions, test_data$WHO_severity)
con$table
ggplot(data = data.frame(con$table, class = rownames(con$table)),
       aes(x = Prediction, y = Reference, fill = Freq)) + 
  geom_tile() + 
  geom_text(aes(label = Freq), size = 10, colour = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  xlab("Predicted Class") + 
  ylab("True Class") + 
  ggtitle("Confusion Matrix_RF1") + 
  theme(plot.title = element_text(hjust = 0.5))

# Sort the features by their importance score and Select the top 50 most important features
sorted_importance <- row.names(importance)[order(importance[,1], decreasing = TRUE)[1:50]]
sorted_importance
top50_features_RF1 <- names(merge_data[,sorted_importance])
print(top50_features_RF1)

#download top50_features_RF1 file
write.table(top50_features_RF1,"top50_features_RF_1.txt",sep="\t")

#################################################################################################################################################
#the classification performance on the three different groups.
#2. Support Vector Machines(SVM)

#creates a character vector
non_numeric_vars <- c("Age", "Sex")

###for categorical data analysis###
#converts each variable in non_numeric_vars to a factor variable for perform categorical data analysis on these variables.
for (var in non_numeric_vars) {
  merge_data[[var]] <- as.factor(merge_data[[var]])
}
for (var in names(merge_data)) {
  if (is.factor(merge_data[[var]]) && length(levels(merge_data[[var]])) < 2) {
    merge_data[[var]] <- NULL
  }
}
#loop converts any character or factor variable in merged_data_sub to a factor variable. 
for (var in names(merge_data)) {
  if (is.character(merge_data[[var]]) || is.factor(merge_data[[var]])) {
    merge_data[[var]] <- factor(merge_data[[var]])
  }
}
set.seed(123)

#This randomly samples 81% of the rows from the merged dataset to be used as the training dataset. The remaining 19% of the rows will be used as the test dataset.
train_idx <- sample(nrow(merge_data), nrow(merge_data) * 0.81)
#creates the training dataset 
train_data <- merge_data[train_idx, ] 
#creates the test dataset 
test_data <- merge_data[-train_idx, ]    

# Train SVM model using the radial kernel
svm_model <- svm(WHO_severity ~ ., data = train_data,kernal ="radial")
svm_predictions <- predict(svm_model, newdata = test_data)    #used trained SVM model to make predictions
confusion_matrix_SVM <- confusionMatrix(svm_predictions, test_data$WHO_severity)
confusion_matrix_SVM$table

#calculates the overall accuracy of the SVM model  [get 0.5 accuracy ]
accuracy_SVM <- sum(diag(confusion_matrix_SVM$table)) / sum(confusion_matrix_SVM$table)
print(paste0("Accuracy: ", round(accuracy_SVM, 3)))
####create important features
train_model <- train(WHO_severity ~ ., data = train_data, method = "svmRadial")
svm_varimp <- varImp(train_model, scale = FALSE)
svm_varimp$importance

# Sort the features by their importance score and Select the top 50 most important features
sorted_importance_SVM <- row.names(svm_varimp$importance)[order(svm_varimp$importance[,1], decreasing = TRUE)[1:50]]
sorted_importance_SVM

#download file
write.table(sorted_importance_SVM,"top50_features_SVM_1.txt",sep="\t")

###plot Confusion Matrix of SVM model for all 3 grp using ggplot###
ggplot(data = data.frame(confusion_matrix_SVM$table, class = rownames(confusion_matrix_SVM$table)),
       aes(x = Prediction, y = Reference, fill = Freq)) + 
  geom_tile() + 
  geom_text(aes(label = Freq), size = 10, colour = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  xlab("Predicted Class") + 
  ylab("True Class") + 
  ggtitle("Confusion Matrix SVM 1") + 
  theme(plot.title = element_text(hjust = 0.5))

###################################################################################################################################################
#classification performance only between mild and critical_severe groups.
#1. Using Random Forest

merged_data_sub <- merge_data[merge_data$WHO_severity %in% c("1", "3"),]
merged_data_sub
merged_data_sub$WHO_severity <- factor(merged_data_sub$WHO_severity)

###builds a random forest model using the training dataset
set.seed(123)
trainIndex <- createDataPartition(merged_data_sub$WHO_severity, p=0.7, list=FALSE)
trainData <- merged_data_sub[trainIndex,]
testData <- merged_data_sub[-trainIndex,]
testData$Sample_ID<-NULL
trainData$Sample_ID<-NULL
rf_model <- randomForest(WHO_severity ~ ., data = trainData, importance = TRUE, ntree = 1000,type="classification")
predicted_severity <- predict(rf_model, newdata=testData)
conf_mat <- confusionMatrix(predicted_severity, testData$WHO_severity)
print(conf_mat$table)

###plot Confusion Matrix of Random Forest model for mild and critical_severe groups using ggplot###
ggplot(data = data.frame(conf_mat$table, class = rownames(conf_mat$table)),
       aes(x = Prediction, y = Reference, fill = Freq)) + 
  geom_tile() + 
  geom_text(aes(label = Freq), size = 10, colour = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  xlab("Predicted Class") + 
  ylab("True Class") + 
  ggtitle("Confusion Matrix RF 2") + 
  theme(plot.title = element_text(hjust = 0.5))

#calculates the overall accuracy of the Random Forest model
accuracy_2 <- sum(diag(conf_mat$table)) / sum(conf_mat$table)
print(paste0("Accuracy: ", round(accuracy_2, 3)))
importance_2 <- importance(rf_model,type=2)
importance_2

# Sort the features by their importance score and Select the top 50 most important features
sorted_importance_2 <- row.names(importance_2)[order(importance_2[,1], decreasing = TRUE)[1:50]]
sorted_importance_2
top50_features_RF_2 <- names(merged_data_sub[,sorted_importance_2])
print(top50_features_RF_2)

write.table(top50_features_RF_2,"top50_features_RF_2.txt",sep="\t")  #download

################################################################################################################
#classification performance only between mild and critical_severe groups.
#2.Using Support Vector Machines(SVM)

merged_data_sub <- merge_data[merge_data$WHO_severity %in% c("1", "3"),]
merged_data_sub
merged_data_sub$WHO_severity <- factor(merged_data_sub$WHO_severity)

#creates a character vector
non_numeric_vars <- c("Age", "Sex")

###for categorical data analysis###
#converts each variable in non_numeric_vars to a factor variable for perform categorical data analysis on these variables.
for (var in non_numeric_vars) {
  merged_data_sub[[var]] <- as.factor(merged_data_sub[[var]])
}
for (var in names(merged_data_sub)) {
  if (is.factor(merged_data_sub[[var]]) && length(levels(merged_data_sub[[var]])) < 2) {
    merged_data_sub[[var]] <- NULL
  }
}
for (var in names(merged_data_sub)) {
  if (is.character(merged_data_sub[[var]]) || is.factor(merged_data_sub[[var]])) {
    merged_data_sub[[var]] <- factor(merged_data_sub[[var]])
  }
}

###builds a SVM model using the training dataset
set.seed(123)
train_idx_svm <- sample(nrow(merged_data_sub), nrow(merged_data_sub) * 0.81)
train_data_svm <- merged_data_sub[train_idx, ]
train_data_svm
test_data_svm <- merged_data_sub[-train_idx, ]
train_data_svm <- train_data_svm[complete.cases(train_data_svm), ]
test_data_svm<-test_data_svm[complete.cases(test_data_svm), ]
any(is.na(test_data_svm))

#trains an SVM model using the radial kernel
svm_model_2 <- svm(WHO_severity ~ ., data = train_data_svm,kernal ="radial")
merged_data_sub$WHO_severity <- factor(merged_data_sub$WHO_severity)
svm_predictions_2 <- predict(svm_model_2, newdata = test_data_svm)
confusion_matrix_SVM_2 <- confusionMatrix(svm_predictions_2, test_data_svm$WHO_severity)
confusion_matrix_SVM_2$table

#calculates the overall accuracy of the SVM model
accuracy_SVM_2 <- sum(diag(confusion_matrix_SVM_2$table)) / sum(confusion_matrix_SVM_2$table)
print(paste0("Accuracy: ", round(accuracy_SVM_2, 3)))

train_model_2 <- train(WHO_severity ~ ., data = train_data_svm, method = "svmRadial")
svm_varimp_2 <- varImp(train_model_2, scale = FALSE)
svm_varimp$importance

# Sort the features by their importance score and Select the top 50 most important features
sorted_importance_SVM <- row.names(svm_varimp$importance)[order(svm_varimp$importance[,1], decreasing = TRUE)[1:50]]
sorted_importance_SVM

#Download file
write.table(sorted_importance_SVM,"top50_features_SVM_2.txt",sep="\t")

###plot Confusion Matrix of SVM model for mild and critical_severe groups using ggplot###
ggplot(data = data.frame(confusion_matrix_SVM_2$table, class = rownames(confusion_matrix_SVM_2$table)),
       aes(x = Prediction, y = Reference, fill = Freq)) + 
  geom_tile() + 
  geom_text(aes(label = Freq), size = 10, colour = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  xlab("Predicted Class") + 
  ylab("True Class") + 
  ggtitle("Confusion Matrix SVM 2") + 
  theme(plot.title = element_text(hjust = 0.5))
### Inference : By Comparing the classification model , Radom Forest shows high accuracy compared to SVM when investigating for two groups of covid_severity and also for three groups of covid_severity

######################################################################################################################################

