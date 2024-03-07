rm(list=ls())
setwd("C:\\Users\\sakshi Mahajan\\Desktop\\hmds_assi1")
getwd()
###load required packages###
library(vegan)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(geosphere)
####upload given data files### 
clr_species_abundance<-read.table("Assignment1_ClrTrans_Species.txt",sep="\t",header=TRUE)
metadata_species<-read.table("Assignment1_Metadata.txt",sep="\t",header=TRUE)
###convert WHO_severity column into numeric values
metadata_species$WHO_severity <- as.numeric(factor(metadata_species$WHO_severity, levels=c("mild", "moderate", "critical_severe")))
rownames(clr_species_abundance)<-clr_species_abundance$Sample_ID
clr_species_abundance<-clr_species_abundance[,-1]
###create empty vector with following coloumn names###############################################################
metadata_vars <- c("Sex", "Age", "BMI", "WHO_severity", "comorbidities_total", 
                   "HTN", "Diabetes", "Respiratory_disease", "Heart_disease", 
                   "Renal_Disease", "Liver_Disease", "Obesity", "Malignancy", 
                   "Immunosuppressive_Disease", "Neurological_disease", "Metabolic_Disease", "Cardiovascular_Disease")
######################################################################################################################
#calculating the Euclidean distance between clr_species_abundance using the dist function.
Euclidean_dist <- dist(clr_species_abundance, method = "euclidean") 

#calculating the bray distance between clr_species_abundance using the vegdist function.
bray_curtis_dist <- vegdist(clr_species_abundance, method = "bray")
#creating an empty data frame
permanova_results=data.frame() 
permanova_euclidean=data.frame()   
#creating an empty data frame
p_values=data.frame()
p_values_euclidean=data.frame()
###############################################################################################################
#PERMANOVA using the Bray-Curtis distance
#loop that will iterate through each column in "metadata_vars"
for (col in metadata_vars) {                         
  formula <- as.formula(paste("bray_curtis_dist ~ metadata_species$", col))
  results <- adonis2(formula, data = metadata_species)
  #extracts the p-value 
  p_value<-results$"Pr(>F)"[1]   
  print(p_value)
  p_values<-cbind(p_value,col)
  print(p_values)
  permanova_results <- rbind(permanova_results, p_values)
}
permanova_results <- permanova_results[order(permanova_results[,1]),]

#selects rows from permanova_results where the p-value is less than or equal to 0.05
significant_results_bray <- permanova_results[permanova_results[, 1] <= 0.05, ]    
print(significant_results_bray)

write.table(significant_results_bray,"significant_metadata_braycurtis.txt",sep="\t")
########################################################################################################################
#PERMANOVA analysis using Euclidean distance
#loop that will iterate through each column in "metadata_vars"
for (col in metadata_vars) {                                
  formula <- as.formula(paste("Euclidean_dist ~ metadata_species$", col))
  results <- adonis2(formula, data = metadata_species)
  p_values_euclidean<-results$"Pr(>F)"[1]     #extract p-value from the PERMANOVA analysis
  print(p_values_euclidean)
  p_values_euclidean<-cbind(p_values_euclidean,col)
  print(p_values_euclidean)
  permanova_euclidean <- rbind(permanova_euclidean, p_values_euclidean)
}
permanova_euclidean <- permanova_euclidean[order(permanova_euclidean[,1]),]

#selects rows from permanova_results where the p-value is less than or equal to 0.05
significant_results_euclidean <- permanova_results[permanova_results[, 1] <= 0.05, ]  
print(significant_results_euclidean)
write.table(significant_results_euclidean,"significant_metadata_euclidean.txt",sep="\t")
##we get same  result using both the distances
#####Inference  : four metadata (HTN,WHO_severity,Malignancy and Respiratory_disease) shows significant results across clr abundance data

























Assi 2
rm(list=ls())
setwd("C:\\Users\\sakshi Mahajan\\Desktop\\hmds\\HMDS_ASSI 2")
getwd()
#df_ForMetaAnalysis <- read.csv("df_ForMetaAnalysis.csv",row.names=1,header=TRUE)
#SelectedStudies <- read.csv("selected_studies.csv",row.names=1,header=TRUE)

#install.packages("metafor")
library(metafor)

# Subset the data frame to include only the selected studies and age >= 60
df <- subset(df_ForMetaAnalysis, study_name %in% SelectedStudies & age >= 60)

# Create an empty data frame to store the results
results <- data.frame(Species = character(length(SelectSpecies)),
                      `RE-Model Summary Estimate` = numeric(length(SelectSpecies)),
                      P_value = numeric(length(SelectSpecies)),
                      `Corrected P-value` = numeric(length(SelectSpecies)),
                      stringsAsFactors = FALSE)

results <- data.frame(matrix(nrow = nrow(df), ncol = 4))
results <- data.frame(matrix(nrow = length(SelectSpecies), ncol = 4))
colnames(results) <- c("Species", "RE-Model Summary Estimate", "P-value", "Corrected P-value")

# Loop through each species and perform random effects model meta-analysis
for (i in 1:length(SelectSpecies)) {
  species <- SelectSpecies[i]
  y <- df[, species]
  n <- df$sample_size
  study <- df$study_name
  res <- rma(y, n, mods = ~ age, data = df, subset = age >= 60)
  results[i, 1] <- species
  results[i, 2] <- res$b
  results[i, 3] <- summary(res)$pval
}

# Add column for corrected p-value using Benjamini-Hochberg correction
results[, 4] <- p.adjust(results[, 3], method = "BH")

# Print the results
print(results)
rm(list=ls())
setwd("C:\\Users\\sakshi Mahajan\\Desktop\\hmds\\HMDS_ASSI 2")
getwd()
load("Assignment2.RData")
df_ForMetaAnalysis <- read.csv("df_ForMetaAnalysis.csv",row.names=1,header=TRUE)
SelectedStudies <- read.csv("selected_studies.csv",row.names=1,header=TRUE)

install.packages("meta")
library(meta)
library(dplyr)
df_ForMetaAnalysis_filtered <- df_ForMetaAnalysis %>%
  filter(age >= 60, study_name %in% SelectedStudies)
results <- data.frame()
#(xtot-x) mean mass conservation i.e. one active and one inactive
#-k11*x*(z^n/(EC50^n + z^n))it is minus of h+
#####
# Subset the data frame to only include samples with age >= 60 and in the selected studies
df_ForMetaAnalysis_filtered <- subset(df_ForMetaAnalysis, age >= 60 & study_name %in% SelectedStudies)

df_ForMetaAnalysis_filtered <- subset(df_ForMetaAnalysis, age >= 60 & study_name %in% SelectedStudies & species %in% SelectSpecies)

df_ForMetaAnalysis_filtered <- df_ForMetaAnalysis %>%
  filter(age >= 60, study_name %in% SelectedStudies)

table(df_ForMetaAnalysis_filtered$species %in% SelectSpecies)
###
for (s in SelectSpecies) {
  df_species <- df_ForMetaAnalysis_filtered %>%
    filter(species == s) %>%
    select(study_name, abundance)
  
  meta <- metagen(abundance, se = 0, studlab = study_name, data = df_species, method.tau = "DL", method.ci = "DL")
  
  estimate <- meta$TE.random
  se <- meta$seTE.random
  weight <- meta$weights
  
  results <- rbind(results, data.frame(species = s, estimate = estimate, se = se, weight = weight))
}
results$pvalue <- 2 * pnorm(-abs(results$estimate / results$se))
results$FDR <- p.adjust(results$pvalue, method = "fdr")
##############################################################################################################
#install.packages("ggrepel")

library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(MASS)
library(gplots)
library(RColorBrewer)
library(sfsmisc)
library(ggplot2)
library(ggrepel)
library(sfsmisc)
library(compositions)





compute_meta_corr_group <- function(data,feature_list,metadata_var,grouping_var,grouping_list)
{
  return_out <- as.data.frame(matrix(NA,length(feature_list),10))
  rownames(return_out) <- feature_list
  colnames(return_out) <- c("beta","pval","ci.ub","ci.lb","tau2","QE","QEp","qval","dir","consistency")
  return_out[,1] <- 0
  return_out[,2] <- 1
  return_out[,3] <- 0
  return_out[,4] <- 0
  return_out[,5] <- 0
  return_out[,6] <- 0
  return_out[,7] <- 1
  return_out[,10] <- 0
  
  for(i in 1:length(feature_list))
  {
    species_name <- feature_list[i]
    #print(species_name)
    tryCatch(               
      expr = {                     
        temp_res <- compute_meta_corr(data,species_name,metadata_var,grouping_var,grouping_list)
        print(species_name)
        return_out[i,"beta"] <- temp_res$model$beta
        return_out[i,"pval"] <- temp_res$model$pval
        return_out[i,"ci.ub"] <- temp_res$model$ci.ub
        return_out[i,"ci.lb"] <- temp_res$model$ci.lb
        return_out[i,"tau2"] <- temp_res$model$tau2
        return_out[i,"QE"] <- temp_res$model$QE
        return_out[i,"QEp"] <- temp_res$model$QEp
        return_out[i,"consistency"] <- length(which(sign(temp_res$df_studies[temp_res$df_studies$ri!=0,"ri"])==sign(as.numeric(temp_res$model$beta))))/length(temp_res$df_studies[temp_res$df_studies$ri!=0,"ri"])
        
      },
      error = function(e){    
        print(e)
        print("Error observed. Moving to next")
      },
      finally = {            
        print("finally Executed")
      }
    )
  }
  return_out$qval <- p.adjust(return_out$pval,method="fdr")
  return_out$dir <- ifelse(return_out$qval <= 0.1,3*sign(return_out$beta),ifelse(return_out$pval <= 0.08,2*sign(return_out$beta),sign(return_out$beta)))
  return_list <- list("model" = temp_res$model,"df_studies"=temp_res$df_studies)
  #return(return_list)
  return(return_out)
}

