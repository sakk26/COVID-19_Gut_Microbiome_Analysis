rm(list=ls())
setwd("C:\\Users\\sakshi Mahajan\\Desktop\\hmds_assi1")
getwd()
###load required packages###
library(vegan)
library(tidyverse)
library(dplyr)
library(ggplot2)
####upload given data files### 
clr_species_abundance<-read.table("Assignment1_ClrTrans_Species.txt",sep="\t",header=TRUE)
metadata_species<-read.table("Assignment1_Metadata.txt",sep="\t",header=TRUE)
###convert WHO_severity column into numeric values
metadata_species$WHO_severity <- as.numeric(factor(metadata_species$WHO_severity, levels=c("mild", "moderate", "critical_severe")))
##creating a dataframe of sample ID and WHO_severity
divesity_severity <- metadata_species %>% select(Sample_ID,WHO_severity )
##merge clr abundance with Covid_severity
merge_data<-merge(clr_species_abundance,divesity_severity,by="Sample_ID")
##########################5% detected ########################################################################
#covert clr abundance data into binary data
binary_data <- ifelse(clr_species_abundance > 0, 1, 0)    #value greater that 0 is assign as 1 else put 0 
col_sums<-colSums(binary_data)                      #overall column sum
sort(col_sums)
col_sums
#filter out only the detected species in at-least 5% of the samples.
at_least_5percent <- which(col_sums >= 0.05 * nrow(binary_data))   
at_least_5percent
filtered_species<-clr_species_abundance[,at_least_5percent]

##############################################################################################################
#linear regression analysis for each species in a data frame 
regression_results <- list()
for (species_name in colnames(filtered_species)) {
  # Skip the Sample_ID column
  if (species_name == "Sample_ID") {
    next
  }
  lm_model <- lm(filtered_species[,species_name] ~ WHO_severity, data = metadata_species)
  # Store the regression result in the list
  regression_results[[species_name]] <- summary(lm_model)
}
### #extract the p-values from the regression results for each species.
p_values <- sapply(regression_results, function(x) coef(x)[2,4]) 
# adjust the p-values for multiple comparisons using the Benjamini-Hochberg method. 
p_values_adjusted <- p.adjust(p_values, method = "BH")  
p_values_adjusted
#sort the adjusted p-values in ascending order.
sort(p_values_adjusted)  
#Extract only the  names of the species with adjusted p-values less than or equal to 0.1.
significant_species <- names(p_values_adjusted[p_values_adjusted <= 0.1])    
print(significant_species)
write.table(significant_species,"significant_species.txt",sep="\t")

########################################################################################################
# Doing for confounders [i.e.WHO_severity + Age + BMI + Sex + comorbidities_total]
#Multiple linear regression model
regression_results_confounder <- list()       
for (species_name in colnames(filtered_species)) {
  # Skip the Sample_ID column
  if (species_name == "Sample_ID") {
    next
  }
  lm_model <- lm(filtered_species[,species_name] ~ WHO_severity + Age + BMI + Sex + comorbidities_total, data = metadata_species)
  # Store the regression result in the list
  regression_results_confounder[[species_name]] <- summary(lm_model)
}
#extract the p-values from the regression results for each species.
p_values_confounder <- sapply(regression_results_confounder, function(x) coef(x)[2,4])  
#adjust the p-values for multiple comparisons using the Benjamini-Hochberg method.
p_values_adjusted_confounder <- p.adjust(p_values, method = "BH")  
#Extract the names of the species with adjusted p-values less than or equal to 0.1.
significant_species_confounder <- names(p_values_adjusted[p_values_adjusted <= 0.1]) 
print(significant_species_confounder)
write.table(significant_species_confounder,"significant_species_confounder.txt",sep="\t")
####As a final result without and with confounder gives same results that means there is no influence of confounders on this  
