###Assignment 1###
#Quesion 1####################################################################################################
rm(list=ls())
setwd("C:\\Users\\sakshi Mahajan\\Desktop\\hmds_assi1")
getwd()
###install required packages###
#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("picante")
#install.packages("ggplot2")
library(vegan)
library(tidyverse)
library(dplyr)
library(picante)
library(ggplot2)

#upload data files 
Assignment_RawCount_Species<-read.table("Assignment_RawCount_Species.txt",sep="\t",header=TRUE)
clr_species_abundance<-read.table("Assignment1_ClrTrans_Species.txt",sep="\t",header=TRUE)
metadata_species<-read.table("Assignment1_Metadata.txt",sep="\t",header=TRUE)

#assign rownames as sample_ID and then removed sample id seperate coloumb
rownames(Assignment_RawCount_Species)<-Assignment_RawCount_Species$Sample_ID
Assignment_RawCount_S<-Assignment_RawCount_Species[,-1]
rownames(metadata_species)<-metadata_species$Sample_ID
#creating a dataframe of sample ID and WHO_severity
divesity_severity <- metadata_species %>% select(Sample_ID,WHO_severity )
###Shannon Index##############################################################################################
#Que 1[a]

#normalizes the matrix
RawCount_Species_N<-Assignment_RawCount_S/rowSums(Assignment_RawCount_S)  
RawCount_Species_N

# calculates the Shannon diversity index for each sample in the normalized matrix 
Shannon_diversity_species<-diversity(RawCount_Species_N,index="shannon")   
Shannon_diversity_species
#calculates Pielou's evenness index for each sample in the normalized matrix
Pielou_evennes_df <- Shannon_diversity_species/log(specnumber(RawCount_Species_N))  
Pielou_evennes_df


###creating a dataframe of shannon,Pileou and severity as alpha diversity###
alpha_diversity_df<-data.frame(Shannon_diversity_species,Pielou_evennes_df)  #combine both in dataframe
alpha_diversity_df[,3]<-rownames(alpha_diversity_df)  
names(alpha_diversity_df)[3]<-"Sample_ID"
names(alpha_diversity_df)[1]<-"Shannon_index"

###merge alpha_diverity and divesity_severity on basis of sample ID###
common_samples_meta = intersect(rownames(alpha_diversity_df),metadata_species$Sample_ID)
alpha_diverity<-alpha_diversity_df[common_samples_meta,]
merge_df<-merge(alpha_diverity,divesity_severity,by = "Sample_ID")

###############################################################################################################
# One-way ANOVA for Shannon index
shannon_aov <- aov(merge_df$Shannon_index~merge_df$WHO_severity,data=merge_df)
summary(shannon_aov)

###Inference: Shannon_aov P-values is greater than 0.05 which means there is no significant difference between the sample and covid_severity
# One-way ANOVA for Peilou index
Peilou_evenness_aov <- aov(merge_df$Pielou_evennes_df~merge_df$WHO_severity,data=merge_df)
summary(Peilou_evenness_aov)
###Inference: Pielou_evenness Anova, P-values is greater than 0.05 which means there is no significant difference between the sample and covid_severity

###Boxplot result###
#Create box plots for each index for Shannon_index
ggplot(merge_df, aes(x = WHO_severity, y = Shannon_index)) + 
  geom_boxplot() +
  labs(title = "Shannon Index by Severity Group", x = "Severity", y = "Shannon Index")

#Create box plots for each index for Pielou_index
ggplot(merge_df, aes(x = merge_df$WHO_severity, y = merge_df$Pielou_evennes_df)) + 
  geom_boxplot() +
  labs(title = "Pielou Index by Severity Group", x = "Severity", y = "Pielou Index")
##############################################################################################################
#Que 1[b] PCOA plot

###beta_diversity###
#Calculate Bray-Curtis index for the RawCount_Species_N data matrix using the vegdist() function
beta_diversity<-(vegdist(RawCount_Species_N,method="bray"))  
beta_pcoa <- pcoa(beta_diversity)

#extracts the vectors from the beta_pcoa 
B_vectors<-beta_pcoa$vectors   
common_samples_meta = intersect(rownames(RawCount_Species_N),metadata_species$Sample_ID)  
metadata_species_c<-B_vectors[common_samples_meta,]
df_pcoa<-data.frame(metadata_species_c,metadata_species)  #combines the PCoA coordinates 
###Pcoa plot###
pcoa_plot <- ggplot(df_pcoa, aes(x=Axis.1, y=Axis.2, color=metadata_species$WHO_severity)) +
  geom_point(size=3) +
  labs(x="PCoA 1", y="PCoA 2", color="WHO_severity") +
  ggtitle("Bray-Curtis PCoA plot")
print(pcoa_plot)
####################################################################################################################
###The test we have to use for this is PERMOANOVA### 
###PERMOANOVA###
metadata<-read.table("Assignment1_Metadata.txt",sep="\t",header=TRUE)
#convert WHO_severity coloumn in numeric format
metadata$WHO_severity <- as.numeric(factor(metadata$WHO_severity, levels=c("mild", "moderate", "critical_severe")))   
as.data.frame(RawCount_Species_N)
#Creating distance between the abundance matrix 
RAW<-RawCount_Species_N[common_samples_meta,]
#calculating the bray distance between Raw data using the vegdist function.
beta_diversity_raw<-(vegdist(RAW,method="bray"))

permanova <- adonis2(beta_diversity_raw ~ metadata$WHO_severity, data=metadata)
summary(permanova)
###Inference: permanova P-values are greater than 0.05 which means there is no significant difference between the sample as a whole and covid_severity


