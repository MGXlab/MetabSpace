8# CCV approach with intensity weights
#Pilleriin Peets (pilleriin.peets@ut.ee,  pilleriin.peets@gmail.com)
# Code for calculating chemical characteristics vectors (CCV) from SIRIUS predicted fingerpints and canopus vectors


#-------libraries------------

library(tidyverse)
library(data.table)
library(ggplot2)
library(viridis)
library(janitor)
library(umap)
library(plotly)
library(rjson)
library(rcdk)
library(scales)
library(sunburstR)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(plotly)
library(vegan)
library(svglite)
library(readxl)
library(tidytable)



#---------reading in data--------------

#Example data for testing in "https://github.com/MGXlab/MetabSpace" -> "Data" -> "Intensity_weighted_CCV"

folder <- "C:/Users/b15157/OneDrive - Tartu Ülikool/UniJena/AcidSinkhole/GitHub_test"
setwd(folder)

#SIRIUS predicted vectors with LC-MS feature ID that matched intensity table
#for CCvectors and MFP vectors from SIRIUS data use functions MFtable_SIRIUS5 and CCtable_SIRIUS5 below
#For SIRIUS6 vectors must be extracted from API into a table. Un and Canop id numbers are from absoluteIndex column in metadata files from SIRIUS.

CCvectors <- fread("20240628_canopus_canoprank.csv")
MFPvectors <- fread("20240628_fingerprint_canoprank.csv")

#MS2 peaks information from LC-MS analysis (e.g. intensity table from MZmine)

#this specific example is for MZmine .csv output.
ms2peaks <- fread("20241106_aligned_main_ms2.csv") %>%
  select(!contains("O2005")) %>%
  dplyr::filter(rt < 16.9) %>%
  dplyr::filter(mz < 1000) %>%
  mutate(peakcheck = height/area) %>%
  select(id,rt,mz,height,area,peakcheck,contains("area")) %>%
  replace(is.na(.), 0)


# selecting only id and intensities and creating long table
intensitytable <- ms2peaks %>%
  mutate(featureId = id) %>%
  select(featureId, contains(":area")) %>%
  gather(key = "sampleID", value = "intensity", -featureId)

#adjusts the sample name, as MZmine column has file extension and are info 
intensitytable$sampleID <- str_extract(intensitytable$sampleID, "(?<=datafile:)[^.]+")


#------------Function continue from other files-----------


samplelist <- unique(intensitytable$sampleID)
sample <- samplelist[1]
sample


IntWeight_CCV_MFP <- tibble()
for (sample in samplelist) {
  filtered_list <- intensitytable %>%
    dplyr::filter(sampleID == sample) %>%  #chooses sample n in list
    na.omit() %>% 
    dplyr::filter(featureId %in% MFPvectors$featureId) %>%
    group_by(sampleID) %>% 
#calculates relative intensity for each peak
    mutate(intsum = sum(intensity)) %>%
    ungroup() %>%
    mutate(intnorm = intensity/intsum) %>%
    select(featureId, sampleID,intensity, intnorm) %>%
# adds molecular fingerprint information for each LC-MS feature from SIRIUS tables
    left_join(MFPvectors, by = "featureId") %>%
    na.omit() %>% 
# adds the weight to each bit in a SIRIUS vector
    mutate_at(vars(starts_with("Un")), ~. * intnorm) %>%
    select(featureId, sampleID, starts_with("Un"))
# sums the intensities. while is actual averaging in not done, the summing is giving the same results as multiplying the peaks and averaging
  averages <- filtered_list %>%
    summarise(across(contains("Un"), sum), .groups = "drop") %>%
    mutate(sampleID = sample) %>%
    select(sampleID, everything())
  
  IntWeight_CCV_MFP  <- bind_rows(IntWeight_CCV_MFP , averages)
  
}

#IntWeight_CCV_MFP  table has now one CCV vector per sample


#for Canop classes the same

IntWeight_CCV_CC <- tibble()
for (sample in samplelist) {
  filtered_list <- intensitytable %>%
    dplyr::filter(sampleID == sample) %>%  #chooses sample n in list
    na.omit() %>% 
    dplyr::filter(featureId %in% CCvectors$featureId) %>%
    group_by(sampleID) %>% 
    #calculates relative intensity for each peak
    mutate(intsum = sum(intensity)) %>%
    ungroup() %>%
    mutate(intnorm = intensity/intsum) %>%
    select(featureId, sampleID,intensity, intnorm) %>%
    # adds molecular fingerprint information for each LC-MS feature from SIRIUS tables
    left_join(CCvectors, by = "featureId") %>%
    na.omit() %>% 
    # adds the weight to each bit in a SIRIUS vector
    mutate_at(vars(starts_with("Canop")), ~. * intnorm) %>%
    select(featureId, sampleID, starts_with("Canop"))
  # sums the intensities. while is actual averaging in not done, the summing is giving the same results as multiplying the peaks and averaging
  averages <- filtered_list %>%
    summarise(across(contains("Canop"), sum), .groups = "drop") %>%
    mutate(sampleID = sample) %>%
    select(sampleID, everything())
  
  IntWeight_CCV_CC  <- bind_rows(IntWeight_CCV_CC , averages)
  
}



#------PCA plots for testing-------------

all_aver_classif <- IntWeight_CCV_CC %>%
  na.omit() %>%
  column_to_rownames("sampleID")

constant_cols <- apply(all_aver_classif, 2, function(col) length(unique(col)) == 1)
all_aver_classif <- all_aver_classif[, !constant_cols]

zero_cols <- apply(all_aver_classif, 2, function(col) all(col == 0))
all_aver_classif <- all_aver_classif[, !zero_cols]

# setwd(paste(admin, "LCMS_measured/Sinkhole/RESULTS/Tables", sep = "/"))
# formeine <- as.data.frame(t(all_aver_classif)) %>%
#   rownames_to_column("Canop_av")
# write_delim(formeine, "20241104_rfe2CanopIntAver_allsamples.csv", delim = ";")

#pca and umap
set.seed(123)
data_pca <- prcomp(all_aver_classif, rank = 3)

pca_results <- as.data.frame(data_pca$x) %>%
  rownames_to_column("sampleID")
PC1var <- round(summary(data_pca)$importance[2,1]*100,1)
PC2var <- round(summary(data_pca)$importance[2,2]*100,1)

graphdata <- pca_results 
#add here metadata related to samples

plot_pca <- ggplot(data = graphdata) +
  geom_point(mapping = aes(x = PC1,
                           #shape = cluster,
                           y = PC2),
             size = 6, alpha = 0.7, color = "black") + # NA points
  theme_minimal() +
  coord_fixed(ratio = 1) +
  scale_color_viridis(option = "H", direction = -1) +
  ggtitle("Name") +
  #facet_wrap(~ SampleName, nrow = 7) +
  xlab(paste("PC1, ", PC1var, "%", sep = "")) +
  ylab(paste("PC2, ", PC2var, "%", sep = ""))
plot_pca






# Functions to extract SIRIUS fingerprints and canopus vectrors !! works only for SIRIUS5

MFtable_SIRIUS5 <- function(folderwithSIRIUSfiles){
  
  setwd(folderwithSIRIUSfiles)
  subfolder_fp_zip <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = "fingerprints")
  
  fp_names_pos <- paste("Un", read_delim(paste(folderwithSIRIUSfiles,"/csi_fingerid.tsv", sep = ""), delim = "\t")$absoluteIndex, sep = "")
  fp_names_neg <- paste("Un", read_delim(paste(folderwithSIRIUSfiles,"/csi_fingerid_neg.tsv", sep = ""), delim = "\t")$absoluteIndex, sep = "")
  
  ranks <- fread("formula_identifications_all.tsv") %>%
    select(id, precursorFormula, adduct, formulaRank, molecularFormula, featureId) %>%
    mutate(adduct = gsub(" ", "", adduct))
  
  
  all_fp <- tibble()
  for(zipF in subfolder_fp_zip){
    
    zipFile <- paste(folderwithSIRIUSfiles, zipF, sep = "/")
    zip_contents <- unzip(zipFile, list = TRUE)
    fp_files <- zip_contents$Name[grep("\\.fpt$", zip_contents$Name)]
    
    fp_data <- lapply(fp_files, function(fp_file) {
      file_content <- readLines(unz(zipFile, fp_file))
      return(file_content)
    })
    
    fp_df <- do.call(rbind, lapply(fp_data, function(x) as.data.frame(t(x)))) 
    
    fp_df <- fp_df %>%
      mutate_all(as.numeric)
    
    fp_treshold = 0.5
    fp_df[fp_df < fp_treshold & fp_df > 0] <- 0
    fp_df[fp_df >= fp_treshold & fp_df <= 1] <- 1
    
    fp_df <- cbind(FPformula = fp_files, fp_df) %>%
      mutate(file = zipF) %>%
      select(file, everything())
    
    all_fp <- all_fp %>%
      bind_rows(fp_df)
    
  }
  
  if(grepl("]-.fpt", all_fp$FPformula[1])) {
    colnames(all_fp) <- c("file", "FPformula", fp_names_neg)
  }
  
  if(grepl("]+.fpt", all_fp$FPformula[1])) {
    colnames(all_fp) <- c("file", "FPformula", fp_names_pos)
  }
  
  final_fp <- all_fp %>%
    separate(file, into = c("id"), sep = "/") %>%
    separate(FPformula, into = c("precursorFormula", "adduct"), sep = "_") %>%
    separate(adduct, into = c("adduct"), sep = ".fp")
  
  final_fp_ranks <- ranks %>%
    inner_join(final_fp, by = c("id", "precursorFormula", "adduct"))
  
  return(final_fp_ranks)
}


CCtable_SIRIUS5 <- function(folderwithSIRIUSfiles){
  
  setwd(folderwithSIRIUSfiles)
  subfolder_fp_zip <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = "canopus")
  subfolder_fp_zip <- subfolder_fp_zip[!grepl("_npc", subfolder_fp_zip)]
  subfolder_fp_zip <- subfolder_fp_zip[!grepl(".tsv", subfolder_fp_zip)]
  
  fp_names_pos <- paste("Canop", read_delim(paste(folderwithSIRIUSfiles,"/canopus.tsv", sep = ""), delim = "\t")$absoluteIndex, sep = "")
  fp_names_pos <- paste("Canop", read_delim(paste(folderwithSIRIUSfiles,"/canopus_neg.tsv", sep = ""), delim = "\t")$absoluteIndex, sep = "")
  
  ranks <- fread("formula_identifications_all.tsv") %>%
    select(id, precursorFormula, adduct, formulaRank, molecularFormula, featureId) %>%
    mutate(adduct = gsub(" ", "", adduct))
  
  
  all_fp <- tibble()
  for(zipF in subfolder_fp_zip){
    
    zipFile <- paste(folderwithSIRIUSfiles, zipF, sep = "/")
    zip_contents <- unzip(zipFile, list = TRUE)
    fp_files <- zip_contents$Name[grep("\\.fpt$", zip_contents$Name)]
    
    fp_data <- lapply(fp_files, function(fp_file) {
      file_content <- readLines(unz(zipFile, fp_file))
      return(file_content)
    })
    
    fp_df <- do.call(rbind, lapply(fp_data, function(x) as.data.frame(t(x)))) 
    
    fp_df <- fp_df %>%
      mutate_all(as.numeric)
    
    fp_treshold = 0.5
    fp_df[fp_df < fp_treshold & fp_df > 0] <- 0
    fp_df[fp_df >= fp_treshold & fp_df <= 1] <- 1
    
    fp_df <- cbind(FPformula = fp_files, fp_df) %>%
      mutate(file = zipF) %>%
      select(file, everything())
    
    all_fp <- all_fp %>%
      bind_rows(fp_df)
    
  }
  
  if(grepl("]-.fpt", all_fp$FPformula[1])) {
    colnames(all_fp) <- c("file", "FPformula", fp_names_neg)
  }
  
  if(grepl("]+.fpt", all_fp$FPformula[1])) {
    colnames(all_fp) <- c("file", "FPformula", fp_names_pos)
  }
  
  final_fp <- all_fp %>%
    separate(file, into = c("id"), sep = "/") %>%
    separate(FPformula, into = c("precursorFormula", "adduct"), sep = "_") %>%
    separate(adduct, into = c("adduct"), sep = ".fp")
  
  final_fp_ranks <- ranks %>%
    inner_join(final_fp, by = c("id", "precursorFormula", "adduct"))
  
  return(final_fp_ranks)
}