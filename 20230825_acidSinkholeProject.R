###Analysing Meine and Milous data

library(tidyverse)
library(data.table)
library(ggplot2)
library(viridis)
library(janitor)
library(umap)
library(plotly)
library(Rdisop)
library(rjson)
library(patRoon)
library(rcdk)
library(scales)


#------------Functions-------------

admin <- "C:/Users/pille/OneDrive - Friedrich-Schiller-Universität Jena/JenaUniversity"
setwd(paste(admin, "RcodesJena", sep = "/"))
source("Graph/my_theme_1.R")
source("Functions_SIRIUS_DataAnalysis.R")


#-----------SIRIUS fp and classes-----------
setwd(paste(admin, "Datasets_fp/GNPS_Kai", sep = "/"))
sir_fp_pos_Kai <- fread("csi_fingerid.tsv") %>%
  mutate(Un = paste("Un", absoluteIndex, sep = ""))
sir_canop_pos_Kai <- fread("canopus.tsv") %>%
  mutate(Canop = paste("Canop", absoluteIndex, sep = ""))

#------------metadata-------------
setwd(paste(admin, "LCMS_measured/Milou_samples", sep = "/"))

sample_info_extra <- read_delim("Sinkhole_16S_metabolomics_matthew_samples_metadata.csv",
                                delim = ";",
                                locale=locale(encoding="latin1"))

sample_info <- read_delim("orbitrapsequence.csv", delim = ",")

sample_info <- sample_info %>%
  left_join(sample_info_extra, by = "SampleName") %>%
  select(`Sample number_Qbit`, FileName, Injection_Type) %>%
  # select(`Sample number_Qbit`, FileName, SampleName, Injection_Type, `Date `, Sinkhole, `Sampling depth (m)`,
  #        Matthew_type, Matthew_complex, Station) %>%
  separate(`Sample number_Qbit`, into = c("v1", "sampleMetag"), sep = "HL-") %>%
  #na.omit %>%
  select(FileName, sampleMetag, Injection_Type)

setwd(paste(admin, "LCMS_measured/Milou_samples/Sinkhole_microbiome", sep = "/"))
sampledata <- fread("sample_data_final.csv", colClasses = "character") %>%
  na.omit %>%
  clean_names()
colnames(sampledata)[1] <- "sampleMetag"

sampledata <- sample_info %>%
  left_join(sampledata)

samples <-sampledata %>%
  dplyr::filter(Injection_Type == "Sample")
samples <- samples$FileName


#---------Microbiome data read in---------------

setwd(paste(admin, "LCMS_measured/Milou_samples/Sinkhole_microbiome", sep = "/"))

otudata <- fread("tax_clean.csv")

microbsCswapped <- fread("OTU_swapped.csv") %>%
  select(OTU, everything()) %>%
  column_to_rownames("OTU")

microsbC <- as.data.frame(t(microbsCswapped)) %>%
  rownames_to_column("sampleMetag")

#---------------Microbiome data processing---------------

phylum_sum_pre <- microbsCswapped %>%
  rownames_to_column("OTU") %>%
  left_join(otudata %>% select(OTU, Class) %>% mutate(OTU = as.character(OTU))) %>%
  select(OTU, Class, everything())


phylum_sum <- phylum_sum_pre %>%
  column_to_rownames("OTU") %>%
  group_by(Class) %>%                        
  summarise_at(vars(-group_cols()), list(sum = sum)) %>%
  ungroup() %>%
  column_to_rownames("Class")

microb_amount <- as.data.frame(t(phylum_sum)) %>%
  rownames_to_column("sample") %>%
  separate(sample, into = c("sampleMetag"), sep = "_s")

row_sums <- rowSums(microb_amount[, -1])

row_sum_table <- as.data.frame(row_sums)

microb_amount_normalized <- microb_amount[, -1] / row_sums
#colnames(microb_amount_normalized) <- colnames(microb_amount)[-1]

microb_amount_normalized <- cbind(microb_amount[, "sampleMetag", drop = FALSE], microb_amount_normalized)


library(reshape2)
melted_data <- melt(microb_amount_normalized, id.vars = "sampleMetag")

melted_data$variable <- factor(melted_data$variable)
levels(melted_data$variable) <- c(levels(melted_data$variable), "Other")

melted_data$variable[melted_data$value < 0.1] <- "Other"



# Create a stacked bar plot
ggplot(data = melted_data, aes(x = sampleMetag, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(title = "Normalized Stacked Bar Plot",
       x = "Sample",
       y = "Normalized Value") +
  scale_fill_brewer(palette = "Set3") +  # You can choose a different color palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#find highest rank for each sample

highest_columns <- data.frame(sampleMetag = character(0), highestColumn = character(0))

# Loop through rows and find highest column
for (i in 1:nrow(microb_amount_normalized)) {
  row_data <- microb_amount_normalized[i, -1]  # Exclude the sampleMetag column
  highest_col <- names(row_data)[which.max(row_data)]
  highest_columns <- rbind(highest_columns, data.frame(sampleMetag = microb_amount_normalized$sampleMetag[i], highestColumn = highest_col))
}

# Print the new dataframe
print(highest_columns)


#--------------Metabolites list for targeted search, ready tables-------------

folder <- paste(admin, "LCMS_measured/Milou_samples/MeineLists", sep = "/")
setwd(folder)

compounds <- fread("compounds.tsv") %>%
  select(cpdID, name, formula, inchikey, smiles, charge) %>%
  left_join(fread("formulamass.tsv")) %>%
  dplyr::filter(exactMass < 1500 & exactMass > 50) %>%
  mutate(exactMass = round(exactMass, 5)) %>%
  select(cpdID, name, formula, exactMass, inchikey, smiles, charge) 


massfile <- fread("station_separated_model_metabolites_curatedmanually_mass.csv") %>%
  #left_join(pubchem, by = "cpdID") %>%
  mutate(exactMass_CID = round(exactMass_CID, digits = 5)) %>%
  mutate(adductH = exactMass_CID + 1.00728, 
         adductNa = exactMass_CID + 22.98922,
         adductM = exactMass_CID - 0.00055) %>%
  select(cpdID, adductH, adductNa, adductM) %>%
  gather(key = "adduct", value = "mz_calc", -c(cpdID)) %>%
  mutate(combine = "combine")

logP <- fread("logP_cpdID.csv")  

massFromNTS <- ms1peaks %>%
  select(id, mz, rt) %>%
  mutate(combine = "combine") %>%
  left_join(massfile, by = "combine") %>%
  mutate(Dmz = (abs(mz_calc-mz)/mz_calc)*1000000) %>%
  dplyr::filter(Dmz < 20) %>%
  left_join(logP)



#logP graph
ggplot(data = massFromNTS) +
  geom_point(mapping = aes(x = rt,
                           #colour = matched,
                           y = logP),
             size = 5, alpha = 0.7) + 
  #scale_color_viridis(option = "H") +
  facet_grid(~adduct) +
  ylab("logP, Water soluble -> Octanol soluble") +
  xlab("RT, Polar -> Non-polar") +
  my_theme

print("ok")

#--------------Getting metabolites list for targeted search-------------

#getting the manually curated list, some part later done in excel as well

folder <- paste(admin, "LCMS_measured/Milou_samples/MeineLists", sep = "/")
setwd(folder)

compounds <- fread("compounds.tsv") %>%
  select(smiles) %>%
  subset(nchar(smiles) > 1)

#write_delim(compounds, "compounds_smiles.csv", delim = ";")

#all metabolites and their chemical info
compounds <- fread("compounds.tsv") %>%
  #left_join(fread("formulamass.tsv")) %>%
  #dplyr::filter(exactMass < 1500 & exactMass > 50) %>%
  #mutate(exactMass = round(exactMass, 5)) %>%
  select(cpdID, formula, inchikey, smiles, charge) #%>% #exactMass
  #dplyr::filter(charge == 0)
  #select(exactMass) %>% unique()

#PubChem lists

pubchem <- compounds %>%
  left_join(fread("PubChem/all_CIDorigsmiles.txt"), by = "smiles") %>%
  left_join(fread("PubChem/all_inchikey.txt"), by = "smiles") %>%
  left_join(fread("PubChem/all_CIDsmiles.txt"), by = "smiles") %>%
  left_join(fread("PubChem/CID_mass_calculated.csv")) %>%
  unique() %>%
  select(cpdID, smiles, exactMass_CID, formula_CID)

# CIDmass <- fread("PubChem/all_CIDorigsmiles.txt") %>%
#   left_join(fread("PubChem/all_CIDsmiles.txt"), by = "smiles") %>%
#   na.omit() %>%
#   mutate(formula_CID = "form",
#          exactMass_CID = "mass")
# 
# for (n in 1:length(CIDmass$smiles_CID)) {
#   CIDmass$exactMass_CID[[n]] = MassFromSmiles(CIDmass$smiles_CID[[n]])
# }



logP_table <- tibble()
for(n in compounds$smiles){

  logP <- tryCatch({
    fn_logP(n)
  }, error = function(e){
    return(NULL)
  })
  
  if(!is.null(logP)){
    logP <- as.data.frame(logP) %>%
    mutate(SMILES = n)
  }
  
  logP_table <- logP_table %>%
    bind_rows(logP)
}



#Meine smaller lists
predcompunds <- fread("station_separated_model_metabolites.tsv") %>%
  dplyr::select(-c(station, pH, depth))
predcompunds <- gather(predcompunds)

selectedcompounds <- tibble()
for(n in predcompunds$value){
  id <- gsub("[{}']", "", n)
  id <- gsub(" ", "", id)
  id <- unlist(strsplit(id, ","))
  
  listtable <- as.data.frame(id)
  selectedcompounds <- selectedcompounds %>%
    bind_rows(listtable) %>%
    unique()
}


selectedcompounds <- selectedcompounds %>%
  left_join(compounds) %>%
  dplyr::filter(exactMass < 1500 & exactMass > 50)



predmasses <- fread("station_separated_model_metabolites_curatedmanually_mass.csv") %>%
  left_join(fread("PubChem/all_inchikey.txt"), by = c("smiles")) %>%
  mutate(InChIKey2D = InChIKey_CID) %>%
  separate(InChIKey2D, into = c("InChIkey2D"), sep = "-") %>%
  unique() %>%
  select(cpdID, InChIkey2D)

  
print("ok")


#---------Reading in LCMS, MZmine files------------------

folder <- paste(admin, "LCMS_measured/Milou_samples/MSmine_results/20240208", sep = "/")
setwd(folder)

#CSV tables

files <- dir(folder, pattern = ".csv")
files

ms1peaks <- fread(files[1]) %>%
  mutate(peakcheck = height/area) %>%
  select(id,rt,mz,height,area,peakcheck,contains("area")) %>%
  dplyr::filter(peakcheck > 1)

ms2peaks <- fread(files[2]) %>%
  mutate(peakcheck = height/area) %>%
  select(id,rt,mz,height,area,peakcheck,contains("area")) %>%
  #filter(peakcheck > 1) %>%
  mutate(mzgroup = "none") %>%
  select(id,mz,mzgroup,contains(":area"))

for(n in 1:length(ms2peaks$id)){
  if(ms2peaks$mz[[n]] <= 300){
    ms2peaks$mzgroup[[n]] = "massgr150-300"
  }
  if(ms2peaks$mz[[n]] <= 400 & ms2peaks$mz[[n]] > 300){
    ms2peaks$mzgroup[[n]] = "massgr300-400"
  }
  if(ms2peaks$mz[[n]] <= 700 & ms2peaks$mz[[n]] > 400){
    ms2peaks$mzgroup[[n]] = "massgr400-700"
  }
}


# MGF sirius files

files <- dir(paste0(folder, "/smallerGroups"), pattern = ".mgf")
setwd(paste0(folder, "/smallerGroups"))
files
filename <- files[1]

combined_peaks <- tibble()

for(filename in files){
  mgf_file_id <- read_delim(filename, delim = ";", col_names = FALSE) %>%
    dplyr::filter(grepl("FEATURE_ID=", X1)) %>%
    dplyr::rename(id = X1)
  
  mgf_file_ms2 <- read_delim(filename, delim = ";", col_names = FALSE) %>%
    dplyr::filter(grepl("MSLEVEL=", X1)) %>%
    dplyr::rename(mslevel = X1)
  
  mgf_file_ret <- read_delim(filename, delim = ";", col_names = FALSE) %>%
    dplyr::filter(grepl("RTINSECONDS=", X1)) %>%
    dplyr::rename(RT = X1)
  
  mgf_file_mass <- read_delim(filename, delim = ";", col_names = FALSE) %>%
    dplyr::filter(grepl("PEPMASS=", X1)) %>%
    dplyr::rename(mass = X1)
  
  mgf_file <- cbind(mgf_file_id, mgf_file_ms2, mgf_file_ret, mgf_file_mass) %>%
    unique() %>%
    dplyr::filter(mslevel == "MSLEVEL=2")
  
  combined_peaks <- combined_peaks %>%
    bind_rows(mgf_file)
  
}



#----------------Analysing MZmine tables---------------------

dataname <- "Aligned, bothdays, MS2, all features"
sample_feature <- ms2peaks %>%
  select(id, contains(":area")) %>%
  column_to_rownames("id") %>%
  replace(is.na(.), 0) %>%
  mutate_all(~ ifelse(. != 0, log10(.), 0))

sample_feature <- as.data.frame(t(sample_feature))

#eliminating zero rows and column and continuing with clustering
constant_cols <- apply(sample_feature, 2, function(col) length(unique(col)) == 1)
sample_feature <- sample_feature[, !constant_cols]

zero_cols <- apply(sample_feature, 2, function(col) all(col == 0))
sample_feature <- sample_feature[, !zero_cols]

set.seed(123)
data_pca <- prcomp(sample_feature, rank = 3)
PC1var <- round(summary(data_pca)$importance[2,1]*100,1)
PC2var <- round(summary(data_pca)$importance[2,2]*100,1)

pca_fp <- as.data.frame(data_pca$x) %>%
  rownames_to_column("FileName")
pca_fp$FileName <- str_extract(pca_fp$FileName, "(?<=datafile:)[^.]+")  
# graphdata <- clustering_p %>%
#   left_join(sample_info) %>%
#   left_join(sampledata)


#------------Metabolomics results from SIRIUS, averaged--------------

setwd(paste(admin, "LCMS_measured/Milou_samples/CalculatedFingerprints", sep = "/"))
dir(paste(admin, "LCMS_measured/Milou_samples/CalculatedFingerprints", sep = "/"))

filename <- "20240208_CANOPUS_aligned_bintres05_aver_INTENSITYnorm.csv"

all_aver_classif <- fread(filename) %>%
  #right_join(highph) %>%
  column_to_rownames("sample") %>%
  na.omit()

# all_aver_classif <- averages_msgroup_ %>% #averages_msgroup_
#   column_to_rownames("sample")

#not from .csv files
all_aver_classif <- averages %>%
  dplyr::filter(sample %in% samples) %>%
  column_to_rownames("sample")


#eliminating zero rows and column and continuing with clustering
constant_cols <- apply(all_aver_classif, 2, function(col) length(unique(col)) == 1)
all_aver_classif <- all_aver_classif[, !constant_cols]

zero_cols <- apply(all_aver_classif, 2, function(col) all(col == 0))
all_aver_classif <- all_aver_classif[, !zero_cols]



#pca and umap
set.seed(123)
pca_fp <- prcomp(all_aver_classif, rank = 3)

PC1var <- round(summary(pca_fp)$importance[2,1]*100,1)
PC2var <- round(summary(pca_fp)$importance[2,2]*100,1)

pca_fp <- as.data.frame(pca_fp$x) %>%
  rownames_to_column("FileName")

#colnames(aver_fp)[1] <- "FileName"

#UMAP
set.seed(123)
umap <- umap(all_aver_classif)  #n_components = 3

pca_fp <- as.data.frame(umap$layout) %>%
  rownames_to_column("FileName")

#-----------Graphs for comparing samples------------

graphdata <- sampledata %>%
  full_join(pca_fp, by = "FileName") %>%
  mutate(silicate_umol_kg = as.numeric(gsub(",", ".", silicate_umol_kg))) %>%
  mutate(p_h = as.numeric(gsub(",", ".", p_h))) %>%
  mutate(depth_m = as.numeric(gsub(",", ".", depth_m))) %>%
  #drop_na(PC1) %>%
  mutate(date = substr(FileName, 1, 5))

# highph <- graphdata %>%
#   filter(p_h > 7 | is.na(p_h)) %>%
#   select(FileName)
# colnames(highph) <- "sample"

#colnames(graphdata)

ggplot(data = graphdata) +
  # geom_point(mapping = aes(x = PC1,
  #                          colour = p_h, #Injection_Type  , Matthew_type, `Sampling depth (m)`
  #                          y = PC2),
  #            size = 8, alpha = 1) +
  geom_point(mapping = aes(x = PC1,
                           shape = date,
                           y = PC2),
             size = 6, alpha = 0.5, color = "gray", data = graphdata[is.na(graphdata$depth_m),]) + # NA points
  geom_point(mapping = aes(x = PC1,
                           colour = p_h,
                           shape = date,
                           y = PC2),
             size = 6, alpha = 1, data = graphdata[!is.na(graphdata$p_h),]) + # Non-NA points
  my_theme +
  scale_color_viridis(option = "H") +
  ggtitle("Aligned, bothdays, CANOPUS average intensities, all LCMS features") +
  #facet_wrap(~ Injection_Type) +
  xlab(paste("PC1, ", PC1var, "%", sep = "")) +
  ylab(paste("PC2, ", PC2var, "%", sep = ""))

ggplot(data = graphdata) +
  geom_point(mapping = aes(x = V1,
                           shape = date,
                           y = V2),
             size = 6, alpha = 0.5, color = "gray", data = graphdata[is.na(graphdata$depth_m),]) + # NA points
  geom_point(mapping = aes(x = V1,
                           colour = p_h,
                           shape = date,
                           y = V2),
             size = 6, alpha = 1, data = graphdata[!is.na(graphdata$p_h),]) + # Non-NA points
  my_theme +
  scale_color_viridis(option = "H") +
  #facet_wrap(~ Injection_Type)
  ggtitle("Aligned, bothdays, CANOPUS average intensities, all LCMS features") +
  xlab(paste("UMAP1", sep = "")) +
  ylab(paste("UMAP2", sep = ""))


PC_aver


#---------------Chemicals distribution in sample-------------

folder <- paste(admin, "LCMS_measured/Milou_samples/MSmine_results/FPtables", sep = "/")
setwd(folder)

files <- dir(folder, pattern = ".csv")

sample_canopus <- fread("CANOPUS_O2006085_20230503.csv")

sample_canopus_pca <- sample_canopus %>%
  select(id, matches("Canop")) %>%
  column_to_rownames("id")

pca <- prcomp(sample_canopus_pca, rank = 3)
#pca <- umap(sample_canopus_pca)

PC1var <- round(summary(pca)$importance[2,1]*100,1)
PC2var <- round(summary(pca)$importance[2,2]*100,1)

#pca_table <- as.data.frame(pca$layout) %>%
pca_table <- as.data.frame(pca$x) %>%
  rownames_to_column("id") %>%
  left_join(sample_canopus %>%
              select(id, rt, mz) %>%
              mutate(id = as.character(id)))


ggplot(data = pca_table) +
  geom_point(mapping = aes(x = PC1,
                           colour = mz,
                           y = PC2),
             size = 5, alpha = 0.7) + 
  scale_color_viridis(option = "H") +
  #facet_grid(~organic_modifier_percentage) +
  xlab(paste("PC1, ", PC1var, "%", sep = "")) +
  ylab(paste("PC2, ", PC2var, "%", sep = "")) +
  my_theme

#----------Important features for separating sinkholes---------------

sample_feature_cluster <- sample_feature %>%
  rownames_to_column("FileName")
sample_feature_cluster$FileName <- str_extract(sample_feature_cluster$FileName, "(?<=datafile:)[^.]+")

cluster <- sample_info %>%
  na.omit() %>%
  dplyr::filter(Injection_Type == "Sample") %>%
  left_join(sampledata %>% select(sampleMetag, p_h)) %>%
  na.omit() %>%
  mutate(p_h = as.numeric(gsub(",", ".", p_h))) %>%
  mutate(cluster = "1") 

for(n in 1:length(cluster$FileName)){
  if(cluster$p_h[n] < 7){
    cluster$cluster[n] = 2
  }
}

sample_feature_cluster <- cluster %>%
  select(FileName, cluster) %>%
  left_join(sample_feature_cluster)


#finding variables
classes <- as.factor(sample_feature_cluster$cluster)
df <- sample_feature_cluster %>%
  select(-cluster, -FileName)
#Continue in code 20231103_RFE_Aris.R


#working with selected LCMS features
folder <- paste(admin, "LCMS_measured/Milou_samples/MSmine_results/20240208", sep = "/")
setwd(folder)
importantfeatures <- fread("20240208_bothdays_MS2_selectedfeatures_Allbest100for10fold.txt") %>%
  #na.omit() %>%
  mutate(FoldMean = rowMeans(select(., contains("Fold")), na.rm = TRUE)) %>%
  mutate(feature = as.character(feature))

feature_analysis <- sample_feature %>%
  select(one_of(importantfeatures$feature)) %>%
  #select(one_of(targets_fromNTS$id.y)) %>%
  rownames_to_column("FileName")
feature_analysis$FileName <- str_extract(feature_analysis$FileName, "(?<=datafile:)[^.]+")
feature_analysis <- feature_analysis %>%
  gather(key = "feature", value = "logInt", -FileName) %>%
  #left_join(graphdata) %>%
  #mutate(p_h = as.numeric(gsub(",", ".", p_h))) %>%
  left_join(cluster) %>%
  na.omit() %>%
  left_join(importantfeatures %>% select(feature, FoldMean))

#Clustering similar features together based on sample similarity intesities

#clustering similar features together
feature_analysis_clusters <- feature_analysis %>%
  select(FileName, feature, logInt) %>%
  spread(key = FileName, value = logInt) %>%
  column_to_rownames("feature")

feature_analysis_clusters <- as.matrix(feature_analysis_clusters)
k <- 7  # Number of clusters
kmeans_result <- kmeans(feature_analysis_clusters, centers = k)
cluster_assignment <- as.data.frame(kmeans_result$cluster) %>%
  rownames_to_column("feature")
colnames(cluster_assignment)[2] <- "FeatureCluster"




#Graphs!
var_plot <- feature_analysis %>%
  #left_join(targets_fromNTS %>% select(feature, name)) %>%
  left_join(cluster_assignment, by = "feature") %>%
  mutate(FoldMean = round(FoldMean, digits=2)) %>%
  #dplyr::filter(FoldMean > 2) %>%
  dplyr::filter(feature %in% c(653,897,391)) %>%
  #dplyr::filter(FeatureCluster == 5) %>%
  ggplot(aes(x = factor(cluster), 
             y = logInt)) +
  geom_point(aes(x = cluster, y = logInt, colour = p_h),
            size = 6,
            position = position_jitter(width = 0.2)) +
  # geom_violin(aes(fill = cluster), 
  #             position = "dodge") +  
  # geom_text(aes(x = methodnr,
  #               label = methodname,
  #               y = value), size = 3, angle = 30) +
  my_theme +
  scale_color_viridis(option = "H") +
  #theme_light() +
  #scale_color_manual(values=color_vector_empo_womix) + 
  #ggtitle(modelname) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~ feature+FoldMean, scales = "free_y") + #scales = "free_y"
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1))
var_plot


#Violin
var_plot <- feature_analysis %>%
  ggplot(aes(x = cluster, 
             y = logInt)) +
  geom_violin(aes(fill = cluster), 
              scale = "width", width = 0.7) +
  facet_wrap(~ feature) +
  theme(legend.position="none") +
  #ggtitle("Cluster one") + 
  #scale_fill_manual(values = color_vector_empo_womix_new) + 
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  #theme_light()
  my_theme
var_plot

#----------Important SIRIUS classes for separating sinkholes---------------

cluster <- sampledata %>%
  na.omit() %>%
  dplyr::filter(Injection_Type == "Sample") %>%
  na.omit() %>%
  mutate(p_h = as.numeric(gsub(",", ".", p_h))) %>%
  mutate(cluster = "1")

for(n in 1:length(cluster$FileName)){
  if(cluster$p_h[n] < 7){
    cluster$cluster[n] = 2
  }
}

cluster <- cluster %>%
  mutate(sample = FileName) %>%
  select(sample, cluster) %>%
  left_join(all_aver_classif %>% rownames_to_column("sample"))


# classes <- as.factor(cluster$cluster)
# df <- cluster %>%
#   select(-cluster, -sample)

selectedfeatures <- fread("20240208_rfe_CANOPUS_aligned_bintres05_aver_INTENSITYnorm_Allbest100for10fold.csv") %>%
  na.omit()
selectedfeatures <- selectedfeatures$var

feature_analysis <- cluster %>%
  gather(key = "feature", value = "average", -c(cluster, sample)) %>%
  dplyr::filter(feature %in% selectedfeatures) %>%
  mutate(FileName = sample)  %>%
  left_join(sampledata %>% select(FileName, p_h)) %>%
  mutate(p_h = as.numeric(gsub(",", ".", p_h))) %>%
  separate(feature, into = c("Canop"), sep = "_") %>%
  left_join(sir_canop_pos_Kai %>% select(Canop, name)) %>%
  group_by(Canop) %>%
  mutate(Canop_aver = mean(average)) %>%
  ungroup()


var_plot <- feature_analysis %>%
  #dplyr::filter(Canop %in% c("Canop1238", "Canop118", "Canop1994", "Canop3450", "Canop3409", "Canop36", "Canop3619")) %>%
  ggplot(aes(x = factor(cluster), 
             y = average)) +
  geom_point(aes(x = cluster, y = average, color = p_h),
             size = 6,
             position = position_jitter(width = 0.2)) +
  my_theme +
  scale_color_viridis(option = "H") +
  #theme_light() +
  ggtitle("CANOPUS_aligned_bintres05_aver_all") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~ Canop+name, scales = "free_y") + #scales = "free_y"
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1))
var_plot


#samples chemical space
Canop_var_decreasing <- feature_analysis %>%
  select(Canop, Canop_aver) %>%
  unique() %>%
  arrange(desc(Canop_aver)) %>%
  mutate(CanVarDesc = row_number())

feature_analysis_ <- feature_analysis %>%
  left_join(Canop_var_decreasing) 

var_plot <- feature_analysis_ %>%
  #dplyr::filter(Canop %in% c("Canop1238", "Canop118", "Canop1994", "Canop3450", "Canop3409", "Canop36", "Canop3619")) %>%
  ggplot(aes(x = CanVarDesc,
             y = average)) +
  geom_point(aes(x = CanVarDesc, y = average, color = p_h),
             size = 6,
             position = position_jitter(width = 0.2)) +
  my_theme +
  scale_color_viridis(option = "H") +
  #theme_light() +
  ggtitle("CANOPUS_aligned_bintres05_aver_all") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~ cluster) + #scales = "free_y"
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1))
var_plot


print("ok")

#--------------Cutting .mgf files into smaller--------------

files <- dir(folder, pattern = ".mgf")
files
filename <- files[1]

mgf_file_id <- read_delim(filename, delim = ";", col_names = FALSE) %>%
  dplyr::filter(grepl("FEATURE_ID=", X1)) %>%
  dplyr::rename(id = X1)

mgf_file_ms2 <- read_delim(filename, delim = ";", col_names = FALSE) %>%
  dplyr::filter(grepl("MSLEVEL=", X1)) %>%
  dplyr::rename(mslevel = X1)

mgf_file_ret <- read_delim(filename, delim = ";", col_names = FALSE) %>%
  dplyr::filter(grepl("RTINSECONDS=", X1)) %>%
  dplyr::rename(RT = X1)

mgf_file_mass <- read_delim(filename, delim = ";", col_names = FALSE) %>%
  dplyr::filter(grepl("PEPMASS=", X1)) %>%
  dplyr::rename(mass = X1)

mgf_file <- cbind(mgf_file_id, mgf_file_ms2, mgf_file_ret, mgf_file_mass) %>%
  unique() %>%
  dplyr::filter(mslevel == "MSLEVEL=2")



#make smaller files, in five groups

input_file <- "20240208_sinkhole_aligned_bothdays_blanksub_MS2.mgf"
output_directory <- "output_files/"

# Create the output directory if it doesn't exist
dir.create(output_directory, showWarnings = FALSE)

# Read the input .mgf file
mgf_content <- readLines(input_file)

# Find the index of the row containing "FEATURE_ID=1292"
start_index <- which(mgf_content == "FEATURE_ID=11")
end_index <- which(mgf_content == "FEATURE_ID=44")

# If "FEATURE_ID=1292" is found, split the content and write out the new .mgf file
if (length(end_index) > 0) {
  # Extract content before the "FEATURE_ID=1292" row
  new_mgf_content <- mgf_content[(start_index[1] - 1):(end_index[1] - 2)]
  
  # Define the path for the output .mgf file
  output_file <- paste0(output_directory, "test.mgf")
  
  # Write out the new .mgf file
  writeLines(new_mgf_content, con = output_file)
  
  # Print message indicating the process completion
  cat("Output file", output_file, "has been created.\n")
} else {
  # If "FEATURE_ID=1292" is not found, print a message
  cat("Row 'FEATURE_ID=1292' not found in the input file.\n")
}

#-----------Getting averaged fp from tables---------------

fppath <- paste(admin, "LCMS_measured/Milou_samples/MSmine_results/FPtables", sep = "/")
setwd(fppath)
files <- dir(fppath, pattern = "CANOPUS_20240208")

all_canopus <- tibble()
for (file in files){
  canopus <- fread(file)
  all_canopus <- all_canopus %>%
    bind_rows(canopus)
}


ms2peakspresence <- ms2peaks %>%
  select(id, mzgroup, contains("datafile:"))
ms2peakspresence <- gather(ms2peakspresence, key = "sample", value = "intensity", -id, -mzgroup) %>%
  na.omit()

ms2peakspresence$sample <- str_extract(ms2peakspresence$sample, "(?<=datafile:)[^.]+")

ms2peakspresence <- ms2peakspresence %>%
  left_join(all_canopus %>% select(id, contains("Canop")), 
            by = "id") %>%
  na.omit() %>%
  group_by(sample) %>%
  mutate(intsum = sum(intensity)) %>%
  ungroup() %>%
  mutate(intnorm = intensity/intsum) %>%
  select(id, sample,intensity, intsum, intnorm, everything())

fp_treshold = 0.5
# ms2peakspresence[ms2peakspresence < fp_treshold & ms2peakspresence > 0] <- 0
# ms2peakspresence[ms2peakspresence >= fp_treshold & ms2peakspresence < 1] <- 1
canop_cols <- grep("^Canop", names(ms2peakspresence), value = TRUE)
# Loop through each "Canop" column and apply the thresholding
for (col in canop_cols) {
  ms2peakspresence[col][ms2peakspresence[col] < fp_treshold & ms2peakspresence[col] > 0] <- 0
  ms2peakspresence[col][ms2peakspresence[col] >= fp_treshold & ms2peakspresence[col] <= 1] <- 1
}


# test <- head(ms2peakspresence, 10) %>%
#   mutate(intsum = sum(intensity)) %>%
#   mutate(intnorm = intensity/intsum)

# get table with intensities 
ms2peaksintens <- ms2peakspresence %>%
  mutate_at(vars(starts_with("Canop")), ~. * intnorm) %>%
  select(id, mzgroup, sample, starts_with("Canop"))
#average by summing the scaled intensities
averages <- ms2peaksintens %>%
  select(sample, contains("Canop")) %>%
  group_by(sample) %>%                        
  summarise_at(vars(contains("Canop")),
               list(av = sum)) %>%
  ungroup()




#just averages
averages <- ms2peaksintens %>%
  select(sample, contains("Canop")) %>%
  group_by(sample) %>%                        
  summarise_at(vars(contains("Canop")),
               list(av = mean)) %>%
  ungroup()


#massgroup averages
samplelist <- unique(ms2peakspresence$sample)

averages_msgroup <- tibble()
for(samplename in samplelist){
  
  fp_table <- ms2peakspresence %>%
    filter(sample == samplename) %>%
    na.omit()
  
    average_fingerprints <- fp_table %>%
    group_by(mzgroup) %>%                        
    summarise_at(vars(contains("Canop")),
                 list(av = mean)) %>%
    ungroup() %>%
    column_to_rownames("mzgroup")
  
  average_fingerprints <- as.matrix(average_fingerprints)
  
  sample_vector <- as.data.frame(t(c(average_fingerprints)))
  rownames(sample_vector) <- c(samplename)
  
  averages_msgroup <- averages_msgroup %>%
    bind_rows(sample_vector)
}

avfp_names <- colnames(average_fingerprints)
additions <- c("100to300", "300to400", "400to700")
msgr_colnames <- character()

for (original_element in avfp_names) {
  for (addition in additions) {
    modified_element <- paste0(original_element, addition)
    msgr_colnames <- c(msgr_colnames, modified_element)
  }
}

averages_msgroup_ <- averages_msgroup
colnames(averages_msgroup_) <- msgr_colnames
averages_msgroup_ <- rownames_to_column(averages_msgroup_, "sample")


#---------PCA for all the features, to find similar----------------

features_per_sample <- ms2peakspresence %>%
  select(id,contains("CANOP")) %>%
  unique() %>%
  rownames_to_column("v1") %>%
  select(-v1) %>%
  column_to_rownames("id")

#pca and umap
set.seed(123)
pca_fp <- prcomp(features_per_sample, rank = 3)

PC1var <- round(summary(pca_fp)$importance[2,1]*100,1)
PC2var <- round(summary(pca_fp)$importance[2,2]*100,1)

pca_fp_table <- as.data.frame(pca_fp$x) %>%
  rownames_to_column("id") %>%
  left_join(ms2peakspresence %>% select(id, sample) %>% mutate(id = as.character(id))) %>%
  mutate(FileName = sample) %>%
  left_join(graphdata %>% select(FileName, p_h, date)) %>%
  na.omit() %>%
  dplyr::filter(FileName %in% c("O2006040", "O2006052", "O2006100", "O2006076"))
  #dplyr::filter(id %in% importantfeatures$feature)

ggplot(data = pca_fp_table) +
  geom_point(mapping = aes(x = PC1,
                           colour = p_h,
                           y = PC2),
             size = 5, alpha = 0.8) + 
  scale_color_viridis(option = "H") +
  facet_wrap(~sample, nrow = 1) +
  xlab(paste("PC1, ", PC1var, "%", sep = "")) +
  ylab(paste("PC2, ", PC2var, "%", sep = "")) +
  my_theme


#---------SIRIUS annotation tables----------------

sir_folder <- "C:/Users/pille/Desktop/MZmine_results/20240208"
setwd(sir_folder)


#Rank 1 structures only
comp_ind_files <- dir(sir_folder, pattern = "compound_identifications", recursive = TRUE)
file <- comp_ind_files[3]

top1annotations <- tibble()
for(file in comp_ind_files){
  id_file <- fread(file) %>%
    mutate(ConfidenceScore = as.numeric(ConfidenceScore))

  top1annotations <- top1annotations %>%
    bind_rows(id_file)
}

top1annotations <- top1annotations %>%
  select(featureId,smiles,name, InChIkey2D, ConfidenceScore, molecularFormula, adduct) %>%
  # right_join(importantfeatures %>% 
  #              mutate(featureId = as.numeric(feature)) %>%
  #              select(featureId, FoldMean)) %>%
  na.omit(name)


all_annotations_targets <- top1annotations %>%
  left_join(predmasses, by = "InChIkey2D")



#getting all annotated structures, filter maybe to top50

folderdata <- "C:/Users/pille/Desktop/MZmine_results/20240208"

selected_structures <- SiriusStructures(folderdata, 
                                        subfoldercommand = "sinkhole")


all_sirius_structures_top <- selected_structures %>%
  dplyr::filter(structureRankPerFormula <= 50)


targets_fromNTS <- predmasses %>%
  inner_join(all_sirius_structures_top, by = "InChIkey2D") %>%
  select(id, cpdID, InChIkey2D) %>%
  left_join(compounds %>% select(cpdID, name))

