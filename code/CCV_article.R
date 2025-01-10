#Pilleriin Peets (pilleriin.peets@gmail.com)

#codes for Chemical characteristics vectors from SIRIUS calculated data

#Data in Zenodo connected to the article
#Zenodo doi: 

library(tidyverse)
library(data.table)
library(umap)
library(caret)
library(scales)
library(fpc)
library(rlist)
library(Rdisop)
library(tibble)
library(rcdklibs)
library(rcdk)

#--------Biome colors--------------
color_vector_empo_womix_new <- c("#9467bd","#d62728","#e377c2","#ff7f0e",
                                 "#FFCC33","#bcbd22","#2ca02c","#CC9966",
                                 "#660000","#999999","#3300CC")

#--------Function to get MFP and CC for samples from SIRIUS5 calculated files-----------

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


#--------Calculating average CCV from SIRIUS5 tables---------------

#MF from data where all samples and fingerprints are in one table

MFdata = characteristics_table #table with columns: sample + rank1 selected SIRIUS characteristics from MFtable_SIRIUS5
CCdata = characteristics_table #table with columns: sample + rank1 selected SIRIUS characteristics from CCtable_SIRIUS5

#Probabilistic SIRIUS values into binary
fp_treshold = 0.5
MFdata[MFdata < fp_treshold & MFdata > 0] <- 0
MFdata[MFdata >= fp_treshold & MFdata < 1] <- 1
CCdata[CCdata < fp_treshold & CCdata > 0] <- 0
CCdata[CCdata >= fp_treshold & CCdata < 1] <- 1


average_MF <- MFdata %>%
  group_by(sample) %>%                        
  summarise_at(vars(contains("Un")),
               list(av = mean)) %>%
  ungroup()

average_CC <- CCdata %>%
  group_by(sample) %>%                        
  summarise_at(vars(contains("Canop")),
               list(av = mean)) %>%
  ungroup()

#for extensive amount of data, use loops to do the averaging one by one

samplelist <- average_MF %>%
  select(sample) %>%
  unique()
samplelist <- as.vector(samplelist$sample)

average_MF <- tibble()
for(samplename in samplelist){
  table <- MFdata %>%
    dplyr::filter(sample == samplename) %>%
    na.omit()
  
  average_fingerprint <- table %>%
    group_by(sample) %>%                        
    summarise_at(vars(contains("Un")),
                 list(av = mean)) %>%
    ungroup()
  
  average_MF <- average_MF %>%
    bind_rows(average_fingerprint)
}


#--------PCA, UMAP, metrics for data-----------------

#add correct folder
setwd("C:/Users/pille/OneDrive - Friedrich-Schiller-Universität Jena/JenaUniversity/Manuscripts/EMP_article/Data_CCV")
file <- "SI2_MF_average_all_variables.csv"

data <- fread(file)

data_for_clustering <- data %>%
  select(sample, contains("_av")) %>%
  column_to_rownames("sample")

#PCA
set.seed(123)
data_pca <- prcomp(data_for_clustering, rank = 2)

plot_data <- as.data.frame(data_pca$x) %>%
  rownames_to_column("sample") %>%
  left_join(data %>% select(sample, empo_4_type))
colnames(plot_data) <- c("sample","Dim1", "Dim2", "Biome")

#UMAP
set.seed(123)
umap <- umap(data_for_clustering) 
plot_data <- as.data.frame(umap$layout) %>%
  rownames_to_column("sample") %>%
  left_join(data %>% select(sample, empo_4_type))
colnames(plot_data) <- c("sample","Dim1", "Dim2", "Biome")

#plot

plot <- ggplot(data = plot_data) +
  geom_point(mapping = aes(x = Dim1,
                           colour = biome,
                           y =Dim2),
             size = 3,
             alpha = 0.7) +
  theme_light() +
  scale_color_manual(values=color_vector_empo_womix_new)

plot


#CH, Dunn and Silhouette indeces

biomes_label <- data %>%
  select(empo_4_type) %>%
  unique() %>%
  arrange(empo_4_type) %>%
  mutate(lable = row_number())

sample_for_distance <- data %>%
  mutate(shortsample = paste(empo_4_type, row_number(), sep = "_")) %>%
  select(-sample, -empo_4_type) %>%
  column_to_rownames("shortsample") 

dist_matrix <- dist(sample_for_distance, method = "euclidean")
similarity_matrix <- as.matrix(dist_matrix)

similarity_matrix_scaled <- as.data.frame(similarity_matrix) %>%
  rownames_to_column("sample1") %>%
  gather(key = "sample2", value = "dist_aver", -sample1) %>%
  mutate(dist_aver = rescale(dist_aver, to = c(0, 1))) %>%
  spread(key = sample2, value = dist_aver) %>%
  column_to_rownames("sample1")

similarity_matrix_scaled <- as.matrix(similarity_matrix_scaled)

cluster <- as.data.frame(similarity_matrix_scaled) %>%
  rownames_to_column("sample") %>%
  select(sample) %>%
  separate(sample, into = c("empo_4_type", "sample"), sep = "_") %>%
  left_join(biomes_label)

metrics <- cluster.stats(similarity_matrix_scaled, cluster$lable)

siluettwidth <- metrics$avg.silwidth
dunn <- metrics$dunn
ch <- metrics$ch



#--------RFE for selecting important features-----------------

#Supplementary files SI2-SI13 have all tested CCV data
file <- "SI2_MF_average_all_variables.csv"

data <- fread(file)

biomes_label <- data %>%
  select(empo_4_type) %>%
  unique() %>%
  arrange(empo_4_type) %>%
  mutate(class = row_number())

data_rfe <- data %>%
  left_join(biomes_label)

df <- data_rfe %>%
  select(-empo_4_type, -class, -sample)

classes <- as.factor(data_rfe$class)

#Recursive Feature Elimination (RFE)

set.seed(123)
ctlr=rfeControl(rfFuncs,method = 'cv',number = 10)
varsize <- c(1:9, seq(10, 100, by = 10), seq(200, 1000, by = 100), seq(2000, 5000, by = 1000)) #change according to data
varsize <- c(1,10,100,1000)

set.seed(123)
results = rfe(x = df , y = classes , rfeControl = ctlr, sizes = varsize, verbose=T)
plot(results, type = "o", main = "averages")


bestvariables <- results$variables %>%
  dplyr::filter(Variables %in% c(10)) %>%  #choose suitable N
  select(Overall, var, Resample, Variables) %>%
  unique()

bestvariables_spread <- bestvariables %>%
  spread(key = Resample, value = Overall) %>%
  na.omit()

#------------END------------------