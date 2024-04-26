library(tidyverse)
library(rjson)
library(janitor)
library(data.table)
library(umap)
library(plotly)
library(cluster)
library(factoextra)
library(caret)
library(caTools)
library(shapr)
library(scatterplot3d)
library(reshape2)
library(Rtsne)
library(ape)
library(cluster)
library(factoextra)
library(fpc)
library(proxy)
library(scales)
library(corrr)


#----------Functions-----------
admin <- "C:/Users/pille/OneDrive - Friedrich-Schiller-Universität Jena/JenaUniversity"
setwd(paste(admin, "RcodesJena", sep = "/"))
source("Graph/my_theme_1.R")
#source("SU_codes/Function_fingerprint_calculations.R")

color_vector_empo_womix_new <- c("#9467bd","#d62728","#e377c2","#ff7f0e",
                                 "#FFCC33","#bcbd22","#2ca02c","#CC9966",
                                 "#660000","#999999","#3300CC")


# color_vector <- c("#17becf", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#3300CC",
#                   "#e377c2", "#7f7f7f", "#bcbd22", "#1f77b4", "#FFCC33", "#9c4a5a", "#000000")
# 
# color_vector_bar <- c("#17becf", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#3300CC",
#                       "#e377c2", "#7f7f7f", "#bcbd22", "#1f77b4", "#FFCC33", "#9c4a5a", "#000000")
# 
# color_vector_empo <- c("#17becf", "light grey", "#ff7f0e", "#d62728", "#9c4a5a", "#2ca02c", "#1f77b4", 
#                        "#3300CC", "#e377c2","#FFCC33" , "#bcbd22", "light grey", "#9467bd")
# 
# color_vector_empo_womix <- c("#17becf", "#ff7f0e", "#d62728", "#9c4a5a", "#2ca02c", "#1f77b4", 
#                              "#3300CC", "#e377c2","#FFCC33" , "#bcbd22", "#9467bd")


#-------------SIRIUS fingerprints and CANOPUS----------

setwd(paste(admin, "Datasets_fp/GNPS_Kai", sep = "/"))
sir_fp_pos_Kai <- fread("csi_fingerid.tsv") %>%
  mutate(Un = paste("Un", absoluteIndex, sep = ""))
sir_canop_pos_Kai <- fread("canopus.tsv") %>%
  mutate(Canop = paste("Canop", absoluteIndex, sep = ""))


#----------Metadata from github and SI for metabolomics analysis, selected data---------------
gnpsid <- "MSV000083475"
setwd(paste(admin, "Pairedomics_projects",gnpsid, sep = "/"))

#type_provider <- fread("41564_2022_1266_MOESM3_ESM_type.csv") 

sample_type <- fread("sample_type_information_extra.csv") %>%
  left_join(fread("41564_2022_1266_MOESM3_ESM_type.csv"), 
            by = "emp500_study_id") %>%
  dplyr::filter(!is.na(emp500_study_id)) %>%
  left_join(fread("samplenames_corrected.csv")) %>%
  na.omit("sample") %>%
  select(sample, samples, everything())

path_calc_fp <- paste(admin, "Datasets_fp/GNPS_Kai/CalculatedFingerprints", sep = "/")
setwd(path_calc_fp)

names <- fread("samplenames_corrected.csv")

dir(path_calc_fp, pattern = "rfe")
selected_features <- fread("20240117_rfe_canopaver_Allbest100for10fold.csv") %>%
  na.omit() %>%
  mutate(var = as.character(var)) %>%
  mutate(FoldMean = rowMeans(select(., contains("Fold")), na.rm = TRUE))
selected_features_vector <- selected_features$var

selected_features_described <- selected_features %>%
  select(var) %>%
  separate(var, sep = "_", into = c("Canop")) %>%
  left_join(sir_canop_pos_Kai) %>%
  select(Canop, name, name)


selected_samples <- fread("20240125_selectedSamples_empo4type.txt")
short_name <- selected_samples %>%
  mutate(shortsample = paste(empo_4_type, row_number(), sep = "_")) %>%
  select(sample, shortsample)



#-----------Reading in intensity data GNSP-----------

gnpsid <- "MSV000083475"
setwd(paste(admin, "Pairedomics_projects",gnpsid, sep = "/"))

features <- fread("1907_EMPv2_INN_GNPS_quant.csv") %>%
  unique() %>%
  clean_names() %>%
  dplyr::filter(row_m_z <= 900) %>%
  select(row_id, selected_samples$sample) %>% #contains("blank")
  column_to_rownames("row_id")

feature_intensities <- as.data.frame(t(features))

feature_intensities_log <- feature_intensities %>%
  mutate_all(~ ifelse(. != 0, log10(.), 0))

constant_cols <- apply(feature_intensities_log, 2, function(col) length(unique(col)) == 1)
feature_intensities_log <- feature_intensities_log[, !constant_cols]

zero_cols <- apply(feature_intensities_log, 2, function(col) all(col == 0))
feature_intensities_log <- feature_intensities_log[, !zero_cols]

feature_intensities[feature_intensities > 0] <- 1

constant_cols <- apply(feature_intensities, 2, function(col) length(unique(col)) == 1)
feature_intensities <- feature_intensities[, !constant_cols]

zero_cols <- apply(feature_intensities, 2, function(col) all(col == 0))
feature_intensities <- feature_intensities[, !zero_cols]


#intensity one from Kai files

path_calc_fp <- paste(admin, "Datasets_fp/GNPS_Kai/CalculatedFingerprints", sep = "/")
setwd(path_calc_fp)

filename <- "20240206_sample_feature_intlog_fromKAIfiles"

feature_intensities_log <- sample_type %>%
  left_join(fread(paste(filename, "csv", sep = ".")), by = "sample") %>%
  select(sample, contains("PLATE")) %>% 
  filter(sample %in% selected_samples$sample) %>%
  column_to_rownames("sample")

constant_cols <- apply(feature_intensities_log, 2, function(col) length(unique(col)) == 1)
feature_intensities_log <- feature_intensities_log[, !constant_cols]

zero_cols <- apply(feature_intensities_log, 2, function(col) all(col == 0))
feature_intensities_log <- feature_intensities_log[, !zero_cols]

feature_intensities <- feature_intensities_log
feature_intensities[feature_intensities > 0] <- 1



#----------Reading in compiled fingerprints/classes from Kai-------------------

setwd(paste(admin,"Datasets_fp/GNPS_Kai", sep = "/"))

peaksKai <- fread("20230803_EMP_sirius_samples_Kai.csv") %>%
  bind_rows(fread("20230803_EMP_sirius_samples_Kai_extra.csv")) %>%
  select(samples, idname) %>%
  unique() %>%
  group_by(samples) %>%
  summarise(count = n()) %>%
  ungroup()

peakcount <- selected_samples %>%
  left_join(sample_type %>% select(sample, samples)) %>%
  left_join(peaksKai)

var_plot <- peakcount %>%
  ggplot(aes(x = empo_4_type, 
             y = count)) +
  geom_violin(aes(fill = empo_4_type), 
              scale = "width", width = 0.7) +
  #facet_wrap(~ Canop + name + FoldMean, scales = "free_y") +
  #theme(legend.position="none") +
  ggtitle("Amount of LC-HRMS features") + 
  scale_fill_manual(values = color_vector_empo_womix_new) + 
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  #theme_light()
  my_theme
var_plot



#---------Reading in SIRIUS data for clustering, comparing etc-------------

path_calc_fp <- paste(admin, "Datasets_fp/GNPS_Kai/CalculatedFingerprints", sep = "/")
setwd(path_calc_fp)

filenames <- dir(path_calc_fp, pattern = "sircanop")
filenames
filename <- "20231020_sircanopus_all_binartresh05_msgroupaver7"

all_aver_classif <- sample_type %>%
  left_join(names, by = "sample") %>%
  left_join(fread(paste(filename, "csv", sep = ".")), by = "filename")

all_aver_classif <- all_aver_classif %>%
  select(sample, grep("Canop", colnames(all_aver_classif))) %>%
  dplyr::filter(sample %in% selected_samples$sample) %>%
  column_to_rownames("sample") %>%
  select(selected_features_vector)

#eliminating zero rows and column and continuing with clustering
constant_cols <- apply(all_aver_classif, 2, function(col) length(unique(col)) == 1)
all_aver_classif <- all_aver_classif[, !constant_cols]

zero_cols <- apply(all_aver_classif, 2, function(col) all(col == 0))
all_aver_classif <- all_aver_classif[, !zero_cols]


forprint <- all_aver_classif %>%
  #select(selected_features_vector) %>%
  rownames_to_column("sample") %>%
  left_join(sample_type %>% select(sample, empo_4_type)) %>%
  select(sample, empo_4_type, everything())

setwd(paste(admin, "Manuscripts/EMP_article/DataTables", sep = "/"))
write_delim(forprint, "CANOPUS_massGr7_selectedlong_variables.csv")

  
print("ok")

#-----------Reading in ready final tables------------------
method_folder <- paste(admin, "Manuscripts/EMP_article/DataTables", sep = "/")
setwd(method_folder)
methodsdata <- dir(method_folder, pattern = ".csv")
methodsdata

all_aver_classif <- fread(methodsdata[13]) %>%
  #left_join(fread(methodsdata[12])) %>%
  #row_to_names(1) %>%
  select(-empo_4_type) %>%
  #left_join(short_name) %>%
  #select(-sample)  %>%
  column_to_rownames("sample")

#------------Clustering with PAM (Aris), data written out in folders---------------

modelname <- "FP mass Gr7 selectedlong variables"

# feature_intensities_log, #all_aver_classif
sample_data_comparison <- all_aver_classif

biomes <- sample_data_comparison %>%
  rownames_to_column("sample") %>%
  left_join(sample_type) %>%
  select(sample, empo_4_type)

biomes_label <- biomes %>%
  select(empo_4_type) %>%
  unique() %>%
  arrange(empo_4_type) %>%
  mutate(lable = row_number())

biomes <- biomes %>%
  left_join(biomes_label)


### using pam function for statistical analysis

pam1 <- pam(x = sample_data_comparison,
            k = as.numeric(length(unique(biomes$empo_4_type))),
            diss = FALSE)

cluster_table <- data.frame(cluster = pam1$clustering) %>%
  bind_cols(biomes)

st = cluster.stats(pam1$clustering,biomes$lable)
print(st)

cluster_plot <- fviz_cluster(pam1, data = sample_data_comparison)
#print(cluster_plot)

cluster_plot_xy <- as.data.frame(cluster_plot$data) %>%
  rownames_to_column("sample") %>%
  left_join(sample_type %>% select(sample, empo_4_type))

ggplot(cluster_plot_xy) +
  geom_point(aes(x=x, 
                 y=y,
                 colour = empo_4_type), #empo_4_type
             size = 5, alpha = 0.7) + 
  ggtitle(paste("PAM cluster, ", modelname)) + 
  scale_color_manual(values=color_vector_empo_womix) + 
  my_theme

#----------------PCA, UMAP, t-SNE-------------------------

# feature_intensities_log, #all_aver_classif
modelname <- "CANOPUS and FP selectedlong combined"

data_for_clustering <- all_aver_classif

set.seed(123)
data_pca <- prcomp(data_for_clustering, rank = 3)

clustering_p <- as.data.frame(data_pca$x) %>%
  rownames_to_column("sample")

PC1var <- round(summary(data_pca)$importance[2,1]*100,1)
PC2var <- round(summary(data_pca)$importance[2,2]*100,1)

#UMAP
set.seed(123)
umap <- umap(data_for_clustering)  #n_components = 3
clustering_u <- as.data.frame(umap$layout) %>%
  rownames_to_column("sample")
colnames(clustering_u) <- c("sample","UMAP1", "UMAP2")

#tSNE

set.seed(123)
data_pca100 <- prcomp(data_for_clustering, rank = 100)
#plot(x=1:100,y=summary(data_pca100)$importance[2,1:100])
pcatsne <- as.data.frame(data_pca100$x)
tsne <- Rtsne(pcatsne, dims = 2)

clustering_t <- data.frame(sample = rownames(data_for_clustering), tSNE1 = tsne$Y[, 1], tSNE2 = tsne$Y[, 2])


###graph
graph_p <- clustering_p %>%
  left_join(sample_type, by = "sample")

graph_u <- clustering_u %>%
  left_join(sample_type, by = "sample")

graph_t <- clustering_t %>%
  left_join(sample_type, by = "sample")


ggplot(data = graph_p) +
  geom_point(mapping = aes(x = PC1,
                           colour = empo_4_type,# empo_4_type, #empo_4_type
                           #shape = empo_4_salinity,
                           y =PC2),
             size = 5,
             alpha = 0.7) +
  my_theme +
  #theme_light() +
  scale_color_manual(values=color_vector_empo_womix_new) + 
  ggtitle(modelname) +
  #facet_wrap(~empo_4_type, ncol = 4) +
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1)) +
  xlab(paste("PC1, ", PC1var, "%", sep = "")) +
  ylab(paste("PC2, ", PC2var, "%", sep = ""))


ggplot(data = graph_u) +
  geom_point(mapping = aes(x = UMAP1,
                           colour = empo_4_type, 
                           #shape = empo_4_salinity,
                           y =UMAP2),
             size = 5,
             alpha = 0.7) +
  my_theme +
  #theme_light() +
  scale_color_manual(values=color_vector_empo_womix_new) + 
  ggtitle(modelname) +
  #facet_wrap(~empo_4_type, ncol = 4) +
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1)) +
  xlab("UMAP1") +
  ylab("UMAP2")



ggplot(data = graph_t) +
  geom_point(mapping = aes(x = tSNE1,
                           colour = empo_4_type, 
                           #shape = empo_4_salinity,
                           y =tSNE2),
             size = 5,
             alpha = 0.7) +
  my_theme +
  #theme_light() +
  scale_color_manual(values=color_vector_empo_womix_new) + 
  ggtitle(modelname) +
  #facet_wrap(~empo_4_type, ncol = 4) +
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1)) +
  xlab("tSNE1") +
  ylab("tSNE2")



#----------------Distance between samples----------------

## feature_intensities_log, #all_aver_classif
sample_for_distance <- all_aver_classif %>%
  rownames_to_column("sample") %>%
  left_join(short_name) %>%
  select(-sample)  %>%
  column_to_rownames("shortsample")

#Euclideal distance
dist_matrix <- dist(sample_for_distance, method = "euclidean")
similarity_matrix <- as.matrix(dist_matrix)
similarity_matrix[similarity_matrix == 0] <- NA
max = max(similarity_matrix, na.rm = TRUE)
min = min(similarity_matrix, na.rm = TRUE)

dev.off()
par(mar = c(2, 2, 2, 2))
ggplot(data = reshape2::melt(similarity_matrix), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "royalblue4", high = "orange", limits = c(0, 8)) + #limits = c(0, 2)
  #labs(x = "Columns", y = "Rows", fill = "Values") +
  ggtitle(modelname) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_c(begin = 1, end = 0, limits = c((min-0.01),(max+0.01))) #limits = c(0,1)

#print("ok")


similarity_matrix_aver <- as.data.frame(similarity_matrix) %>%
  rownames_to_column("sample1") %>%
  gather(key = "sample2", value = "distance", -sample1) %>%
  na.omit() %>%
  separate(sample1, into = c("type1", "sample1"), sep = "_") %>%
  separate(sample2, into = c("type2", "sample2"), sep = "_") %>%
  group_by(type1, type2) %>%
  mutate(dist_aver = mean(distance)) %>%
  ungroup() %>%
  select(type1, type2, dist_aver) %>%
  unique() %>%
  #mutate(dist_aver = rescale(dist_aver, to = c(1, 2))) %>%
  spread(key = type2, value = dist_aver) %>%
  column_to_rownames("type1")

max = max(similarity_matrix_aver, na.rm = TRUE)
min = min(similarity_matrix_aver, na.rm = TRUE)

similarity_matrix_avertype <- as.matrix(similarity_matrix_aver)

max = max(similarity_matrix_avertype, na.rm = TRUE)
min = min(similarity_matrix_avertype, na.rm = TRUE)

dev.off()
par(mar = c(2, 2, 2, 2))
ggplot(data = reshape2::melt(similarity_matrix_avertype), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "royalblue4", high = "orange", limits = c(0, 8)) + #limits = c(0, 2)
  #labs(x = "Columns", y = "Rows", fill = "Values") +
  ggtitle(modelname) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_c(begin = 1, end = 0, limits = c((min-0.01),(max+0.01))) #limits = c(0,1)


#scaled between 1 and 2

#-------------Obtaining statistic values for comparing methods-----------------


## feature_intensities_log, #all_aver_classif
sample_for_distance <- all_aver_classif %>%
  rownames_to_column("sample") %>%
  left_join(short_name) %>%
  select(-sample)  %>%
  column_to_rownames("shortsample")

#Euclideal distance
dist_matrix <- dist(sample_for_distance, method = "euclidean")
similarity_matrix <- as.matrix(dist_matrix)
#similarity_matrix[similarity_matrix == 0] <- NA
max = max(similarity_matrix, na.rm = TRUE)
min = min(similarity_matrix, na.rm = TRUE)

#scaling values

similarity_matrix_scaled <- as.data.frame(similarity_matrix) %>%
  rownames_to_column("sample1") %>%
  gather(key = "sample2", value = "dist_aver", -sample1) %>%
  #mutate(dist_aver = dist_aver/max) %>%
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

# Print the metrics
print(metrics)

methodname <- "CANOPUSandFPmassGr7SelLong" 
aver_within <- metrics$average.within
aver_between <- metrics$average.between
siluettwidth <- metrics$avg.silwidth
pearsongamma <- metrics$pearsongamma
dunn <- metrics$dunn
ch <- metrics$ch  #bigger better

#statistics <- data.frame(matrix(ncol = 7, nrow = 1))
#colnames(statistics) <- c("methodname", "aver_within", "aver_between",
#                          "siluettewidth","pearsongamma", "dunn", "ch")

datarow <- c(methodname, aver_within, aver_between,siluettwidth, pearsongamma, dunn, ch)
datarow
statistics <- rbind(statistics, datarow)

setwd(paste(admin, "Datasets_fp/GNPS_Kai/Statistics", sep = "/"))
#write_delim(statistics, "20240129_clusterstats_scaleddata.txt", delim = ";")

print("ok")

separationmatrix <- metrics$separation.matrix
max = max(separationmatrix, na.rm = TRUE)
min = min(separationmatrix, na.rm = TRUE)

dev.off()
par(mar = c(2, 2, 2, 2))
ggplot(data = reshape2::melt(separationmatrix), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "royalblue4", high = "orange", limits = c(0, 8)) + #limits = c(0, 2)
  #labs(x = "Columns", y = "Rows", fill = "Values") +
  ggtitle(modelname) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_c(begin = 1, end = 0, limits = c(0,(max+0.01))) #limits = c(0,1)



#-----------comparing methods with already calculated values----------------------

setwd(paste(admin, "Datasets_fp/GNPS_Kai/Statistics", sep = "/"))
statistics <- fread("20240129_clusterstats_scaleddata.txt") %>%
  gather(key = "statmethod", value = "value", -c(methodname)) %>%
  #dplyr::filter(grepl("aver", statmethod)) %>%
  dplyr::filter(statmethod %in% c("CHindex", "DunnIndex", "SilhouetteIndex"))


plotbar <- statistics %>%
  dplyr::filter(!grepl("MS2", methodname)) %>%
  dplyr::filter(!grepl("and", methodname)) %>%
  dplyr::filter(methodname %in% c("MS1", "MS1SelLong","MS1logInt", "MS1logSelLong",
                                  "FPaver", "FPaverSelLong", "FPaverSel",
                                  "FPmassGr7","FPmassGr7SelLong","FPmassGr7Sel",
                                  "CANOPUSaverAll", "CANOPUSaverSelLong", "CANOPUSaverSel",
                                  "CANOPUSmassGr7","CANOPUSmassGr7SelLong", "CANOPUSmassGr7Sel")) %>%
  mutate(methodname = fct_relevel(methodname,
                                  "MS1", "MS1SelLong","MS1logInt", "MS1logSelLong",
                                  "FPaver", "FPaverSelLong", "FPaverSel",
                                  "FPmassGr7","FPmassGr7SelLong","FPmassGr7Sel",
                                  "CANOPUSaverAll", "CANOPUSaverSelLong", "CANOPUSaverSel",
                                  "CANOPUSmassGr7","CANOPUSmassGr7SelLong", "CANOPUSmassGr7Sel")) %>%
  ggplot(aes(x = factor(methodname),
             y = value,
             fill = methodname)) +
  geom_bar(stat = "identity", position = 'dodge') +
  labs(title = "Comparison of methods",
       x = "Method",
       y = "Metric value") +
  #theme(legend.position="none") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values=c("#C62828","#E57373","#F57C00","#FFA726",
                             "#311B92", "#7E57C2","#D1C4E9","#01579B","#03A9F4", "#B3E5FC",
                             "#1B5E20","#66BB6A", "#C8E6C9","#9E9D24","#D4E157","#F0F4C3")) +
  facet_wrap(~statmethod, scales = "free_y", nrow = 1) +
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1))
plotbar  

setwd(paste(admin, "Manuscripts/EMP_article/Graphs", sep = "/"))
ggsave("legend2.png", width = 40, height = 20, units = "cm")

#------------Analysing more important variables--------------

setwd(paste(admin, "Manuscripts/EMP_article/DataTables", sep = "/"))

all_aver_classif <- fread("CANOPUS_average_selected_variables.csv")

#clustering similar features together
for_clustering_variables <- as.data.frame(t(all_aver_classif))

for_clustering_variables_m <- as.matrix(for_clustering_variables)
k <- 6  # Number of clusters
kmeans_result <- kmeans(for_clustering_variables_m, centers = k)
cluster_assignment <- as.data.frame(kmeans_result$cluster) %>%
  rownames_to_column("var")
colnames(cluster_assignment)[2] <- "Cluster"

#hierarchical clustering

for_clustering_variables <- as.data.frame(t(all_aver_classif %>%
                                              select(one_of(selected_features_vector))))
dist_matrix <- dist(for_clustering_variables)
hclust_result <- hclust(dist_matrix, method = "complete")
sample_order <- as.vector(cutree(hclust_result, k = nrow(for_clustering_variables)))
cluster_assignment <- cbind(for_clustering_variables, Cluster = sample_order) %>%
  rownames_to_column("var") %>%
  select(var, Cluster)


#graphs for variables
variables_graph <- all_aver_classif %>%
  #select(one_of(selected_features_vector)) %>%
  #rownames_to_column("sample") %>%
  #left_join(biomes) %>%
  #select(sample, empo_4_type, everything(), -lable) %>%
  gather(key = "Canop", value = "average", -c(empo_4_type, sample)) %>%
  left_join(selected_features %>%
              select(var, FoldMean) %>%
              mutate(Canop = var)) %>%
  separate(Canop, into = c("Canop"), sep = "_") %>%
  left_join(selected_features_described) %>%
  #left_join(cluster_assignment) %>%
  group_by(Canop) %>%
  mutate(Canop_aver = mean(average)) %>%
  ungroup() %>%
  group_by(Canop, empo_4_type) %>%
  mutate(Canop_aver_biome = mean(average)) %>%
  ungroup() 

Canop_var_decreasing <- variables_graph %>%
  select(Canop, Canop_aver) %>%
  unique() %>%
  arrange(desc(Canop_aver)) %>%
  mutate(CanVarDesc = row_number())

variables_graph <- variables_graph %>%
  left_join(Canop_var_decreasing) 



#violin
var_plot <- variables_graph %>%
  mutate(FoldMean = round(FoldMean, digits = 2)) %>%
  #dplyr::filter(CanVarDesc < 20) %>%
  arrange(CanVarDesc) %>%
  #dplyr::filter(FoldMean > 2.78) %>%
  #dplyr::filter(Cluster %in% c(1:7)) %>%
  #dplyr::filter(Canop == "Canop516") %>%
  ggplot(aes(x = CanVarDesc, 
             y = average)) +
  geom_violin(aes(fill = Canop), 
              scale = "width", width = 0.7) +
  facet_wrap(~ empo_4_type, nrow = 2) + #scales = "free_y" 
  theme(legend.position="none") +
  xlab("CANOPUS classes reordered based on their overall average value") +
  #ggtitle("Cluster one") + 
  #scale_fill_manual(values = color_vector_empo_womix_new) + 
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  #theme_light()
  my_theme
var_plot


var_plot <- variables_graph %>%
  mutate(FoldMean = round(FoldMean, digits = 2)) %>%
  arrange(FoldMean) %>%
  #dplyr::filter(FoldMean > 2.78) %>%
  #dplyr::filter(Cluster %in% c(1:7)) %>%
  #dplyr::filter(Canop == "Canop516") %>%
  ggplot(aes(x = empo_4_type, 
             y = average)) +
  geom_violin(aes(fill = empo_4_type), 
              scale = "width", width = 0.7) +
  facet_wrap(~ Canop + name + FoldMean, scales = "free_y") +
  theme(legend.position="none") +
  #ggtitle("Cluster one") + 
  scale_fill_manual(values = color_vector_empo_womix_new) + 
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  #theme_light()
  my_theme
var_plot



#graph class, average in class, but classes ranged from low to high

var_graph_line <- variables_graph %>%
  group_by(Canop, empo_4_type) %>%
  mutate(biomeAver = mean(average),
         biomeMedian = median(average),
         biomeMin = min(average[average != 0], na.rm = TRUE),
         biomeMax = max(average[average != 0], na.rm = TRUE)) %>%
  ungroup() %>%
  select(Canop, empo_4_type, biomeAver, biomeMedian, biomeMin, biomeMax, CanVarDesc) %>%
  unique() %>%
  mutate_all(~ ifelse(is.infinite(.), 0, .))



ggplot(data = var_graph_line) +
  geom_line(aes(x = CanVarDesc, y = biomeMedian, colour = empo_4_type),
             size = 1) +
  # geom_text(aes(x = methodnr,
  #               label = methodname,
  #               y = value), size = 3, angle = 30) +
  my_theme +
  #theme_light() +
  scale_color_manual(values=color_vector_empo_womix) + 
  #ggtitle(modelname) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  #facet_wrap(~statmethod, scales = "free_y", nrow = 1) +
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1))
  

#-----------------------END------------------------

