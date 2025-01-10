###Code by Pilleriin Peets (pilleriin.peets@gmail.com)

##Figure data and code for article:
##"Chemistry-based vectors map the chemical space of natural biomes from untargeted mass spectrometry data"
##by Peets et al.
#Data for these Figures is in Zenodo

#----------Packages---------
library(tidyverse)
library(janitor)
library(data.table)
library(umap)
library(plotly)
library(viridis)
library(broom)


#----------Files reading in-------
setwd("C:/Users/pille/OneDrive - Friedrich-Schiller-Universität Jena/JenaUniversity/Manuscripts/EMP_article/FinalZenodoData")

emp_sample_info <- fread("EMP_sample_data.csv")

color_vector_empo_womix_new <- c("#9467bd","#d62728","#e377c2","#ff7f0e",
                                 "#FFCC33","#bcbd22","#2ca02c","#CC9966",
                                 "#660000","#999999","#3300CC")

#----------Figure2----------

#Figure 2a
sample_peaks_int <- fread("MS1_log_all_variables.csv")

peakcount <- gather(sample_peaks_int, key = "feature", value = "area", -c(sample, empo_4_type)) %>%
  dplyr::filter(area > 0) %>%
  select(sample, empo_4_type, feature) %>%
  unique() %>%
  group_by(sample, empo_4_type) %>%
  summarise(peakcount = n()) %>%
  ungroup() 


# Figure 2b

peakcount <- fread("emp_sample_peak_connector.csv") %>%
  inner_join(fread("EMP_MFP_all.csv")) %>%
  select(samples, idname, empo_4_type) %>%
  group_by(samples, empo_4_type) %>%
  summarise(peakcount = n()) %>%
  ungroup()


var_plot <- peakcount %>%
  ggplot(aes(x = empo_4_type, 
             y = peakcount)) +
  geom_violin(aes(fill = empo_4_type), 
              scale = "width", width = 0.7) +
  theme(legend.position="none") +
  scale_fill_manual(values = color_vector_empo_womix_new) #+ 
  #my_theme
var_plot


#Figure 2c

pca <- fread("Compounds_PCA_MFP.csv") %>%
  mutate(RT = round(rt, digits = 0)) %>%
  group_by(RT) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  arrange(desc(count))

plot <- ggplot(data = pca) +
  geom_point(mapping = aes(x = PC1,
                           colour = rt,
                           #size = mz,
                           y =PC2),
             size = 1,
             alpha = 0.7) +
  #my_theme +
  facet_wrap(~Containing_nitrogen, nrow = 1) +
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1)) +
  scale_color_viridis(option = "H") +
  #scale_color_manual(values= c("#000096","#FF9933", "lightgrey")) +
  xlab(paste("PC1, ", "13.3", "%", sep = "")) +
  ylab(paste("PC2, ", "9.9", "%", sep = ""))
plot

#Figure 2d

pca <- fread("Compounds_PCA_CC.csv") %>%
  mutate(RT = round(rt, digits = 0)) %>%
  group_by(RT) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  arrange(desc(count))

plot <- ggplot(data = pca) +
  geom_point(mapping = aes(x = PC1,
                           colour = rt,
                           #size = mz,
                           y =PC2),
             size = 1,
             alpha = 0.7) +
  #my_theme +
  #facet_wrap(~Containing_nitrogen, nrow = 1) +
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1)) +
  scale_color_viridis(option = "H") +
  #scale_color_manual(values= c("#000096","#FF9933", "lightgrey")) +
  xlab(paste("PC1, ", "13.3", "%", sep = "")) +
  ylab(paste("PC2, ", "9.9", "%", sep = ""))
plot


#Figure 2e

plot <- ggplot(pca, aes(x = rt)) +
  #geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  geom_density(color = "black", adjust = 1, linewidth = 2) +
  #my_theme +
  labs(title = "Distribution of RT", x = "Retention time (min)", y = "Peak density")
plot


#Figure 2h
pca$Containing_nitrogen[is.na(pca$Containing_nitrogen)] <- 'Unknown'
pca$Containing_nitrogen <- factor(pca$Containing_nitrogen, 
                                      levels = c('yes', 'no', 'Unknown'))
graph_p_proportions <- as.data.frame(prop.table(table(pca$Containing_nitrogen)))

yesnocount <- pca %>%
  group_by(Containing_nitrogen) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  select(Containing_nitrogen, count) %>% unique()

plot <- ggplot(graph_p_proportions, aes(x = Var1, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Containing Nitrogen (Including NA as 'Unknown')", 
       x = "Containing Nitrogen", 
       y = "Proportion") +
  scale_fill_manual(values = c('yes' = "#FF9933", 'no' = "#000096", 'Unknown' = 'lightgray'))  +
  #my_theme +
  ylim(0,0.65) 
plot


#----------Figure3----------

#Figure 3a-c, also all others in SI. Just choose the dataset.

data <- fread("SI6_MFP_mass_rfe1_variables.csv") %>%
  select(-empo_4_type) %>%
  column_to_rownames("sample")

set.seed(123)
data_pca <- prcomp(data, rank = 2)

pca <- as.data.frame(data_pca$x) %>%
  rownames_to_column("sample")

PC1var <- round(summary(data_pca)$importance[2,1]*100,1)
PC2var <- round(summary(data_pca)$importance[2,2]*100,1)

graph_p <- pca %>%
  left_join(fread("EMP_sample_data.csv"), by = "sample") %>%
  select(sample, PC1, PC2, empo_4_type)

write_delim(graph_p, "print.csv", delim = ";")

plot <- ggplot(data = graph_p) +
  geom_point(mapping = aes(x = PC1,
                           colour = empo_4_type,
                           y =PC2),
             size = 3,
             alpha = 0.7) +
  #my_theme +
  theme_light() +
  scale_color_manual(values=color_vector_empo_womix_new) +
  #theme(legend.position="none") +
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1)) +
  xlab(paste("PC1, ", PC1var, "%", sep = "")) +
  ylab(paste("PC2, ", PC2var, "%", sep = ""))
plot


#Figure 3d-f, also all others in SI. Just choose the dataset.

data <- fread("SI3_MFP_average_rfe1_variables.csv") %>%
  select(-empo_4_type) %>%
  column_to_rownames("sample")

set.seed(123)
umap <- umap(data)  #n_components = 3
umap <- as.data.frame(umap$layout) %>%
  rownames_to_column("sample")
colnames(umap) <- c("sample","UMAP1", "UMAP2")

graph_u <- umap %>%
  left_join(fread("EMP_sample_data.csv"), by = "sample") %>%
  select(sample, UMAP1, UMAP2, empo_4_type)


plot <- ggplot(data = graph_u) +
  geom_point(mapping = aes(x = UMAP1,
                           colour = empo_4_type,
                           y =UMAP2),
             size = 3,
             alpha = 0.7) +
  #my_theme +
  theme_light() +
  scale_color_manual(values=color_vector_empo_womix_new) +
  #theme(legend.position="none") +
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1))
plot


#----------Figure4----------

#UMAP from 3f, mass grouped MFP rfe1

data <- fread("SI6_MFP_mass_rfe1_variables.csv") %>%
  select(-empo_4_type) %>%
  column_to_rownames("sample")

set.seed(123)
umap <- umap(data)  #n_components = 3
umap <- as.data.frame(umap$layout) %>%
  rownames_to_column("sample")
colnames(umap) <- c("sample","UMAP1", "UMAP2")

graph_u <- umap %>%
  left_join(fread("EMP_sample_data.csv"), by = "sample") %>%
  select(sample, UMAP1, UMAP2, empo_4_type, emp500_principal_investigator, emp500_title)

forprint <- graph_u %>% 
  dplyr::filter(empo_4_type == "Animal distal gut")

write_delim(forprint, "print.csv", delim = ";")

graph_u_mod <- graph_u %>%
  mutate(Origin = paste(emp500_title, " (", emp500_principal_investigator, ")", sep = ""))

plot_u <- ggplot(data = graph_u_mod) +
  geom_point(mapping = aes(x=UMAP1, 
                           y= UMAP2), 
             color = "lightgrey", 
             size = 3,
             alpha = 0.1) +
  geom_point(data = graph_u_mod %>% dplyr::filter(empo_4_type == "Animal distal gut"),  #choose biome to check
             mapping = aes(x = UMAP1,
                           colour = Origin, #emp500_title,#emp500_principal_investigator,#empo_4_type, 
                           y =UMAP2),
             size = 3,
             alpha = 0.5) +
  my_theme +
  #theme_light() +
  scale_color_manual(values = c("#33A02C","#E31A1C","#1F78B4","#6A3D9A","#FF7F00", "#66C2A5",
                                "#FFD92F","#E78AC3",
                                "#CAB2D6","#FDBF6F","#9467bd",
                                "#8DA0CB","#A6D854","#E5C494", "#B15928",
                                "#B2DF8A","#FF7F00", "black", "#FB9A99")) +
  #theme(legend.position="none") +
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1)) +
  xlab("UMAP1") +
  ylab("UMAP2")

plot_u



#----------Figure5----------

#CCV values for compound classes from rfe2 filtered dataset

data <- fread("SI10_CC_average_rfe2_variables.csv")

compclasses <- fread("canopus.tsv") %>%
  mutate(Canop = paste("Canop", absoluteIndex, sep = ""))

variables_graph <- data %>%
  gather(key = "Canop", value = "average", -c(empo_4_type, sample)) %>%
  separate(Canop, into = c("Canop"), sep = "_") %>%
  left_join(compclasses %>% select(Canop, name)) %>%
  group_by(Canop, empo_4_type) %>%
  mutate(Canop_aver_biome = mean(average)) %>%
  ungroup() 

variables_graph <- variables_graph %>%
  dplyr::filter(Canop %in% c("Canop4138"))
  dplyr::filter(Canop %in% c("Canop254", "Canop1540", "Canop3439", "Canop259", "Canop258", "Canop1444", "Canop4138"))

write_delim(variables_graph, "print.csv", delim = ";")

#graph
var_plot <- ggplot(data = variables_graph, 
                   aes(x = empo_4_type,
                       y = average)) +
  geom_violin(aes(fill = empo_4_type),
              scale = "width", width = 0.7) +
  
  geom_point(aes(y = Canop_aver_biome),
             color = "black", shape = 21, size = 1, fill = "black") +
  facet_wrap(~ Canop + name, nrow = 2, scales = "free_y") + #CanVarDesc +  # ~ var + name , scales = "free_y"
  theme(legend.position="none") +
  scale_fill_manual(values = color_vector_empo_womix_new) +
  theme_light()
var_plot


#----------Figure6----------

#Fig 6a-d

compclasses <- fread("canopus.tsv") %>%
  mutate(Canop = paste("Canop", absoluteIndex, sep = ""))


selected_cc <- fread("distalgut_groupanalysis.csv") %>%
  dplyr::filter(Canop %in% c("Canop2", "Canop321", "Canop278", "Canop331"))

p_value_results <- selected_cc %>%
  group_by(Canop) %>%
  do(tidy(wilcox.test(average ~ cluster, data = ., exact = FALSE))) %>%
  ungroup()



write_delim(selected_cc, "print.csv", delim = ";")


#Fig 6e-h

selected_cc <- fread("animalsecretion_groupanalysis.csv") %>%
  mutate(cluster = as.character(cluster)) %>%
  dplyr::filter(Canop %in% c("Canop128", "Canop4734", "Canop1830", "Canop12"))

p_value_results <- selected_cc %>%
  group_by(Canop) %>%
  do(tidy(wilcox.test(average ~ cluster, data = ., exact = FALSE))) %>%
  ungroup()

write_delim(selected_cc, "print.csv", delim = ";")



#----------SI Figure S4---------

folder <- "C:/Users/pille/OneDrive - Friedrich-Schiller-Universität Jena/JenaUniversity/Manuscripts/EMP_article/FinalZenodoData"
setwd(folder)
dir(folder, pattern = ".csv")

compclasses <- fread("canopus.tsv") %>%
  mutate(Canop = paste("Canop", absoluteIndex, sep = ""))

data <- fread("SI10_CC_average_rfe2_variables.csv") 

rfe <- fread("20240625_rfe_canop_120var.csv") %>%
  na.omit() %>%
  mutate(var = as.character(var)) %>%
  mutate(FoldMean = rowMeans(select(., contains("Fold")), na.rm = TRUE))
selected_features_vector <- rfe$var


#clustering similar features together
for_clustering_variables <- as.data.frame(t(data %>% 
                                              select(-empo_4_type) %>% 
                                              column_to_rownames("sample")))

for_clustering_variables_m <- as.matrix(for_clustering_variables)
k <- 10  # Number of clusters
kmeans_result <- kmeans(for_clustering_variables_m, centers = k)
cluster_assignment <- as.data.frame(kmeans_result$cluster) %>%
  rownames_to_column("var")
colnames(cluster_assignment)[2] <- "Cluster"


#hierarchical clustering
for_clustering_variables <- as.data.frame(t(for_clustering_variables))

dist_matrix <- dist(for_clustering_variables)
hclust_result <- hclust(dist_matrix, method = "complete")
sample_order <- as.vector(cutree(hclust_result, k = nrow(for_clustering_variables)))
cluster_assignment <- cbind(for_clustering_variables, Cluster = sample_order) %>%
  rownames_to_column("var") %>%
  select(var, Cluster)


#graphs for variables
variables_graph <- data %>%
  gather(key = "Canop", value = "average", -c(empo_4_type, sample)) %>%
  left_join(rfe %>%
              select(var, FoldMean) %>%  #FoldMean
              mutate(Canop = var)) %>%
  separate(Canop, into = c("Canop"), sep = "_") %>%
  left_join(compclasses %>%
              select(Canop, name)) %>%
  group_by(var) %>%
  mutate(Canop_aver = mean(average)) %>%
  ungroup() %>%
  group_by(var, empo_4_type) %>%
  mutate(Canop_aver_biome = mean(average)) %>%
  ungroup() 

Canop_var_decreasing <- variables_graph %>%
  select(var, Canop_aver) %>%
  unique() %>%
  arrange(desc(Canop_aver)) %>%
  mutate(CanVarDesc = row_number())


#data for final violin graph
variables_graph_ <- variables_graph %>%
  left_join(Canop_var_decreasing) %>%
  mutate(FoldMean = round(FoldMean, digits = 2)) %>%
  left_join(fread("cc_order_delete.csv"))
  #mutate(CanVarDesc = as.factor(CanVarDesc)) %>%
  #arrange(CanVarDesc) %>%
  #dplyr::filter(Canop_aver > 0.1)
  #dplyr::filter(Canop %in% c("Canop254", "Canop1540", "Canop3439", "Canop259", "Canop258", "Canop1444", "Canop4138"))

names <- variables_graph_ %>%
  select(row, Canop, name) %>%
  unique()


#graph
var_plot <- ggplot(data = variables_graph_, 
                   aes(x = empo_4_type,
                       y = average)) +
  geom_violin(aes(fill = empo_4_type),
              scale = "width", width = 0.7) +
  
  geom_point(aes(y = Canop_aver_biome),
             color = "black", shape = 21, size = 1, fill = "black") +
  
  # geom_crossbar(aes(y = Canop_aver_biome, ymin = Canop_aver_biome, ymax = Canop_aver_biome),
  #               width = 1, color = "black", size = 0.3) +
  
  
  facet_wrap(~ row + Canop + FoldMean, ncol = 4, scales = "free_y") + #CanVarDesc +  # ~ var + name , scales = "free_y"
  theme(legend.position="none") +
  #ggtitle("Cluster one") + 
  scale_fill_manual(values = color_vector_empo_womix_new) + 
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  #theme_light()
  my_theme
var_plot


setwd(paste(admin, "Manuscripts/EMP_article/Graphs/svg", sep = "/"))
svg("CC_S4_orig.svg", width = 15, height = 40, pointsize = 1)
print(var_plot)
dev.off()


#----------SI Figure S5------------

folder <- "C:/Users/pille/OneDrive - Friedrich-Schiller-Universität Jena/JenaUniversity/Manuscripts/EMP_article/FinalZenodoData"
setwd(folder)
dir(folder, pattern = ".csv")

compclasses <- fread("canopus.tsv") %>%
  mutate(Canop = paste("Canop", absoluteIndex, sep = ""))

data <- fread("SI11_CC_mass_all_variables.csv") %>%
  gather(key = "Canop", value = "average", -c(empo_4_type, sample)) %>%
  separate(Canop, into = c("Canop", "massgroup"), sep = "_") %>%
  left_join(compclasses %>%
              select(Canop, name)) %>%
  dplyr::filter(Canop %in% c("Canop254", "Canop1540", "Canop3439", "Canop259", "Canop258", "Canop1444", "Canop4138")) %>%
  group_by(Canop, empo_4_type, massgroup) %>%
  mutate(Canop_aver_biome = mean(average)) %>%
  ungroup()

#graph
var_plot <- ggplot(data = data, 
                   aes(x = empo_4_type,
                       y = average)) +
  geom_violin(aes(fill = empo_4_type),
              scale = "width", width = 0.7) +
  
  geom_point(aes(y = Canop_aver_biome),
             color = "black", shape = 21, size = 1, fill = "black") +
  
  # geom_crossbar(aes(y = Canop_aver_biome, ymin = Canop_aver_biome, ymax = Canop_aver_biome),
  #               width = 1, color = "black", size = 0.3) +
  
  
  facet_wrap(~ name + massgroup, ncol = 7, scales = "free_y") + #CanVarDesc +  # ~ var + name , scales = "free_y"
  theme(legend.position="none") +
  #ggtitle("Cluster one") + 
  my_theme +
  scale_fill_manual(values = color_vector_empo_womix_new)

var_plot

setwd(paste(admin, "Manuscripts/EMP_article/Graphs/svg", sep = "/"))
svg("CC_S5_orig.svg", width = 15, height = 20, pointsize = 1)
print(var_plot)
dev.off()


#-----------Figure S6,S8-----------------

#Figure S6
folder <- "C:/Users/pille/OneDrive - Friedrich-Schiller-Universität Jena/JenaUniversity/Manuscripts/EMP_article/FinalZenodoData"
setwd(folder)
dir(folder, pattern = ".csv")

data <- fread("distalgut_groupanalysis.csv") %>%
  mutate(cluster = as.character(cluster))

var_plot <- ggplot(data = data, 
                   aes(x = cluster,
                       y = average)) +
  geom_violin(aes(fill = cluster),
              scale = "width", width = 0.7) +
  facet_wrap(~ Canop + name, scales = "free_y", ncol = 5) + #CanVarDesc +  # ~ var + name , scales = "free_y"
  theme(legend.position="none") +
  #scale_fill_manual(values = color_vector_empo_womix_new) + 
  scale_fill_manual(values = c("#6A3D9A","#E31A1C")) + #"#ff7f0e","#1F78B4"
  #theme_light()
  my_theme
var_plot

nametable <- data %>%
  select(Canop, name) %>%
  unique()

setwd(paste(admin, "Manuscripts/EMP_article/Graphs/svg", sep = "/"))
svg("CC_S6_origsmall.svg", width = 15, height = 17, pointsize = 1)
print(var_plot)
dev.off()


#Figure S8
folder <- "C:/Users/pille/OneDrive - Friedrich-Schiller-Universität Jena/JenaUniversity/Manuscripts/EMP_article/FinalZenodoData"
setwd(folder)
dir(folder, pattern = ".csv")

data <- fread("animalsecretion_groupanalysis.csv") %>%
  mutate(cluster = as.character(cluster))

var_plot <- ggplot(data = data, 
                   aes(x = cluster,
                       y = average)) +
  geom_violin(aes(fill = cluster),
              scale = "width", width = 0.7) +
  facet_wrap(~ Canop + name, scales = "free_y", ncol = 5) + #CanVarDesc +  # ~ var + name , scales = "free_y"
  theme(legend.position="none") +
  #scale_fill_manual(values = color_vector_empo_womix_new) + 
  scale_fill_manual(values = c("#ff7f0e","#1F78B4")) + #"#ff7f0e","#1F78B4"
  #theme_light()
  my_theme
var_plot

nametable <- data %>%
  select(Canop, name) %>%
  unique()

setwd(paste(admin, "Manuscripts/EMP_article/Graphs/svg", sep = "/"))
svg("CC_S8_origbig.svg", width = 15, height = 30, pointsize = 1)
print(var_plot)
dev.off()
