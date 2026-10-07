#Author: Pilleriin Peets, pilleriin.peets@gmail.com, pilleriin.peets@ut.ee Start: 2023/08/25
###Analysing Meine and Milous data

library(tidyverse)
library(data.table)
library(ggplot2)
library(viridis)
library(janitor)
library(umap)
library(plotly)
#library(Rdisop)
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
library(rstatix)

#find("p.adjust")

admin <- "DATA"

#---*---------Functions-------------

Meine_theme <- theme(
  plot.background = element_blank(),
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray88", size = 0.5),
  panel.grid.minor = element_line(color = "gray95", size = 0.25),
  #axis.line = element_line(size = 0.5, color = basecolor),
  plot.title = element_text(color = basecolor,
                            size = fontsize,
                            face = "bold",
                            hjust = 0.5, vjust = 1),
  text = element_text(family = font,
                      size = fontsize,
                      color = basecolor),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(family = font,
                            size = fontsize,
                            color = basecolor),
  legend.text = element_text(family = font,
                             size = fontsize,
                             color = basecolor),
  axis.text = element_text(family = font,
                           size = fontsize-2,
                           color = basecolor),
  axis.ticks = element_blank(),
  aspect.ratio = 1,
  axis.title.x = element_text(hjust = 0.5, vjust = 1),
  axis.title.y = element_text(hjust = 0.5, vjust = 1)
)

colors <- c("#039c75", "#f1e345", "#58b2e9", "#e79e02", "#0271b1")


group_c <- c(
  Surface = "#009E73",Platform = "#ddcc77",Shallow_Sinkhole = "#a2e7f2",Deep_Sinkhole = "#0271b1",
  Acid_Lake = "#e79e02",Mesopelagic = "plum",Exclude = "#000000")




#-----------SIRIUS fp and classes-----------

setwd(paste(admin, "AcidSinkhole/SIRIUS", sep = "/"))
sir_fp_pos_Kai <- fread("csi_fingerid.tsv") %>%
  mutate(Un = paste("Un", absoluteIndex, sep = ""))
sir_canop_pos_Kai <- fread("canopus.tsv") %>%
  mutate(Canop = paste("Canop", absoluteIndex, sep = "")) %>%
  select(Canop, name, everything())


#---*---------metadata-------------

setwd(paste(admin, "AcidSinkhole", sep = "/"))

newgrouping <- fread("Sinkholedata_attempt6_CraigClusters.csv") %>%
  #select(FileName, CraigClusters, CraigClusters_noExclude)
  select(FileName, Meine_Clusters, Craig_Clusters)

sampledata <- fread("Sinkholedata_attempt3_sheet1.csv") %>%
  select(-V32) %>%
  dplyr::filter(grepl("O2006", FileName))

cluster <- sampledata %>%
  select(FileName, p_h) %>%
  na.omit() %>%
  mutate(p_h = as.numeric(gsub(",", ".", p_h))) %>%
  mutate(phcluster = "1")

for(n in 1:length(cluster$FileName)){
  if(cluster$p_h[n] < 7){
    cluster$phcluster[n] = 2
  }
}

sampledata <- sampledata %>%
  left_join(cluster %>% select(FileName, phcluster)) %>%
  left_join(newgrouping) %>%
  #mutate(cluster = CraigClusters_noExclude)  #mute if old clustering needed
  dplyr::filter(Meine_Clusters != "Exclude") %>%
  mutate(cluster = Meine_Clusters)  #choose which clustering is used



#---*------Reading in LCMS, MZmine files------------------

folder <- paste(admin, "AcidSinkhole/MSmine/20241106", sep = "/")
setwd(folder)

#CSV tables
files <- dir(folder, pattern = ".csv")
files

ms1peaks <- fread(files[2]) %>%
  #select(!contains("O2005314")) %>%
  select(!contains("O2005")) %>%
  dplyr::filter(rt < 16.9) %>%
  dplyr::filter(mz < 1000) %>%
  mutate(peakcheck = height/area) %>%
  select(id,rt,mz,height,area,peakcheck,contains("area")) %>%
  replace(is.na(.), 0)

ms2peaks <- fread(files[3]) %>%
  #select(!contains("O2005314")) %>%
  select(!contains("O2005")) %>%
  dplyr::filter(rt < 16.9) %>%
  dplyr::filter(mz < 1000) %>%
  mutate(peakcheck = height/area) %>%
  select(id,rt,mz,height,area,peakcheck,contains("area")) %>%
  replace(is.na(.), 0)


aligned_GNPS <- fread(files[1]) %>%
  select(`row ID`, `row m/z`, `row retention time`) %>%
  mutate(id = `row ID`) %>%
  left_join(ms2peaks %>% select(id, rt, mz))
  

#---*-------------Analysing MZmine tables---------------------

#PCA

dataname <- "Aligned, MS2, all features"
sample_feature <- ms2peaks %>%
  #dplyr::filter(id %in% importantfeatures$feature) %>%
  select(id, contains(":area")) %>%
  column_to_rownames("id") %>%
  #replace(is.na(.), 0) %>%
  #mutate_all(~ ifelse(. != 0, 1, 0)) %>%
  mutate_all(~ ifelse(. != 0, log10(.), 0))

#sample_feature <- t(apply(sample_feature, 1, function(x) x / max(x)))

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

#loadings
loadings <- as.data.frame(data_pca$rotation) %>%
  rownames_to_column("featureId")


#UMAP
set.seed(123)
umap <- umap(sample_feature)  #n_components = 3

umap_fp <- as.data.frame(umap$layout) %>%
  rownames_to_column("FileName")
umap_fp$FileName <- str_extract(umap_fp$FileName, "(?<=datafile:)[^.]+") 


#-------------Instrumental chemical space------------

featuregraph <- ms2peaks %>%
  select(id, rt, mz, contains(":area")) %>% #, contains(":area")
  #dplyr::filter(id %in% importantfeatures$feature) %>%
  gather(key = sample, value = "intensity", -c(rt, mz, id))
featuregraph$sample <- str_extract(featuregraph$sample, "(?<=datafile:)[^.]+")

featuregraph <- featuregraph %>%
  left_join(sampledata %>% mutate(sample = FileName) %>% select(sample, phcluster, p_h, cluster)) %>%
  na.omit() %>%
  dplyr::filter(intensity > 0) %>%
  dplyr::filter(sample %in% c("O2006040", "O2006076"))


#RT vs intensity slopes
ggplot(featuregraph, aes(x = rt, y = intensity, color = p_h)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_smooth(aes(group = cluster), method = "lm", se = FALSE, color = "black") + #linetype = factor(phcluster), group = cluster
  stat_cor(aes(group = cluster), method = "pearson", label.x = 7, label.y = 2) +
  labs(#title = "Scatter Plot of Precursor mass vs Intensity", 
       #x = "Precursor mass (mz)", 
       x = "Retention time (RT)",
       y = "Intensity", 
       color = "pH") +
  scale_color_viridis(option = "H", direction = -1) +
  scale_y_log10() +  # Use log scale if intensities span several orders of magnitude
  facet_wrap(~ cluster, nrow = 1) +
  #my_theme
  theme_minimal()


slopes <- featuregraph %>%
  group_by(sample) %>%
  summarize(
    slope = coef(lm(log10(intensity) ~ rt))[2]
  ) %>%
  left_join(featuregraph %>% select(sample, p_h, phcluster, cluster)) %>%
  unique()

ggplot(slopes, aes(y = slope, x = cluster, color = p_h)) +
  geom_point(alpha = 0.7, size = 4) +
  scale_color_viridis(option = "H", direction = -1) +
  #scale_color_manual(values = c("#990000","#0066CC","#339900","#FFCC33","#000066")) +
  theme_minimal() +
  labs(title = "slopes mz-log(Int)") +
  my_theme


##chemical space
#all peaks
data1 <- ms1peaks
data2 <- ms2peaks

#for only Acid Lakes
data1 <- ms1peaks %>%
  select(id, rt, mz, contains(":area")) %>%
  gather(key = sample, value = "intensity", -c(rt, mz, id))
data1$sample <- str_extract(data1$sample, "(?<=datafile:)[^.]+")

data1 <- data1 %>%
  left_join(sampledata %>% mutate(sample = FileName) %>% select(sample, phcluster)) %>%
  dplyr::filter(phcluster == 2)  %>%
  select(-phcluster) %>%
  dplyr::filter(intensity > 100) %>%
  spread(key = "sample", value = "intensity")

data2 <- ms2peaks %>%
  select(id, rt, mz, contains(":area")) %>%
  gather(key = sample, value = "intensity", -c(rt, mz, id))
data2$sample <- str_extract(data2$sample, "(?<=datafile:)[^.]+")

data2 <- data2 %>%
  left_join(sampledata %>% mutate(sample = FileName) %>% select(sample, phcluster)) %>%
  dplyr::filter(phcluster == 2)  %>%
  select(-phcluster) %>%
  dplyr::filter(intensity > 100) %>%
  spread(key = "sample", value = "intensity")


ggplot(data1, aes(x = rt, y = mz)) +
  geom_point(alpha = 0.6, size = 2, color = "lightgray") +
  geom_point(data = data2, aes(x = rt, y = mz), alpha = 0.6, color = "darkred", size = 2) + 
  #geom_smooth(aes(group = sample), method = "lm", se = FALSE, color = "black") + #linetype = factor(phcluster), group = cluster
  #stat_cor(aes(group = sample), method = "pearson", label.x = 7, label.y = 2) +
  labs(#title = "Scatter Plot of Precursor mass vs Intensity", 
    #x = "Precursor mass (mz)", 
    x = "Retention time (min)",
    y = "Ion mass", 
    color = "pH") +
  #scale_color_viridis(option = "H", direction = -1) +
  #scale_y_log10() +  # Use log scale if intensities span several orders of magnitude
  #facet_wrap(~ p_h + sample, nrow = 1) +
  my_theme
  #theme_minimal()


#heatmap

featuregraphspread <- featuregraph %>%
  arrange(rt) %>%
  select(id, sample, intensity) %>%
  spread(key = "id", value = "intensity")

heatmap(as.matrix(ms2standardsnorm),Rowv = NA, scale = "none", margins = c(5, 1)) 

#---*---------Metabolomics results from SIRIUS, averaged--------------

#ready tables
setwd(paste(admin, "AcidSinkhole/CalculatedFingerprints", sep = "/"))
dir(paste(admin, "AcidSinkhole/CalculatedFingerprints", sep = "/"), pattern = "20241106")

filename <- "20241106_canopus_canoprank_intensity_average.csv"


all_aver_classif <- fread(filename) %>%
  #dplyr::filter(sample != "O2005314") %>%
  #right_join(highph) %>%
  column_to_rownames("sample") %>%
  #select(one_of(selected_features)) %>%
  na.omit()


#eliminating zero rows and column and continuing with clustering
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

PC1var <- round(summary(data_pca)$importance[2,1]*100,1)
PC2var <- round(summary(data_pca)$importance[2,2]*100,1)
PC3var <- round(summary(data_pca)$importance[2,3]*100,1)

pca_fp <- as.data.frame(data_pca$x) %>%
  rownames_to_column("FileName")

loadings <- as.data.frame(data_pca$rotation) %>%
  #rownames_to_column("featureId") %>%
  rownames_to_column("Canop") %>%
  separate(Canop, into = c("Canop"), sep = "_")
  

#UMAP
set.seed(123)
umap <- umap(all_aver_classif)  #n_components = 3

umap_fp <- as.data.frame(umap$layout) %>%
  rownames_to_column("FileName")


#----------RDA----------------

#MS1 data
data <- sample_feature %>%
  rownames_to_column("FileName")
data$FileName <- str_extract(data$FileName, "(?<=datafile:)[^.]+") 
data <- sampledata %>%
  #mutate(cluster = Meine_Clusters) %>%  #mute if old clustering needed
  select(FileName, cluster) %>%
  left_join(data)

grouping <- factor(data$cluster)

data <- data %>%
  select(-cluster) %>%
  column_to_rownames("FileName")


#SIRIUS data  
data <- all_aver_classif %>%
  rownames_to_column("FileName") %>%
  right_join(sampledata %>% select(FileName, cluster))

grouping <- factor(data$cluster)

data <- data %>%
  select(-cluster) %>%
  column_to_rownames("FileName")


#rda
set.seed(123)
data_rda <- rda(data ~ grouping)
#data_rda <- rda(data)

#abiotic factors on top

colnames(sampledata)
abiotic <- data %>%
  rownames_to_column("FileName") %>%
  select(FileName) %>%
  left_join(sampledata) %>%
  select(p_h, depth_m, salinity_ppt, temperature_deg_c, phosphate_umol_kg, chlorophyll_ug_kg, oxygen_umol_kg,
         ammonium_umol_kg, nitrite_umol_kg, nitrate_umol_kg, transmissivity_percent)
colnames(abiotic) <- c("pH", "Depth", "Salinity", "Temp", "Phosphate", "Chlorophyll", "Oxygen", "Ammonium", "Nitrite", "Nitrate", "Transmissivity")

envfit_result <- envfit(data_rda ~ ., data = abiotic)
arrows_df <- as.data.frame(scores(envfit_result, "vectors")) %>%
  rownames_to_column("AbioticFactor")

# Scale the arrow coordinates
scaling_factor <- 5 # You can adjust this value to make the arrows longer
arrows_df <- arrows_df %>%
  mutate(RDA1 = RDA1 * scaling_factor,
         RDA2 = RDA2 * scaling_factor)


#back to calculating RDA with intensities
var_explained <- summary(data_rda)$cont$importance
RDA1var <- round(var_explained[2,1] * 100, 1)
RDA2var <- round(var_explained[2,2] * 100, 1)


rda_scores <- as.data.frame(scores(data_rda, display = "sites")) %>%
  rownames_to_column("FileName") %>%
  left_join(sampledata)

graphdata <- rda_scores


plot_pca <- ggplot(data = graphdata) +
  geom_point(mapping = aes(x = RDA1,
                           #shape = cluster,
                           y = RDA2),
             shape = 21, stroke = 1,
             size = 4, alpha = 0.5, color = "gray", data = graphdata[is.na(graphdata$p_h),]) + # NA points
  geom_point(mapping = aes(x = RDA1,
                           fill = cluster,
                           #shape = cluster,
                           y = RDA2),
             size = 4, alpha = 1, 
             shape = 21, color = "black", stroke = 0.8,
             data = graphdata[!is.na(graphdata$cluster),]) + # Non-NA points
  #my_theme +
  Meine_theme +
  #scale_color_viridis(option = "H", direction = -1) +
  #scale_color_manual(values = colors) +
  scale_color_manual(values = group_c) +
  scale_fill_manual(values = group_c) +
  ggtitle("CanopClasses, CraigClusters_noExclusion") +
  #facet_wrap(~ SampleName, nrow = 7) +
  xlab(paste("RDA1 (", RDA1var, "%)", sep = "")) +
  ylab(paste("RDA2 (", RDA2var, "%)", sep = ""))
plot_pca


plot_pca <- plot_pca +
  geom_segment(data = arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")),
               size = 0.5, color = "red") + # Set arrow color to red
  geom_text(data = arrows_df,
            aes(x = RDA1, y = RDA2, label = AbioticFactor),
            color = "red", vjust = -0.5, hjust = 0.5, size = 4) # Add text labels in red

plot_pca


#---*--------Graphs for comparing samples------------

graphdata <- sampledata %>%
  full_join(pca_fp, by = "FileName") %>%
  #mutate(silicate_umol_kg = as.numeric(gsub(",", ".", silicate_umol_kg))) %>%
  #mutate(p_h = as.numeric(gsub(",", ".", p_h))) %>%
  #mutate(depth_m = as.numeric(gsub(",", ".", depth_m))) %>%
  #mutate(Depth = as.numeric(gsub(",", ".", Depth))) %>%
  #drop_na(PC1) %>%
  #mutate(date = substr(FileName, 1, 5)) %>%
  drop_na(PC1)


plot_pca <- ggplot(data = graphdata) +
  geom_point(mapping = aes(x = PC1,
                           #shape = cluster,
                           y = PC2),
             size = 6, alpha = 0.5, color = "gray", data = graphdata[is.na(graphdata$p_h),]) + # NA points
  geom_point(mapping = aes(x = PC1,
                           colour = p_h,
                           #shape = cluster,
                           y = PC2),
             size = 6, alpha = 1, data = graphdata[!is.na(graphdata$p_h),]) + # Non-NA points
  my_theme +
  scale_color_viridis(option = "H", direction = -1) +
  ggtitle("Compound class, intensity") +
  #facet_wrap(~ SampleName, nrow = 7) +
  xlab(paste("PC1, ", PC1var, "%", sep = "")) +
  ylab(paste("PC2, ", PC2var, "%", sep = ""))
plot_pca

#Meine clusters
plot_pca <- ggplot(data = graphdata) +
  geom_point(mapping = aes(x = PC1,
                           #shape = cluster,
                           y = PC2),
             size = 6, alpha = 0.5, color = "gray", data = graphdata[is.na(graphdata$p_h),]) + # NA points
  geom_point(mapping = aes(x = PC1,
                           colour = cluster,
                           #shape = cluster,
                           y = PC2),
             size = 6, alpha = 1, data = graphdata[!is.na(graphdata$p_h),]) + # Non-NA points
  my_theme +
  #scale_color_viridis(option = "H", direction = -1) +
  #scale_color_manual(values = c("#990000","#0066CC","#339900","#FFCC33","#000066")) +
  scale_color_manual(values = colors) +
  ggtitle("Aligned, MS2 intensity") +
  #facet_wrap(~ SampleName, nrow = 7) +
  xlab(paste("PC1, ", PC1var, "%", sep = "")) +
  ylab(paste("PC2, ", PC2var, "%", sep = ""))
plot_pca


###UMAP

graphdata <- sampledata %>%
  full_join(umap_fp, by = "FileName") %>%
  drop_na(V1)

ggplot(data = graphdata) +
  geom_point(mapping = aes(x = V1,
                           #shape = date,
                           y = V2),
             size = 6, alpha = 0.5, color = "gray", data = graphdata[is.na(graphdata$depth_m),]) + # NA points
  geom_point(mapping = aes(x = V1,
                           colour = cluster,
                           #shape = date,
                           y = V2),
             size = 6, alpha = 1, data = graphdata[!is.na(graphdata$p_h),]) + # Non-NA points
  my_theme +
  #scale_color_viridis(option = "H", direction = -1) +
  scale_color_manual(values = c("#990000","#0066CC","#339900","#FFCC33","#000066")) +
  #facet_wrap(~ Injection_Type)
  ggtitle("Aligned, bothdays, CANOPUS average intensities, all LCMS features") +
  xlab(paste("UMAP1", sep = "")) +
  ylab(paste("UMAP2", sep = ""))



#loadings

biplot(data_pca, scale = 0, cex = 0.6)

N <- 50  # Number of top features to select
top_loadings <- loadings %>%
  mutate(abs_PC1 = abs(PC1), abs_PC2 = abs(PC2)) %>%
  rowwise() %>%
  mutate(max_abs_loading = max(abs_PC1, abs_PC2)) %>%
  ungroup() %>%
  arrange(desc(max_abs_loading)) %>%
  dplyr::slice(1:N) %>%
  #left_join(sir_canop_pos_Kai %>% select(Canop, name))
  left_join(ms2peaks %>% select(id, rt, mz) %>% mutate(featureId = as.character(id)))
  #left_join(top1annotations %>% select(featureId, adduct, `NPC#pathway`) %>% mutate(featureId = as.character(featureId))) %>%
  #replace(is.na(.), "unknown")

ggplot(top_loadings, aes(x = PC1, y = PC2, label = Canop)) +  #label = featureId, color =  rt #color = `NPC#pathway`
  geom_point() +
  geom_text(size = 3) +
  labs(title = "PCA Loadings Plot",
       x = paste0("PC1 (", PC1var, "% variance)"),
       y = paste0("PC2 (", PC2var, "% variance)")) +
  scale_color_viridis(option = "H") +
  my_theme
  theme_minimal()


#---------------Chemicals distribution in sample-------------

folder <- paste(admin, "LCMS_measured/Sinkhole/MSmine/FPtables", sep = "/")
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

sample_feature_cluster <- sample_feature %>% #sample_feature %>%
  rownames_to_column("FileName")
sample_feature_cluster$FileName <- str_extract(sample_feature_cluster$FileName, "(?<=datafile:)[^.]+")

sample_feature_cluster <- sampledata %>%
  select(FileName, phcluster) %>%
  left_join(sample_feature_cluster) %>%
  na.omit()

#finding variables
classes <- as.factor(sample_feature_cluster$phcluster)
df <- sample_feature_cluster %>%
  select(-phcluster, -FileName)
#Continue in code 20231103_RFE_Aris.R


#working with selected LCMS features
folder <- paste(admin, "AcidSinkhole/CalculatedFingerprints", sep = "/")
setwd(folder)
importantfeatures <- fread("20241106_rfe_ms2log10int_100var.csv") %>%
  #na.omit() %>%
  #dplyr::filter(count > 5) %>%
  mutate(mean = round(mean, digits=2)) %>%
  #mutate(FoldMean = rowMeans(select(., contains("Fold")), na.rm = TRUE)) %>%
  mutate(featureId = as.character(var))


#important features annotated
folder <- paste(admin, "AcidSinkhole/RESULTS/Tables", sep = "/")
setwd(folder)
importantfeatures <- fread("20241106_all_important_features_intensities.csv") %>%
  #select(featureId) %>%
  mutate(featureId = as.character(featureId)) %>%
  dplyr::filter(countFold >= 5) %>%
  unique() 

#matched with InChIKey2D
folder <- paste(admin, "LCMS_measured/Sinkhole/RESULTS/Tables", sep = "/")
setwd(folder)
importantfeatures <- fread("matches_1000.csv") %>%
  select(featureId) %>%
  mutate(featureId = as.character(featureId)) %>%
  unique() 


#Selected features from other steps, lists

feature_analysis <- sample_feature %>%
  rownames_to_column("FileName") %>%
  #select(FileName, any_of(as.character(canop_compounds$featureId)))
  select(c(FileName, "30653", "31967"))
  #select(one_of(top_loadings$featureId)) %>%
  select(FileName, any_of(as.character(comparison$featureId)))
feature_analysis$FileName <- str_extract(feature_analysis$FileName, "(?<=datafile:)[^.]+")


feature_analysis <- feature_analysis %>%
  gather(key = "featureId", value = "logInt", -FileName) %>%
  left_join(sampledata %>% select(FileName, p_h, phcluster, cluster)) %>%
  # left_join(top1annotations %>% 
  #             select(featureId, InChIkey2D, ConfidenceScore) %>% 
  #             mutate(featureId = as.character(featureId))) %>%
  # left_join(massFromNTS %>%
  #             mutate(featureId = as.character(id)) %>% select(featureId, cpdID))
  #mutate(p_h = as.numeric(gsub(",", ".", p_h))) %>%
  #left_join(cluster %>% select(FileName, cluster)) %>%
  na.omit() %>%
  left_join(ms2peaks %>% 
              mutate(mz = round(mz, digits = 4)) %>%
              mutate(rt = round(rt, digits = 1)) %>%
              mutate(featureId = as.character(id)) %>% 
              select(featureId, rt, mz)) %>%
  left_join(iodides %>% 
              dplyr::filter(significant == "Significant for 2") %>%
              dplyr::filter(log2FoldChange > 5) %>%
              select(featureId, molecularFormula) %>%
              mutate(featureId = as.character(featureId))) %>%
  na.omit()



#Graphs!
var_plot <- feature_analysis %>%
  #mutate(id = as.numeric(featureId)) %>%
  #left_join(ms2peaks %>% select(id, rt, mz)) %>%
  #left_join(targets_fromNTS %>% select(feature, name)) %>%
  #left_join(cluster_assignment, by = "feature") %>%
  #mutate(mean = round(mean, digits=2)) %>%
  #dplyr::filter(FoldMean > 1.5) %>%
  #dplyr::filter(featureId %in% c(24181,17629,28244,35471)) %>%
  #dplyr::filter(FeatureCluster == 5) %>%
  ggplot(aes(x = cluster, 
             y = logInt)) +
  geom_point(aes(x = cluster, y = logInt, fill = cluster),
            size = 5, position = position_jitter(width = 0.3),
            shape = 21, stroke = 1
            ) +
  my_theme +
  ggtitle("Enriched in Acid Lake, Log2 Fold Change > 5") +
  #scale_color_viridis(option = "H", direction = -1) +
  scale_fill_manual(values = group_c) +
  facet_wrap(~ featureId + mz + rt, scales = "free_y", nrow = 3) + #scales = "free_y"
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1))
var_plot


#annotating 
setwd(paste(admin, "LCMS_measured/Sinkhole/RESULTS/Figures", sep = "/"))

annotated <- feature_analysis %>%
  select(feature) %>%
  unique() %>%
  mutate(featureId = as.numeric(feature)) %>%
  left_join(top1annotations)
write_delim(annotated, "most_important_features_structures.csv", delim = ";")



#boxplot
var_plot <- feature_analysis %>%
  ggplot(aes(x = cluster, 
             y = logInt)) +
  # geom_violin(aes(fill = cluster), 
  #             scale = "width", width = 0.7) +   # Violin plot
  geom_boxplot(width = 0.2, 
               outlier.shape = NA, 
               color = "orange", alpha = 0.6) +  # Boxplot
  geom_jitter(width = 0.05,
              size = 1.2,
              alpha = 0.7,
              color = p_h) +
  scale_color_viridis(option = "H", direction = -1) +
  facet_wrap(~ feature, scales = "free_y", nrow = 4) + #scales = "free_y"
  theme(legend.position="none") +
  my_theme  

var_plot


#----------Important SIRIUS classes for separating sinkholes---------------


#written out file

setwd(paste(admin, "AcidSinkhole/RESULTS/Tables", sep = "/"))
feature_analysis <- fread("20241106_canopusclass_intens_averages_all.csv") %>%
  na.omit() %>%
  dplyr::filter(countFold >=1)


#high in sinkhole
highsinkhole <- feature_analysis %>%
  group_by(phcluster, Canop) %>%
  mutate(aver_cluster = mean(average)) %>%
  ungroup() %>%
  select(Canop, phcluster, aver_cluster) %>%
  mutate(phcluster = paste("ph",phcluster, sep = "")) %>%
  unique() %>%
  spread(key = "phcluster", value = aver_cluster) %>%
  mutate(difference = log2(ph2/ph1)) %>%
  dplyr::filter(difference < 2)

highsinkholehalogen <- highsinkhole %>%
  dplyr::filter(Canop %in% c("Canop35", "Canop1028", "Canop1517", "Canop3092", "Canop483", "Canop4801", "Canop4800", "Canop3096", "Canop3094", "Canop3881")) %>%
  select(Canop)

feature_analysis <- feature_analysis %>%
  dplyr::filter(Canop %in% highsinkhole$Canop) %>%
  select(-cluster) %>%
  left_join(sampledata %>% select(FileName, cluster))


names <- feature_analysis %>%
  select(Canop, name, meanFold) %>% unique()

#from original data#frmeanFoldom original data
cc_cluster <- all_aver_classif %>%
  rownames_to_column("FileName") %>%
  left_join(sampledata %>%
              select(FileName, phcluster, cluster)) %>%
  na.omit()

classes <- as.factor(cc_cluster$phcluster)
df <- cc_cluster %>%
  select(-phcluster, -FileName)


#meine cluster color

var_plot <- feature_analysis %>%
  dplyr::filter(Canop %in% c("Canop4800", "Canop3094", "Canop226")) %>%
  ggplot(aes(x = p_h, 
             y = average)) +
  geom_point(aes(x = p_h, y = average, color = cluster),
             size = 4) +
  my_theme +
  #scale_color_viridis(option = "H", direction = -1) +
  scale_color_manual(values = c(colors)) +
  #theme_light() +
  ggtitle("CANOPUS aligned, 2x higher in sinkhole") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~ meanFold+Canop+name, scales = "free_y", ncol = 7) + #scales = "free_y"
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1))
var_plot


setwd(paste(admin, "AcidSinkhole/RESULTS/Figures", sep = "/"))
svg("canopus_high_in_sinkhole.svg",width = 20, height = 20, pointsize = 10)
print(var_plot)
dev.off()


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


# Distribution per sample type for each Canopus class

feature_analysis_canopus <- feature_analysis %>%
  select(-cluster) %>%
  left_join(sampledata %>% select(FileName, cluster)) %>%
  group_by(Canop) %>%
  mutate(canop_aver = mean(average)) %>%
  ungroup()

canop_order <- feature_analysis_canopus %>%
  select(Canop, canop_aver) %>%
  unique() %>%
  arrange(desc(canop_aver)) %>%
  rownames_to_column("canop_order")

feature_analysis_canopus <- feature_analysis_canopus %>%
  left_join(canop_order) %>%
  mutate(canop_order = as.numeric(canop_order))

var_plot <- feature_analysis_canopus %>%
  #dplyr::filter(cluster == "Deep Acidic") %>%
  #dplyr::filter(Canop %in% c("Canop4800", "Canop3094", "Canop226")) %>%
  ggplot(aes(x = canop_order, 
             y = average)) +
  geom_point(aes(x = canop_order, y = average, color = cluster),
             size = 4, alpha = 1) +
  my_theme +
  #scale_color_viridis(option = "H", direction = -1) +
  scale_color_manual(values = c(colors)) +
  #theme_light() +
  theme(legend.position="none") +
  ggtitle("CANOPUS") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~ cluster, ncol = 5) + #scales = "free_y"
  theme(panel.grid.minor = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1))
var_plot


#-------------FP and CANOPUS from SIRIUS-----------------

# #getting SIRIUS fingerprints into table
folderwithSIRIUSfiles <- "C:/Users/pille/Desktop/SIRIUS/sinkhole/20241106/part5"
setwd(folderwithSIRIUSfiles)

CCresults <- CCtable_SIRIUS5(folderwithSIRIUSfiles)
FPresults <- MFtable_SIRIUS5(folderwithSIRIUSfiles)
setwd(paste(admin, "LCMS_measured/Sinkhole/MSmine/FPtables", sep = "/"))
write_delim(FPresults, "20241106_fingerprint_allranks_part5.csv", delim = ";")
write_delim(CCresults, "20241106_canopus_allranks_part5.csv", delim = ";")


fppath <- paste(admin, "AcidSinkhole/MSmine/FPtables", sep = "/")
setwd(fppath)
files <- dir(fppath, pattern = "20241106_canopus_can")
files

all_sirius <- tibble()
for (file in files){
  sirius <- fread(file)
  all_sirius <- all_sirius %>%
    bind_rows(sirius)
}

check <- all_sirius %>% 
  select(featureId) %>% unique()

all_sirius_rank1f <- all_sirius %>%
  left_join(top1annotations %>% select(id, adduct, precursorFormula, `NPC#pathway`)) %>%
  select(id, adduct, precursorFormula, `NPC#pathway`, everything()) %>%
  dplyr::filter(!is.na(`NPC#pathway`))

all_sirius_rank1f <- all_sirius %>%
  group_by(id) %>%
  arrange(id, formulaRank) %>%
  mutate(newRank = dense_rank(formulaRank)) %>%
  ungroup() %>%
  dplyr::filter(newRank == 1) %>%  #some rows did not have rank 1 existing
  select(-newRank)

#filter m+h for those with several adducts

all_sirius_rank1f <- all_sirius_rank1f %>%
  group_by(id) %>%
  arrange(id, desc(adduct == "[M+M]+")) %>% # Arrange by preference for [M+M]+
  slice(1) %>% # Take the first row of each group
  ungroup()

check <- all_sirius_rank1f %>% select(featureId) %>% unique()

# setwd(paste(admin, "LCMS_measured/Sinkhole/MSmine/FPtables", sep = "/"))
# write_delim(all_sirius_rank1f, "20241106_canopus_canoprank.csv", delim = ";")

setwd(paste(admin, "AcidSinkhole/MSmine/FPtables", sep = "/"))
all_mfp <- fread("20241106_fingerprint_canoprank.csv")


#-------------FP and CANOPUS from SIRIUS6-----------------

fpfolder <- paste(admin, "AcidSinkhole/SIRIUS/SIRIUS6_20241106/results", sep = "/")
setwd(fpfolder)

files <- dir(fpfolder, pattern = "API", include.dirs = TRUE, recursive = TRUE)
files

all_sirius <- tibble()
for (file in files){
  sirius <- fread(file)
  all_sirius <- all_sirius %>%
    bind_rows(sirius)
}

files <- dir(fpfolder, pattern = "canopus_formula_summary", include.dirs = TRUE, recursive = TRUE)
files

formularank <- tibble()
for (file in files){
  sirius <- fread(file) %>%
    select(mappingFeatureId, molecularFormula, adduct, formulaRank, alignedFeatureId, formulaId) %>%
    group_by(alignedFeatureId) %>%
    mutate(newrank = min(formulaRank)) %>%
    ungroup() %>%
    mutate(chosen = formulaRank == newrank) %>%
    dplyr::filter(chosen == TRUE)
  formularank <- formularank %>%
    bind_rows(sirius)
}

rank1form_mfp <- formularank %>%
  left_join(all_sirius) %>%
  select(mappingFeatureId, molecularFormula, adduct, contains("Un")) %>%
  unique()


#structurerank

files <- dir(fpfolder, pattern = "structure_identifications", include.dirs = TRUE, recursive = TRUE)
files
file <- files[3]

formularank <- tibble()
for (file in files){
  sirius <- fread(file) %>%
    #select(mappingFeatureId, molecularFormula, adduct, formulaRank, alignedFeatureId, formulaId) %>%
    #mutate(chosen = formulaRank == newrank) %>%
    dplyr::filter(structurePerIdRank == 1)
  formularank <- formularank %>%
    bind_rows(sirius)
}

rank1form_mfp <- formularank %>%
  left_join(all_sirius) %>%
  select(mappingFeatureId, molecularFormula, adduct, contains("Un")) %>%
  unique()



#-----------Getting averaged fp from tables---------------

fp_folder <- paste(admin, "AcidSinkhole/MSmine/FPtables", sep = "/")
setwd(fp_folder)
dir(fp_folder, pattern = "20241106")


sirius_data <- fread("20241106_canopus_canoprank.csv")


ms2peakspresence <- ms2peaks %>%
  #select(id, mz) %>%
  mutate(featureId = id) %>%
  select(featureId, contains("datafile:"))

ms2peakspresence <- gather(ms2peakspresence, key = "sample", value = "intensity", -featureId) %>%
  dplyr::filter(intensity > 0) %>%
  na.omit()

ms2peakspresence$sample <- str_extract(ms2peakspresence$sample, "(?<=datafile:)[^.]+")

ms2peakspresence <- ms2peakspresence %>%
  left_join(sirius_data %>% select(featureId, contains("Canop")), 
            by = "featureId") %>%
  na.omit() %>%  #until here for regular average
  group_by(sample) %>%
  mutate(intsum = sum(intensity)) %>%
  ungroup() %>%
  mutate(intnorm = intensity/intsum) %>%
  select(featureId, sample,intensity, intsum, intnorm, everything())


peakcount <- ms2peakspresence %>%
  select(featureId, sample, intensity) %>%
  group_by(sample) %>%
  summarise(count = n()) %>%
  ungroup()


# get table with intensities 
ms2peaksintens <- ms2peakspresence %>%
  mutate_at(vars(starts_with("Canop")), ~. * intnorm) %>%
  select(featureId, sample, starts_with("Canop"))

#average by summing the scaled intensities
averages <- ms2peaksintens %>%
  select(sample, contains("Canop")) %>%
  group_by(sample) %>%                        
  summarise_at(vars(contains("Canop")),
               list(av = sum)) %>%
  ungroup()


#just averages for binary qualitative
averages <- ms2peakspresence %>%  
  select(sample, contains("Canop")) %>%
  group_by(sample) %>%                        
  summarise_at(vars(contains("Canop")),
               list(av = mean)) %>%
  ungroup()


setwd(paste(admin, "LCMS_measured/Sinkhole/CalculatedFingerprints", sep = "/"))
#write_delim(averages, "20241106_canopus_canoprank_intensity_average.csv", delim = ";")

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


sir_folder <- paste(admin, "AcidSinkhole/SIRIUS/20241106", sep = "/")
setwd(sir_folder)


#Rank 1 structures only  #compound_identifications, canopus_compound_summary
comp_ind_files <- dir(sir_folder, pattern = "compound_identifications.tsv", recursive = TRUE)
comp_ind_files
#file <- comp_ind_files[7]

top1annotations <- tibble()
for(file in comp_ind_files){
  id_file <- fread(file) %>%
    mutate(ConfidenceScore = as.numeric(ConfidenceScore))
  
  top1annotations <- top1annotations %>%
    bind_rows(id_file) %>%
    unique()
}

top1annotations <- top1annotations %>%
  mutate(adduct = gsub(" ", "", adduct)) %>%
  select(featureId, everything()) %>%
  mutate(RT = retentionTimeInSeconds/60) %>%
  select(confidenceRank,featureId,RT,ionMass,smiles,name, InChIkey2D, ConfidenceScore, molecularFormula, adduct, xlogp) %>%
  # right_join(importantfeatures %>%
  #              mutate(featureId = as.numeric(feature)) %>%
  #              select(featureId, FoldMean)) %>%
  drop_na(name) %>%
  dplyr::filter(featureId %in% ms2peaks$id)


sirius_structure_annotation <- top1annotations


#formula
#Rank 1 structures only  #compound_identifications, canopus_compound_summary
comp_ind_files <- dir(sir_folder, pattern = "formula_identifications.tsv", recursive = TRUE)
comp_ind_files
#file <- comp_ind_files[7]

top1annotations <- tibble()
for(file in comp_ind_files){
  id_file <- fread(file) 
  
  top1annotations <- top1annotations %>%
    bind_rows(id_file) %>%
    unique()
}

top1annotations <- top1annotations %>%
  mutate(adduct = gsub(" ", "", adduct)) %>%
  select(featureId, everything()) %>%
  mutate(RT = retentionTimeInSeconds/60) %>%
  select(featureId,RT,ionMass, molecularFormula)

sirius_formula_annotation <- top1annotations


#formula
#Rank ALL formulas  #compound_identifications, canopus_compound_summary
comp_ind_files <- dir(sir_folder, pattern = "formula_identifications_all.tsv", recursive = TRUE)
comp_ind_files
#file <- comp_ind_files[7]

top1annotations <- tibble()
for(file in comp_ind_files){
  id_file <- fread(file) 
  
  top1annotations <- top1annotations %>%
    bind_rows(id_file) %>%
    unique()
}

top1annotations <- top1annotations %>%
  mutate(adduct = gsub(" ", "", adduct)) %>%
  select(featureId, everything()) %>%
  mutate(RT = retentionTimeInSeconds/60) %>%
  select(featureId,RT,ionMass, molecularFormula)

sirius_formula_annotation_all <- top1annotations





ggplot(data = top1annotations) +
  geom_point(mapping = aes(x = RT,
                           colour = confidenceRank,
                           y = xlogp),
             size = 3, alpha = 0.7) + 
  geom_smooth(aes(x=RT, y=xlogp), method = "lm", se = FALSE, color = "black") +
  stat_cor(aes(x = RT, y = xlogp), 
           method = "pearson", 
           label.x = 2, 
           label.y = -11) +
  my_theme



all_annotations_targets <- top1annotations %>%
  inner_join(targetlist, by = "InChIkey2D") %>%
  select(featureId, cpdID, name, InChIkey2D, ConfidenceScore, molecularFormula, adduct)


#getting all annotated structures, filter to top50

sir_folder <- paste(admin, "AcidSinkhole/SIRIUS/20241106", sep = "/")
setwd(sir_folder)

#Rank 1-50 structures
comp_ind_files <- dir(sir_folder, pattern = "compound_identifications_all.tsv", recursive = TRUE)
comp_ind_files
#file <- comp_ind_files[1]

top50annotations <- tibble()
for(file in comp_ind_files){
  id_file <- fread(file)  %>%
    mutate(ConfidenceScore = as.numeric(ConfidenceScore)) %>%
    dplyr::filter(structurePerIdRank <= 1000)
  
  top50annotations <- top50annotations %>%
    bind_rows(id_file) %>%
    unique()
}


top50annotations <- top50annotations %>%
  mutate(adduct = gsub(" ", "", adduct)) %>%
  select(featureId, everything()) %>%
  mutate(RT = retentionTimeInSeconds/60) %>%
  select(confidenceRank,featureId,RT,ionMass,smiles,name, InChIkey2D, ConfidenceScore, molecularFormula, adduct, xlogp) %>%
  # right_join(importantfeatures %>%
  #              mutate(featureId = as.numeric(feature)) %>%
  #              select(featureId, FoldMean)) %>%
  drop_na(name) %>%
  dplyr::filter(featureId %in% ms2peaks$id)

sirius_structure_annotation <- top1annotations




all_annotations_targets <- top50annotations %>%
  inner_join(targetlist, by = "InChIkey2D") %>%
  select(featureId, cpdID, name, InChIkey2D, ConfidenceScore, molecularFormula, adduct) %>%
  mutate(id = featureId) %>%
  select(id, cpdID, name, InChIkey2D,adduct) %>%
  #left_join(targetmassfile %>% select(cpdID, adduct, mz_calc), by = "cpdID") %>%
  left_join(ms2peaks %>% select(id, contains(":area")))

#change file name to short

colnames(all_annotations_targets) <- gsub("datafile:([^.]+)\\.mzXML.*", "\\1", colnames(all_annotations_targets))




count <- all_annotations_targets %>%
  select(id) %>%
  unique()

setwd(paste(admin, "AcidSinkhole/RESULTS/Tables", sep = "/"))

data <- fread("most_important_features_structures.csv") %>%
  select(featureId, InChIkey2D, everything()) 

#write_delim(all_annotations_targets, "20241106_target_structure_match_top1000.csv", delim = ";")
#write_delim(massFromNTS, "20241106_all_mass_match_ms2.csv", delim = ";")


#---------GNPS Molecular Network results------------------

folder <- paste(admin, "AcidSinkhole/GNPS/20241106/20250402_GNPS2", sep = "/")
setwd(folder)

gnps_annot <- fread("gnps2_library_results.tsv")
clean_names <- function(names_vec) {
  names_vec <- gsub("[^A-Za-z0-9]+", "_", names_vec)  # Replace special characters with _
  names_vec <- gsub("^_+|_+$", "", names_vec)  # Remove leading/trailing underscores
  names_vec <- make.names(names_vec, unique = TRUE)  # Ensure valid R names
  return(names_vec)
}

setnames(gnps_annot, clean_names(names(gnps_annot)))

gnps_annot <- gnps_annot[!grepl("negative", IonMode, ignore.case = TRUE)]

gnps_annot[, ranking := frank(-MQScore, ties.method = "min"), by = Scan]

gnps_annot <- gnps_annot[ranking == 1]

gnps_annot <- gnps_annot[, .(
  Scan, Compound_Name, SharedPeaks, Adduct,
  InChIKey, InChIKey_Planar,
  MQScore, superclass, npclassifier_pathway,ranking,
  SpecMZ, Precursor_MZ
)]

count <- sirius_gnps_match %>% select(Scan) %>% unique()


compare <- ms2peaks %>%
  mutate(featureId = id) %>%
  select(featureId, rt) %>%
  full_join(gnps_annot %>% mutate(featureId = Scan)) %>%
  full_join(sirius_structure_annotation %>%
              select(featureId, InChIkey2D, name)) %>%
  full_join(sirius_class_annotation) %>%
  #na.omit() %>%
  unique()

clean_names_symbols <- function(names_vec) {
  names_vec <- gsub("[^A-Za-z0-9 ]", "", names_vec)  # Remove symbols but keep letters, numbers, and spaces
  return(names_vec)
}
setnames(compare, clean_names_symbols(names(compare)))

colnames(compare)

compare[, sameclass := ifelse(superclass == ClassyFiresuperclass, "yes", "no")]
compare[, sameNPclass := ifelse(npclassifierpathway == NPCpathway, "yes", "no")]
compare[, samestructure := ifelse(InChIKeyPlanar == InChIkey2D, "yes", "no")]


sirius_gnps_match <- compare[samestructure == "yes"]


#----------Analysing instrument blanks----------

folder <- paste(admin, "LCMS_measured/Sinkhole/MSmine/20240809_instrumentblank", sep = "/")
setwd(folder)

#CSV tables
files <- dir(folder, pattern = ".csv")
files

ms1peaks <- fread(files[2]) %>%
  dplyr::filter(rt < 15) %>%
  mutate(peakcheck = height/area) %>%
  dplyr::filter(height > 65000000) %>%
  select(id,rt,mz,height,area,peakcheck,contains(":area")) %>%
  select(!contains("O2005"))


intenscheck <- ms1peaks %>%
  select(id,contains(":area")) %>%
  column_to_rownames("id")
intenscheck <- t(apply(intenscheck, 1, function(x) x / max(x)))
intenscheck <- as.data.frame(intenscheck)

intenscheck <- intenscheck %>%
  rownames_to_column("id") %>%
  gather(key = "sample", value = "intensity", -id) #%>%
  #dplyr::filter(grepl("O2006", sample))
intenscheck$sample <- str_extract(intenscheck$sample, "(?<=datafile:)[^.]+")  

intenscheck_changes <- intenscheck %>%
  group_by(id) %>%
  arrange(sample) %>%
  mutate(change_in_intensity = intensity - lag(intensity)) %>%
  dplyr::filter(!is.na(change_in_intensity))

ggplot(intenscheck_changes, aes(x = sample, y = change_in_intensity, color = as.factor(id), group = id)) +
  geom_line() + 
  #facet_wrap(~ id, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #geom_smooth(aes(group = id), method = "lm", se = FALSE, color = "black")+
  theme(legend.position="none")


#---------------Boxplots with p_values------------

library(ggplot2)
library(tidyverse)
library(ggpubr)

my_comp <- feature_analysis %>%
  select(cluster, logInt, FileName, feature, p_h)

unique_clusters <- unique(my_comp$cluster)

# Perform pairwise comparison between the two clusters
comparisons <- list(c(unique_clusters[1], unique_clusters[2]))

# Create the boxplot with actual p-values
ggboxplot(my_comp,
          x = "cluster", 
          y = "logInt",
          fill = "cluster", 
          palette = "Dark2") +
  facet_wrap(~ feature, scales = "free_y", nrow = 4) +
  stat_compare_means(label = "p.format",    # Use p.format for actual p-value display
                     comparisons = comparisons,
                     method = "t.test") +  # Remove symnum.args to show actual p-values
  theme(legend.position = "none") +
  theme_light()




#----------Volcano plots---------------

mslong <- ms2peaks %>%
  pivot_longer(cols = contains(":area"), names_to = "FileName",values_to = "intensity") %>%
  mutate(
    intensity = as.numeric(intensity),           # ensure numeric
    intensity = if_else(intensity == 0, 100, intensity)
  )
mslong$FileName <- str_extract(mslong$FileName, "(?<=datafile:)[^.]+")

mslong <- mslong %>%
  right_join(sampledata %>% select(FileName, phcluster)) %>%
  select(-FileName)

#volcano for sinkhole vs other
df <- mslong

data_groups <- split(df, df$phcluster)
data_keys <- names(data_groups)

i = 1
data_keys[i]
j = 2
data_keys[j]

label1 <- data_keys[i]
label2 <- data_keys[j]
group1 <- data_groups[[data_keys[i]]] %>% select(id, intensity)
group2 <- data_groups[[data_keys[j]]] %>% select(id, intensity)

#NEW!!!!!
group1 <- data_groups[[data_keys[i]]] %>% select(id, intensity) %>% mutate(group = label1)
group2 <- data_groups[[data_keys[j]]] %>% select(id, intensity) %>% mutate(group = label2)
#NEW END!!!!

pseudocount <- 0  # Small value to avoid division by zero

#comparison <- full_join(group1, group2, by = c("id"), suffix = c("_g1", "_g2")) %>%
comparison <- bind_rows(group1, group2) %>%  #NEW!!!!!!
  group_by(id) %>%
  
  dplyr::mutate(
    mean_intensity_g1 =
      10^(mean(log10(intensity[group == "1"] + pseudocount))),
    mean_intensity_g2 =
      10^(mean(log10(intensity[group == "2"] + pseudocount))),
    
    p_value = tryCatch(
      t.test(
        log10(intensity[group == "1"] + pseudocount),
        log10(intensity[group == "2"] + pseudocount),
        var.equal = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
  ) %>%
  ungroup() %>%
  select(-intensity, -group) %>% unique() %>%
  
  mutate(
    log2FoldChange = log2((mean_intensity_g2 + pseudocount) / (mean_intensity_g1 + pseudocount)),
    adj_p_value = p.adjust(p_value, method = "BH"),
    significant = case_when(
      adj_p_value > 0.05 & abs(log2FoldChange) < 1 ~ "Insignificant",
      adj_p_value < 0.05 & log2FoldChange < -1 ~ paste("Significant for", label1, sep = " "),
      adj_p_value < 0.05 & log2FoldChange > 1 ~ paste("Significant for", label2, sep = " "),
      TRUE ~ "Insignificant")) %>%
  #some p-values got zero value fro some reason, replace with lowest value
  
  mutate(exclusive_to_one_group = 
           (abs(mean_intensity_g1 - 100) < 1e-6 |
              abs(mean_intensity_g2 - 100) < 1e-6)) %>%

    # fortext = case_when(
    #   adj_p_value < 0.05 & abs(log2FoldChange) > 4 ~ "Significant",  # Adjusted to consider more stringent criteria
    #   TRUE ~ "Insignificant"  # Default to "Insignificant" if no other condition is met
    # ),
  left_join(ms2peaks %>% select(id, rt, mz))


# comp_selected <- comparison %>%
#   dplyr::filter(id %in% c(17629, 24181, 28244, 49197, 63718))

# 
comparison <- comparison %>%
  # select(-intensity_g1, -intensity_g2) %>%
  # unique() %>%
  left_join(sirius_structure_annotation %>%
              mutate(id = featureId) %>%
              mutate(SIRIUSannotation = "yes") %>%
              select(id, SIRIUSannotation, molecularFormula, InChIkey2D, ConfidenceScore)) %>%
  mutate(SIRIUSannotation = replace_na(SIRIUSannotation, "no")) %>%
  # left_join(sirius_formula_annotation %>%
  #             mutate(id = featureId,
  #                    formula = molecularFormula) %>%
  #             select(id, formula)) %>%
  left_join(rankAny_halides)


insignificant <- comparison %>%
  dplyr::filter(significant == "Insignificant")


# Generate volcano plot
plot <- ggplot(comparison, aes(x = log2FoldChange, y = -log10(p_value))) +
  geom_point(aes(color = mz, shape = exclusive_to_one_group), alpha = 0.8, size = 3) +
  geom_point(data = insignificant, shape = 21, fill = "lightgray", color = "gray", size = 3, stroke = 0.5) +
  #geom_point(data = comp_selected, shape = 21, fill = "lightgray", color = "red", size = 5, stroke = 0.5) +
  # geom_text_repel(aes(label = ifelse(!is.na(rt), as.character(id), "")),
  #                 max.overlaps = 10, size = 2) +
  scale_color_viridis(option = "H", direction = 1) +
  #scale_color_manual(values = c("red")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +# Different shape for exclusive features
  #scale_shape_manual(values = c("yes" = 16, "no" = 17)) +
  labs(
    #title = paste(label1, "vs", label2),
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    #color = "Mass-to-charge ratio (m/z)",
    color = "Retention time (min)",
    shape = "Exclusive to group"
  ) +
  #facet_wrap(~ SIRIUSannotation , nrow = 2) +
  Meine_theme
  #my_theme
plot


#RT and mz comparison from volcano plot

setwd(paste(admin, "AcidSinkhole/RESULTS/Tables", sep = "/"))
comparison <- fread("20260127_volcano_all.csv")
  # dplyr::filter(abs(log2FoldChange) > 5) %>%
  # dplyr::filter(-log10(p_value) > 5) %>%
  # dplyr::filter(log2FoldChange > 5)

rt_importance <- comparison %>%
  select(id, significant, rt, mz)%>%
  mutate(significant = as.factor(significant))

lvls <- levels(rt_importance$significant)
my_comparisons <- combn(lvls, 2, simplify = FALSE)

#p-values

# Pairwise Welch t-tests (+ BH adjustment) via rstatix
pairwise_results <- rt_importance %>%
  t_test(rt ~ significant, var.equal = FALSE) %>%  # Welch t-test
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  # Optional: rename/format columns for clarity
  rename(
    group1 = group1,
    group2 = group2,
    p_value = p,
    p_value_adj_BH = p.adj
  )



p <- ggplot(rt_importance, aes(x = significant, y = rt, fill = significant)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(alpha = 0.35, width = 0.15, size = 1) +
  labs(
    x = "Significance group",
    y = "Retention time (min)",
    title = "mz comparison across significance groups"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  expand_limits(y = max(rt_importance$rt, na.rm = TRUE) * 1.12) + # space for brackets
  
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    method.args = list(var.equal = FALSE),  # Welch t-test
    label = "p.format",                     # << numeric p-values
    p.adjust.method = "BH",
    tip.length = 0.01,
    step.increase = 0.06
  ) +
  scale_fill_manual(values = c("lightgray", "#58b2e9", "#e79e02")) +
  Meine_theme
  
  # stat_compare_means(
  #   method = "anova",
  #   label = "p.format",
  #   label.y = max(rt_importance$rt, na.rm = TRUE) * 1.10
  # )
 p
setwd(paste(admin, "AcidSinkhole/RESULTS/Figures", sep = "/"))
#ggsave("ttest_enriched_rt.svg", plot = p, device = svglite::svglite, width = 6, height = 6)




iodides <- comparison %>%
  mutate(iod = grepl("F", formula)) %>%
  select(id, significant, molecularFormula, formula, iod, InChIkey2D, SIRIUSannotation, log2FoldChange)



#----------volcano with compound classes-------------

df <- all_aver_classif %>%
  rownames_to_column("FileName") %>%
  right_join(sampledata %>% select(FileName, phcluster)) %>%
  gather(key = "Canop", value = "average", -c(phcluster, FileName))

data_groups <- split(df, df$phcluster)
data_keys <- names(data_groups)

i = 1
data_keys[i]
j = 2
data_keys[j]

label1 <- data_keys[i]
label2 <- data_keys[j]
# group1 <- data_groups[[data_keys[i]]] %>% select(Canop, average)
# group2 <- data_groups[[data_keys[j]]] %>% select(Canop, average)
group1 <- data_groups[[data_keys[i]]] %>% select(Canop, average) %>% mutate(group = label1)
group2 <- data_groups[[data_keys[j]]] %>% select(Canop, average) %>% mutate(group = label2)

pseudocount <- 0.000000001  # Small value to avoid division by zero

comparison <- bind_rows(group1, group2) %>%  #NEW!!!!!!
  group_by(Canop) %>%
  dplyr::mutate(
    mean_intensity_g1 =
      mean(average[group == "1"], na.rm = TRUE),
    mean_intensity_g2 =
      mean(average[group == "2"], na.rm = TRUE),
      
    p_value = tryCatch(
      t.test(
        (average[group == "1"] + pseudocount),
        (average[group == "2"] + pseudocount),
        var.equal = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
  ) %>%
  ungroup() %>%
  select(-average, -group) %>% unique() %>%
  
  #mutate(test = log2((mean_intensity_g2+pseudocount)/(mean_intensity_g1+pseudocount)))
  
  mutate(
    log2FoldChange = log2((mean_intensity_g2 + pseudocount) / (mean_intensity_g1 + pseudocount)),
    adj_p_value = p.adjust(p_value, method = "BH"),
    significant = case_when(
      adj_p_value > 0.05 & abs(log2FoldChange) < 1 ~ "Insignificant",
      adj_p_value < 0.05 & log2FoldChange < -1 ~ paste("Significant for", label1),
      adj_p_value < 0.05 & log2FoldChange > 1 ~ paste("Significant for", label2),
      TRUE ~ "Insignificant"
    ),
    fortext = case_when(
      adj_p_value < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
      TRUE ~ "Insignificant"
    ),
    exclusive_to_one_group = (mean_intensity_g1 == 0 | mean_intensity_g2 == 0)
  )


insignificant <- comparison %>%
  dplyr::filter(significant == "Insignificant")

comparison <- comparison %>%
  separate(Canop, into = c("Canop", "av")) %>%
  left_join(sir_canop_pos_Kai %>% select(Canop, name)) %>%
  #left_join(feature_analysis %>% select(Canop, name)) %>%
  unique()

# Generate volcano plot
plot <- ggplot(comparison, aes(x = log2FoldChange, y = -log10(p_value))) +
  geom_point(aes(color = significant, shape = exclusive_to_one_group), alpha = 0.8, size = 4) +
  #geom_point(data = insignificant, shape = 21, fill = "lightgray", color = "gray", size = 5, stroke = 0.5) +
  geom_text_repel(aes(label = ifelse(fortext=="Significant", as.character(name), "")),
                  max.overlaps = 20, size = 1) +
  #scale_color_viridis(option = "H", direction = 1) +
  scale_color_manual(values = c("lightgray", "#58b2e9", "#e79e02")) + 
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +  # Different shape for exclusive features
  labs(
    #title = paste(label1, "vs", label2),
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Retention time",
    shape = "Exclusive"
  ) +
  Meine_theme
plot


#---------Analyzing individual compounds from important compound classes------------

fppath <- paste(admin, "AcidSinkhole/MSmine/FPtables", sep = "/")
setwd(fppath)
files <- dir(fppath, pattern = "20241106_canopus_can")

#canoplist <- highsinkholehalogen$Canop

canoplist <- c("Canop35", "Canop1028", "Canop1517", "Canop483")

canop_compounds <- fread(files[1]) %>%
  select(featureId, any_of(canoplist)) %>%
  gather(key = "Canop", value = "Canop_value", -featureId) %>%
  dplyr::filter(Canop_value == 1) %>%
  left_join(feature_analysis %>% mutate(featureId = as.numeric(featureId)), by = "featureId") %>%
  group_by(Canop, FileName) %>%
  mutate(intsum = sum(logInt)) %>%
  ungroup()

countclassperfeature <- canop_compounds %>%
  left_join(sir_canop_pos_Kai %>% mutate(name = absoluteIndex) %>% select (Canop, name)) %>%
  select(featureId, name) %>%
  distinct() %>%
  dplyr::group_by(featureId) %>%
  dplyr::summarise(
    count = dplyr::n(),
    classes = paste(name, collapse = ", ")
  ) %>%
  dplyr::ungroup() %>%
  select(featureId, classes) %>%
  unique()

classgroup <- countclassperfeature %>%
  select(classes) %>%
  unique() %>%
  rownames_to_column("classgroup")


cluster_levels <- c("Surface", "Platform", "Shallow_Sinkhole",
                    "Deep_Sinkhole", "Acid_Lake", "Mesopelagic")

canop_compounds <- canop_compounds %>%
  left_join(countclassperfeature) %>%
  left_join(classgroup) %>%
  #dplyr::filter(classgroup %in% c(8, 6, 9, 11)) %>%
  dplyr::group_by(featureId) %>%
  mutate(maxvalue = (1/max(logInt))) %>%
  dplyr::ungroup() %>%
  mutate(cluster = factor(cluster, levels = cluster_levels)) %>%
  left_join(sir_canop_pos_Kai %>% select (Canop, name))

#boxplot for individual features

featurebox <- ggplot(canop_compounds, aes(x = factor(cluster), y = log10(logInt), color = factor(cluster))) +
  geom_boxplot(outlier.shape = NA, aes(color = factor(cluster), fill = factor(cluster)), alpha = 0.7) +  # Avoid double-plotting outliers
  geom_jitter(width = 0.05, alpha = 0.7, size = 0.5) +  # Add points
  facet_wrap(~ maxvalue + classes + featureId, scales = "free_y", ncol = 8) +
  labs(x = "pH Cluster", y = "log10 Intensity") +
  Meine_theme + 
  scale_color_manual(values = group_c) +
  scale_fill_manual(values = group_c)
featurebox

setwd(paste(admin, "AcidSinkhole/RESULTS/Figures", sep = "/"))
svg("feature_int_plots_fourhalogen.svg",width = 25, height = 50, pointsize = 6)
print(featurebox)
dev.off()

#boxplot for classes

featurebox <- ggplot(canop_compounds, aes(x = factor(cluster), y = log10(intsum), color = factor(cluster))) +
  geom_boxplot(outlier.shape = NA, aes(color = factor(cluster), fill = factor(cluster)), alpha = 0.7) +  # Avoid double-plotting outliers
  geom_jitter(width = 0.05, alpha = 0.7, size = 0.5) +  # Add points
  facet_wrap(~ Canop, scales = "free_y", nrow = 1) +
  labs(x = "pH Cluster", y = "log10 peak area") +
  Meine_theme + 
  scale_color_manual(values = group_c) +
  scale_fill_manual(values = group_c)
featurebox


#just checking how many peaks in each group
canop_sub <- canop_compounds[, ..canop_cols]

# Now filter rows where at least one of those columns has a 1
canop_filtered <- canop_compounds[rowSums(canop_sub == 1) > 0, ]

print <- canop_filtered[, .(name = featureId, halogen = "yes")]


#---------CANOPUS class levels network in sinkhole data-----------

sir_folder <- paste(admin, "AcidSinkhole/SIRIUS/20241106", sep = "/")
setwd(sir_folder)


#Rank 1 structures only  #compound_identifications, canopus_compound_summary
comp_ind_files <- dir(sir_folder, pattern = "canopus_compound_summary.tsv", recursive = TRUE)
comp_ind_files
#file <- comp_ind_files[7]

top1annotations <- tibble()
for(file in comp_ind_files){
  id_file <- fread(file) # %>%
    #mutate(ConfidenceScore = as.numeric(ConfidenceScore))
  
  top1annotations <- top1annotations %>%
    bind_rows(id_file) %>%
    unique()
}

sirius_class_annotation <- top1annotations %>%
  select(id, molecularFormula, `NPC#pathway`, `ClassyFire#superclass`, `ClassyFire#all classifications`)

classyfire_all <- sirius_class_annotation %>%
  select(id, `ClassyFire#all classifications`) %>%
  rename(classyfire = `ClassyFire#all classifications`) %>%
  separate_wider_delim(classyfire, delim = ";")


classyfire_all <- classyfire_all %>%
  pivot_longer(cols = starts_with("classyfire"),
               names_to = "classification_number",
               values_to = "classification",
               values_drop_na = TRUE) %>%
  select(classification) %>%
  unique() %>%
  mutate(name = str_trim(classification)) %>%
  left_join(sir_canop_pos_Kai)



library(igraph)
library(ggraph)
# Example data frame


setwd("D:/Onedrive/OneDrive - Tartu Ülikool/Programs")

df <- fread("canopus.tsv") %>%
  mutate(Canop = paste("Canop", absoluteIndex, sep = "")) %>%
  mutate(across(where(is.character), ~ ifelse(nchar(.) < 1, NA, .))) %>%
  #dplyr::filter(name != "Canop0") %>%
  dplyr::filter(id %in% c("CHEMONT:0000035", "CHEMONT:0002279", "CHEMONT:0002448", "CHEMONT:0000000", "CHEMONT:0001028", 
                             "CHEMONT:0002867", "CHEMONT:0000267", "CHEMONT:0000483", "CHEMONT:0000262", 
                             "CHEMONT:0003909", "CHEMONT:0000012", "CHEMONT:0001518")) %>%
  select(Canop, name, id, parentId) %>%
  left_join(canop_compounds %>% select(Canop, intsum, FileName, cluster, featureId, logInt) %>% unique())

  

df <- fread("classyfire_taxonomy_parentId.csv")
#write_delim(df, "classyfire_taxonomy_parentId.csv", delim = ";")

df <- df %>%
  # mutate(id = trimws(as.character(id)),
  #        parentId = trimws(as.character(parentId))) %>%
  select(id, parentId, name)


# Create edges only where parentId exists in id
edges1 <- df %>%
  filter(!is.na(parentId) & parentId %in% id) %>%
  select(from = parentId, to = id)

edges2 <- df %>%
  filter(!is.na(featureId)) %>%
  select(from = featureId, to = id) %>%
  unique()

edges <- rbind(edges1, edges2) %>%
  unique()



# Create graph for just Classyfire network
g <- graph_from_data_frame(d = edges1, vertices = df, directed = TRUE)


plot <- ggraph(g, layout = "tree") +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_label(aes(label = name), repel = TRUE) +
  theme_minimal()
plot

setwd(paste(admin, "AcidSinkhole/RESULTS/Figures", sep = "/"))
svg("classyfire_taxonomy.svg",width = 1000, height = 300, pointsize = 6)
print(plot)
dev.off()


#chat

# Lae andmed
df <- fread("classyfire_parent_level_connections_groups.csv") %>%
  select(id, parentId, name)

# Veendu, et id ja parentId on tekstina
df <- df %>%
  mutate(id = trimws(as.character(id)),
         parentId = trimws(as.character(parentId)))

# Loo servad
edges <- df %>%
  filter(!is.na(parentId) & parentId %in% id) %>%
  select(from = parentId, to = id)

# Loo graaf
g <- graph_from_data_frame(d = edges, vertices = df, directed = TRUE)

# Leia paigutus (layout)
layout <- layout_as_tree(g)

# Koosta andmed sõlmede jaoks
nodes <- data.frame(layout)
nodes$id <- V(g)$name
nodes$name <- V(g)$name
nodes$label <- df$name[match(nodes$id, df$id)]

# Koosta andmed servade jaoks
edge_list <- get.edgelist(g)
edge_df <- data.frame(
  x = layout[match(edge_list[,1], V(g)$name), 1],
  y = layout[match(edge_list[,1], V(g)$name), 2],
  xend = layout[match(edge_list[,2], V(g)$name), 1],
  yend = layout[match(edge_list[,2], V(g)$name), 2]
)

# Plotly graafik
fig <- plot_ly(type = 'scatter', mode = 'markers+text') %>%
  add_segments(data = edge_df,
               x = ~x, y = ~-y, xend = ~xend, yend = ~-yend,
               line = list(color = 'gray')) %>%
  add_trace(data = nodes,
            x = ~X1, y = ~-X2,
            text = ~label,
            hovertext = ~id,
            mode = 'markers+text',
            marker = list(size = 5, color = 'blue'),
            textposition = 'top center') %>%
  layout(title = "Keemiliste ühendite taksonoomiapuu (interaktiivne)",
         xaxis = list(showgrid = FALSE, zeroline = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE),
         hovermode = 'closest')

fig

#------------Feature + class network------------

# df has: name (not unique), id (unique), parentId (points to id)
setwd("C:/Users/b15157/OneDrive - Tartu Ülikool/Programs")

df <- fread("canopus.tsv") %>%
  mutate(Canop = paste("Canop", absoluteIndex, sep = "")) %>%
  mutate(across(where(is.character), ~ ifelse(nchar(.) < 1, NA, .))) %>%
  #dplyr::filter(name != "Canop0") %>%
  dplyr::filter(id %in% c("CHEMONT:0000035", "CHEMONT:0002279", "CHEMONT:0002448", "CHEMONT:0000000", "CHEMONT:0001028", 
                          "CHEMONT:0002867", "CHEMONT:0000267", "CHEMONT:0000483", "CHEMONT:0000262", 
                          "CHEMONT:0003909", "CHEMONT:0000012", "CHEMONT:0001518")) %>%
  select(Canop, name, id, parentId) %>%
  left_join(canop_compounds %>% select(Canop, intsum, FileName, cluster, featureId, logInt) %>% unique())

parentname <- df %>%
  select(id, name) %>%
  mutate(parentname = name, 
         parentId = id) %>%
  select(parentId, parentname)

sinkholeint <- df %>%
  dplyr::filter(cluster == "Acid_Lake") %>%
  select(featureId, cluster, logInt) %>%
  unique() %>%
  dplyr::group_by(featureId, cluster) %>%
  mutate(sinkholeint = log10(mean(logInt))) %>%
  dplyr::ungroup() %>%
  select(featureId, sinkholeint) %>%
  unique()

canopintlake <- df %>%
  dplyr::filter(cluster == "Acid_Lake") %>%
  select(Canop, cluster, intsum) %>%
  unique() %>%
  dplyr::group_by(Canop) %>%
  mutate(canopintlake = log10(mean(intsum))) %>%
  dplyr::ungroup() %>%
  select(Canop, canopintlake) %>%
  unique()

df <- df %>%
  left_join(parentname) %>%
  mutate(id = name,
         parentId = parentname) %>%
  left_join(sinkholeint) %>%
  mutate(intsum = log10(intsum),
         logInt = log10(logInt)) %>%
  left_join(canopintlake)


#new!!!!

# 1. Numeric → Numeric edges (if featureId points to another featureId)
edges_numeric <- df %>% 
  filter(!is.na(featureId) & !is.na(id) & grepl("^[0-9]+$", featureId)) %>% # & grepl("^[0-9]+$", parentId)
  transmute(from = as.character(featureId), to = as.character(id)) %>%
  unique()

# 2. Numeric → String edges (numeric feature → Canop class)
edges_to_canop <- df %>%
  filter(!is.na(id) & !is.na(parentId)) %>%
  transmute(from = as.character(id), to = as.character(parentId)) %>%
  unique() %>%
  dplyr::filter(from != "CHEMONT:0000000")

# 3. Combine all
edges <- bind_rows(edges_numeric, edges_to_canop) %>% distinct()


all_nodes <- unique(c(edges$from, edges$to))

vertices <- tibble(name = all_nodes) %>%
  mutate(
    label = name,
    type = ifelse(grepl("^[0-9]+$", label), "numeric", "string")
  )

g <- graph_from_data_frame(edges, vertices = vertices, directed = TRUE)


#for int values size point

# Non-string vertices: names are numeric featureIds (as character); use sinkholeint
numeric_sizes <- df %>%
  distinct(featureId, sinkholeint) %>%
  filter(!is.na(featureId)) %>%
  transmute(name = as.character(featureId), size_raw = sinkholeint)

# String vertices: names are class names; aggregate intsum per class
# (Use mean by default; switch to sum(intsum, na.rm=TRUE) if you prefer)
string_sizes <- df %>%
  distinct(name, canopintlake) %>%
  dplyr::filter(!is.na(name)) %>%
  transmute(name = as.character(name), size_raw = canopintlake)
  # filter(!is.na(name)) %>%
  # group_by(name) %>%
  # summarise(size_raw = sum(intsum, na.rm = TRUE), .groups = "drop")

# Combine and align to graph vertex order
size_map <- bind_rows(numeric_sizes, string_sizes) %>%
  distinct(name, .keep_all = TRUE)

vsize_raw <- size_map$size_raw[match(V(g)$name, size_map$name)]
vsize_raw[!is.finite(vsize_raw)] <- NA

# Replace missing sizes with the median of available sizes
if (all(is.na(vsize_raw))) {
  vsize_plot <- rep(12, length(vsize_raw))  # fallback if nothing available
} else {
  x <- vsize_raw
  x[is.na(x)] <- median(x, na.rm = TRUE)
  
  # Optional: compress heavy-tailed string-node sizes
  # (e.g., intsum can be large; sqrt keeps relative differences without huge circles)
  #x <- ifelse(is_string, sqrt(pmax(x, 0)), x)
  
  # Linearly map to [6, 24] points (adjust to taste)
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    vsize_plot <- rep(12, length(x))
  } else {
    vsize_plot <- 6 + (x - rng[1]) * (24 - 6) / (rng[2] - rng[1])
  }
}

# Assign to graph
V(g)$size <- vsize_plot


is_string <- V(g)$type == "string"


manual_colors <- c(
  "Organic compounds" = "#FF6633",
  "Lipids and lipid-like molecules" = "#FF9900",
  "Fatty acids and conjugates" = "#F0E442",
  "Halobenzenes" = "#CC0033",
  "Organohalogen compounds" = "#FF9900",
  "Halogenated fatty acids" = "#003399",
  "Alkyl iodides" = "#339900",
  "Organoiodides" = "#993399",
  "Benzene and substituted derivatives" = "#FFCC66",
  "Benzenoids" = "#FF9900",
  "Alkyl halides" = "#FFCC66",
  "Fatty Acyls" = "#FFCC66"
)

# Assign colors to string nodes
string_colors <- manual_colors[V(g)$name[is_string]]

# Set V(g)$color
V(g)$color[is_string] <- string_colors


# Shape strong points
V(g)$shape[is_string] <- "circle"

# Optional: default shape for others
V(g)$shape[!is_string] <- "circle"


# propagate colors to numeric nodes
for (v in V(g)[!is_string]) {
  reachable <- subcomponent(g, v, mode = "out")
  reachable_strings <- reachable[is_string[reachable]]
  if (length(reachable_strings) > 0) {
    V(g)$color[v] <- V(g)$color[reachable_strings[1]]
  } else {
    V(g)$color[v] <- "#CCCCCC"
  }
}
#for all gray
V(g)$color[!is_string] <- "#CCCCCC"

to_vids <- ends(g, E(g), names = FALSE)[, 2]

for (i in seq_along(to_vids)) {
  target <- to_vids[i]
  reachable <- subcomponent(g, target, mode = "out")
  s <- reachable[is_string[reachable]]
  if (length(s) > 0) {
    E(g)$color[i] <- V(g)$color[s[1]]
  } else {
    E(g)$color[i] <- "#BBBBBB"
  }
}

set.seed(85)
plot(
  g,
  layout = layout_with_fr(g),
  vertex.color = V(g)$color,
  edge.color = E(g)$color,
  vertex.label = V(g)$label,
  vertex.label.cex = 1,
  vertex.size = V(g)$size,   # <-- sizes from sinkholeint / i
  vertex.size = 12,
  edge.arrow.size = 0.2,
  edge.width = 3,
  main = "Network colored by Canop class propagation"
)

setwd(paste(admin, "AcidSinkhole/RESULTS/Figures", sep = "/"))
ggsave(
  filename = "canopus_network_nodes.svg",
  plot = last_plot(),
  device = svglite::svglite,
  width = 12,
  height = 10
)


#---------Iodine and Bromine compounds ALL------------


rank1_halides <- sirius_formula_annotation %>%
  mutate(id = featureId) %>%
  mutate(halide = grepl("I|Br", molecularFormula)) %>%
  select(id, halide)


rankAny_halides <- sirius_formula_annotation_all %>%
  mutate(id = featureId) %>%
  mutate(halide = grepl("I|Br", molecularFormula)) %>%
  select(id, halide) %>%
  dplyr::filter(halide == TRUE) %>%
  unique()
