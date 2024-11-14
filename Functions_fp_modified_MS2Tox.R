#fingerprint table with all fingerprints from SIRIUS
library(dplyr)
library(rlist)
library(Rdisop)
library(tibble)
library(readr)
library(magrittr)
library(xgboost)
library(tidyselect)
library(stringr)
library(rcdklibs)
library(rcdk)

#------------------------Functions NEW with unzipping--------------------------------

# folderwithSIRIUSfiles <- "C:/Users/pille/Desktop/SIRIUS/sinkhole/20240628/test"
# setwd(folderwithSIRIUSfiles)


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
  
  
  
  if(grepl("-.fpt", all_fp$FPformula[1])) {
    colnames(all_fp) <- c("file", "FPformula", fp_names_neg)
  }
  
  if(grepl("+.fpt", all_fp$FPformula[1])) {
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
  
  if(grepl("-.fpt", all_fp$FPformula[1])) {
    colnames(all_fp) <- c("file", "FPformula", fp_names_neg)
  }
  
  if(grepl("+.fpt", all_fp$FPformula[1])) {
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





#------------------------Functions OLD with unzipping--------------------------------
#' FpTableForPredictionsMetabol <- function(folderwithSIRIUSfiles){
#'   
#'   subfolder <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = ".fpt")
#'   subfolder_score <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = ".info")
#'   
#'   fp_names_pos <- paste("Un", read_delim(paste(folderwithSIRIUSfiles,"/csi_fingerid.tsv", sep = ""), delim = "\t")$absoluteIndex, sep = "")
#'   fp_names_neg <- paste("Un", read_delim(paste(folderwithSIRIUSfiles,"/csi_fingerid_neg.tsv", sep = "" ), delim = "\t")$absoluteIndex, sep = "")
#'   fp_names_common <- as.data.frame(fp_names_pos)  %>%
#'     mutate(fp_names_neg = fp_names_pos) %>%
#'     inner_join(as.data.frame(fp_names_neg))
#'   fp_names_common <- as.vector(fp_names_common$fp_names_pos)
#'   
#'   fb_table <- FingerPrintTable(subfolder, fp_names_pos, fp_names_neg, fp_names_common, folderwithSIRIUSfiles)
#'   
#'   ranks <- SiriusScoreRank1(subfolder_score, folderwithSIRIUSfiles)
#'   
#'   final_table_mass <- ranks %>%
#'     inner_join(fb_table, by = c("id", "foldernumber", "predion"))
#'   
#' 
#'   return(final_table_mass)
#' }
#' 
#' 
#' FingerPrintTable <- function(subfolder, fp_names_pos, fp_names_neg, fp_names_common, folderwithSIRIUSfiles){
#'   
#'   subfolder_pos <- c()
#'   subfolder_neg <- c()
#'   for(subfold in subfolder){
#'     if ((grepl("]2+", subfold, fixed=TRUE) | grepl("]+", subfold, fixed=TRUE)) & grepl("/fingerprint", subfold, fixed=TRUE)){
#'       subfolder_pos <- subfolder_pos %>%
#'         list.append(subfold)
#'     }
#'     
#'     if ((grepl("]2-", subfold, fixed=TRUE) | grepl("]-", subfold, fixed=TRUE)) & grepl("/fingerprint", subfold, fixed=TRUE)){
#'       subfolder_neg <- subfolder_neg %>%
#'         list.append(subfold)
#'     }
#'   }
#'   
#'   fp_pos <- FingerPrintTablePOS(subfolder_pos, folderwithSIRIUSfiles)
#'   
#'   if(nrow(fp_pos) != 0){
#'     colnames(fp_pos) <- c(fp_names_pos, "predion", "id", "foldernumber", "predform")
#'   }
#'   
#'   fp_neg <- FingerPrintTableNEG(subfolder_neg, folderwithSIRIUSfiles)
#'   
#'   if(nrow(fp_neg) != 0){
#'     colnames(fp_neg) <- c(fp_names_neg, "predion", "id", "foldernumber", "predform")
#'   }
#'   
#'   if (!is.null(subfolder_neg) & !is.null(subfolder_pos)){
#'     
#'     fp_all <- fp_pos %>%
#'       select(id, foldernumber, predform, predion, all_of(fp_names_common)) %>%
#'       bind_rows(fp_neg %>% select(id, foldernumber, predform, predion, all_of(fp_names_common))) %>%
#'       unique()
#'   }
#'   
#'   
#'   if (is.null(subfolder_neg) & !is.null(subfolder_pos)){
#'     
#'     fp_all <- fp_pos %>%
#'       select(id, foldernumber, predform, predion, all_of(fp_names_common)) %>%
#'       unique()
#'   }
#'   
#'   if (!is.null(subfolder_neg) & is.null(subfolder_pos)){
#'     
#'     fp_all <- fp_neg %>%
#'       select(id, foldernumber, predform, predion, all_of(fp_names_common)) %>%
#'       unique()
#'   }
#'   
#'   if (is.null(subfolder_neg) & is.null(subfolder_pos)){
#'     print("NO data found in subfolders")
#'   }
#'   
#'   return(fp_all)
#' }
#' 
#' #' @export
#' FingerPrintTablePOS <- function(subfolder, folderwithSIRIUSfiles){
#'   fingerprint_data <- tibble()
#'   for(direct in subfolder){
#'     # subfolder with data must be on a form where id is on second place after nr_ (0_Einc270001_Einc270001)
#'     file_name <- str_split(direct, "_")
#'     comp_name <- str_split(direct, "/")
#'     
#'     sir_fold <- file_name[[1]][1]
#'     #id_this <- file_name[[1]][2] # if id in some other position, this can be changed
#'     id_this <- comp_name[[1]][1]
#'     pred_ion <- comp_name[[1]][3]
#'     
#'     filedata <- read_delim(paste(folderwithSIRIUSfiles, direct, sep = "/"), delim = " ", col_names = FALSE)
#'     filedata <- as_tibble(t(filedata))
#'     filedata <- filedata %>%
#'       mutate(predion = pred_ion) %>%
#'       mutate(predion = sub("\\..*", "", predion)) %>%
#'       mutate(id = id_this) %>%
#'       mutate(sir_fol_nr = sir_fold) %>%
#'       mutate(predform = sub("\\_.*", "", predion))
#'     fingerprint_data <- fingerprint_data %>%
#'       bind_rows(filedata)
#'   }
#'   return(fingerprint_data)
#' }
#' 
#' #' @export
#' FingerPrintTableNEG <- function(subfolder, folderwithSIRIUSfiles){
#'   fingerprint_data <- tibble()
#'   for(direct in subfolder){
#'     # subfolder with data must be on a form where id is on second place after nr_ (0_Einc270001_Einc270001)
#'     file_name <- str_split(direct, "_")
#'     comp_name <- str_split(direct, "/")
#'     
#'     sir_fold <- file_name[[1]][1]
#'     id_this <- comp_name[[1]][1]
#'     pred_ion <- comp_name[[1]][3]
#'     
#'     filedata <- read_delim(paste(folderwithSIRIUSfiles, direct, sep = "/"), delim = " ", col_names = FALSE)
#'     filedata <- as_tibble(t(filedata))
#'     filedata <- filedata %>%
#'       mutate(predion = pred_ion) %>%
#'       mutate(predion = sub("\\..*", "", predion)) %>%
#'       mutate(id = id_this) %>%
#'       mutate(sir_fol_nr = sir_fold) %>%
#'       mutate(predform = sub("\\_.*", "", predion))
#'     fingerprint_data <- fingerprint_data %>%
#'       bind_rows(filedata)
#'     
#'   }
#'   return(fingerprint_data)
#' }
#' 
#' #' @export
#' SiriusScoreRank1 <- function(subfolder_score, folderwithSIRIUSfiles){
#'   scores_table <- tibble()
#'   scores_table <- scores_table %>%
#'     mutate(id = "id") %>%
#'     mutate(foldernumber = "sirfolnr") %>%
#'     mutate(siriusscore = "siriusscore")
#'   
#'   for (filename in subfolder_score) {
#'     if (grepl("/scores", filename, fixed=TRUE)){
#'       
#'       file_name_score <- str_split(filename, "_")
#'       comp_name_score <- str_split(filename, "/")
#'       filename2_exeption <- 
#'       
#'       foldernumber <- file_name_score[[1]][1]
#'       id <- comp_name_score[[1]][1] #if id is not in second place, [2] must be changed
#'       pred_st <- comp_name_score[[1]][3]
#'       
#'       fileConnection <- file(paste(folderwithSIRIUSfiles, filename, sep = "/"))
#'       record <- readLines(fileConnection)
#'       close(fileConnection)
#'       
#'       SiriusScore <- substring(grep('sirius.scores.SiriusScore', record, value = TRUE, fixed = TRUE),27)
#'       
#'       filedata <- data.frame(id , foldernumber)
#'       
#'       filedata <- filedata %>%
#'         mutate(predion = pred_st) %>%
#'         mutate(predion = sub("\\..*", "", predion)) %>%
#'         mutate(siriusscore = SiriusScore)
#'       
#'       scores_table <- scores_table%>%
#'         bind_rows(filedata)
#'     }
#'   }
#'   data_scores <- scores_table %>%
#'     unique() %>%
#'     select(-predion) %>% #if two compound have same siriusscore they get same rank
#'     unique() %>%
#'     mutate(siriusscore = as.numeric(siriusscore)) %>%
#'     group_by(id, foldernumber) %>%
#'     arrange(desc(siriusscore)) %>%
#'     mutate(rank = row_number()) %>%
#'     ungroup() %>%
#'     left_join(scores_table %>% unique() %>%  mutate(siriusscore = as.numeric(siriusscore)), 
#'               by = c("id", "foldernumber", "siriusscore"),
#'               multiple = "all") %>%
#'     #dplyr::filter(rank == 1) %>%
#'     #select(-rank, -siriusscore) %>%
#'     unique()
#'   
#'   return(data_scores)
#' }
#' 
#' UnZip_SIRIUS5_modified <- function(folderwithSIRIUSfiles) {
#'   
#'   subfolder_fp_zip <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = "fingerprints")
#'   
#'   for (zipF in subfolder_fp_zip){
#'     outfolder <- str_split(zipF, "/")[[1]][1]
#'     outDir <- paste(folderwithSIRIUSfiles, outfolder, "fingerprints1",
#'                     sep = "/")
#'     zipFile <- paste(folderwithSIRIUSfiles, zipF, sep = "/")
#'     unzip(zipFile, exdir = outDir)
#'   }
#'   
#'   
#'   subfolder_scores_zip <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = "scores")
#'   
#'   for (zipF in subfolder_scores_zip) {
#'     outfolder <- str_split(zipF, "/")[[1]][1]
#'     outDir <- paste(folderwithSIRIUSfiles, outfolder, "scores1",
#'                     sep = "/")
#'     zipFile <- paste(folderwithSIRIUSfiles, zipF, sep = "/")
#'     unzip(zipFile, exdir = outDir)
#'   }
#'   
#'   
#'   subfolder_canopus_zip <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = "canopus")
#'   subfolder_canopus_zip <- subfolder_canopus_zip[!grepl("_npc", subfolder_canopus_zip)]
#'   subfolder_canopus_zip <- subfolder_canopus_zip[!grepl(".tsv", subfolder_canopus_zip)]
#'   
#'   for (zipF in subfolder_canopus_zip) {
#'     outfolder <- str_split(zipF, "/")[[1]][1]
#'     outDir <- paste(folderwithSIRIUSfiles, outfolder, "canopus1",
#'                     sep = "/")
#'     zipFile <- paste(folderwithSIRIUSfiles, zipF, sep = "/")
#'     unzip(zipFile, exdir = outDir)
#'   }
#'   
#'   
#'   subfolder_canopusnpc_zip <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = "canopus_npc")
#'   subfolder_canopusnpc_zip <- subfolder_canopusnpc_zip[!grepl(".tsv", subfolder_canopusnpc_zip)]
#'   
#'   for (zipF in subfolder_canopusnpc_zip) {
#'     outfolder <- str_split(zipF, "/")[[1]][1]
#'     outDir <- paste(folderwithSIRIUSfiles, outfolder, "canopus_npc1",
#'                     sep = "/")
#'     zipFile <- paste(folderwithSIRIUSfiles, zipF, sep = "/")
#'     unzip(zipFile, exdir = outDir)
#'   }
#'   
#' }
#' 
#' #Canopus
#' 
#' CANOPUSTableForPredictionsMetabol <- function(folderwithSIRIUSfiles){
#'   
#'   subfolder <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = ".fpt")
#'   subfolder_score <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = ".info")
#'   
#'   fp_names_pos <- paste("Canop", read_delim(paste(folderwithSIRIUSfiles,"/canopus.tsv", sep = ""), delim = "\t")$absoluteIndex, sep = "")
#'   fp_names_neg <- paste("Canop", read_delim(paste(folderwithSIRIUSfiles,"/canopus_neg.tsv", sep = "" ), delim = "\t")$absoluteIndex, sep = "")
#'   fp_names_common <- as.data.frame(fp_names_pos)  %>%
#'     mutate(fp_names_neg = fp_names_pos) %>%
#'     inner_join(as.data.frame(fp_names_neg))
#'   fp_names_common <- as.vector(fp_names_common$fp_names_pos)
#'   
#'   fp_table <- CanopusTable(subfolder, fp_names_pos, fp_names_neg, fp_names_common, folderwithSIRIUSfiles)
#'   
#'   ranks <- SiriusScoreRank1(subfolder_score, folderwithSIRIUSfiles)
#'   
#'   final_table_mass <- ranks %>%
#'     inner_join(fp_table, by = c("id", "foldernumber", "predion"))
#'   
#'   
#'   return(final_table_mass)
#' }
#' 
#' 
#' CanopusTable <- function(subfolder, fp_names_pos, fp_names_neg, fp_names_common, folderwithSIRIUSfiles){
#'   
#'   subfolder_pos <- c()
#'   subfolder_neg <- c()
#'   for(subfold in subfolder){
#'     if ((grepl("]2+", subfold, fixed=TRUE) | grepl("]+", subfold, fixed=TRUE)) & grepl("/canopus1", subfold, fixed=TRUE)){
#'       subfolder_pos <- subfolder_pos %>%
#'         list.append(subfold)
#'     }
#'     
#'     if ((grepl("]2-", subfold, fixed=TRUE) | grepl("]-", subfold, fixed=TRUE)) & grepl("/canopus1", subfold, fixed=TRUE)){
#'       subfolder_neg <- subfolder_neg %>%
#'         list.append(subfold)
#'     }
#'   }
#'   
#'   fp_pos <- FingerPrintTablePOS(subfolder_pos, folderwithSIRIUSfiles)
#'   
#'   if(nrow(fp_pos) != 0){
#'     colnames(fp_pos) <- c(fp_names_pos, "predion", "id", "foldernumber", "predform")
#'   }
#'   
#'   fp_neg <- FingerPrintTableNEG(subfolder_neg, folderwithSIRIUSfiles)
#'   
#'   if(nrow(fp_neg) != 0){
#'     colnames(fp_neg) <- c(fp_names_neg, "predion", "id", "foldernumber", "predform")
#'   }
#'   
#'   if (!is.null(subfolder_neg) & !is.null(subfolder_pos)){
#'     
#'     fp_all <- fp_pos %>%
#'       select(id, foldernumber, predform, predion, all_of(fp_names_common)) %>%
#'       bind_rows(fp_neg %>% select(id, foldernumber, predform, predion, all_of(fp_names_common))) %>%
#'       unique()
#'   }
#'   
#'   
#'   if (is.null(subfolder_neg) & !is.null(subfolder_pos)){
#'     
#'     fp_all <- fp_pos %>%
#'       select(id, foldernumber, predform, predion, all_of(fp_names_common)) %>%
#'       unique()
#'   }
#'   
#'   if (!is.null(subfolder_neg) & is.null(subfolder_pos)){
#'     
#'     fp_all <- fp_neg %>%
#'       select(id, foldernumber, predform, predion, all_of(fp_names_common)) %>%
#'       unique()
#'   }
#'   
#'   if (is.null(subfolder_neg) & is.null(subfolder_pos)){
#'     print("NO data found in subfolders")
#'   }
#'   
#'   return(fp_all)
#' }
#' 
