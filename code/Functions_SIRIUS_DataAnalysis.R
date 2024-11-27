#Functions for SIRIUS data analysis


#selected_structures <- SiriusStructures(folderdata, subfoldercommand = "SIRIUScommand5", subfoldercommand_exl = "allpeaks")


SiriusStructures <- function(folder, subfoldercommand, subfoldercommand_exl){
  
  setwd(folder)
  subfolders_structure <- dir(folder, all.files = TRUE, recursive = TRUE, pattern = "structure_candidates")
  subfolders_structure <- subfolders_structure[grepl(subfoldercommand, subfolders_structure)]
  #subfolders_structure <- subfolders_structure[!grepl(subfoldercommand_exl, subfolders_structure)]
  
  
  structuretable <- tibble()
  for(subfolder in subfolders_structure){
    subfolder_name_split <- str_split(subfolder, "/")
    subfolder_foldername <- str_split(subfolder_name_split[[1]][2], "_")
    
    #sample_name <- subfolder_foldername[[1]][1]
    foldernr <- subfolder_foldername[[1]][1]
    #aq_mode <- subfolder_name_split[[1]][6]  #[[1]][1] for one subfolder at a time
    id <- subfolder_foldername[[1]][length(subfolder_foldername[[1]])]
    
    structures <- read_delim(subfolder, delim = "\t", show_col_types = FALSE) %>%
      select(molecularFormula, formulaRank, structureRankPerFormula, adduct, `CSI:FingerIDScore`, 
             ConfidenceScore, smiles, InChIkey2D, InChI) %>%
      mutate(ConfidenceScore = as.numeric(ConfidenceScore)) %>%
      mutate(foldernr = foldernr) %>%
      mutate(id = id)
    
    if(nrow(structures) != 0){
      structuretable <- structuretable %>%
        bind_rows(structures)
    }
  }
  
  selected_structures <- structuretable
  # filter(formulaRank <= 3) %>%
  # filter(rank <= 10)
  
  # sel_H <- selected_structures %>%
  #   filter(adduct == "[M + H]+") %>%
  #   mutate(adduct = "[M+H]+")
  # sel_Na <- selected_structures %>%
  #   filter(adduct == "[M + Na]+") %>%
  #   mutate(adduct = "[M+Na]+")
  # sel_Hneg <- selected_structures %>%
  #   filter(adduct == "[M - H]-") %>%
  #   mutate(adduct = "[M-H]-")
  # 
  # selected_structures <- rbind(sel_H, sel_Na, sel_Hneg) #to combine with chemical info
}
