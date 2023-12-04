source("Scripts/Function_read_in_EAD_v1.R")
source("Scripts/Function_normalization_v1.R")
source("Scripts/Metcalfa/EAD_readin_withNAs.R")
source("Scripts/Function_read_in_MS_v1.R")


#READING DATA WITHOUT NA ANSWERS
  Data_readin_AIL <- rEADin(File_Dir_In = "Data/metcalfa_ailanthus_SebféleTemplate.xlsx") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "random") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "unknown")
Data_readin_ARIST <- rEADin(File_Dir_In = "Data/metcalfa_aristolochia_SebféleTemplate.xlsx") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "random") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "unknown")
Data_readin_TAG <- rEADin(File_Dir_In = "Data/metcalfa_tagetes_Seb.xlsx") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "random") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "unknown")


#READING DATA WITH NA ANSWERS
NAs_EAD_readIN_AIL <- rEADinNAs(File_Dir_In = "Data/metcalfa_ailanthus_SebféleTemplate.xlsx") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "random") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "unknown") %>% 
  mutate(MV_fisher = ifelse(is.na(MV), 0, 1))
NAs_EAD_readIN_ARIST <- rEADinNAs(File_Dir_In = "Data/metcalfa_aristolochia_SebféleTemplate.xlsx") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "random") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "unknown") %>% 
  mutate(MV_fisher = ifelse(is.na(MV), 0, 1))
NAs_EAD_readIN_TAG <- rEADinNAs(File_Dir_In = "Data/metcalfa_tagetes_Seb.xlsx") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "random") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "unknown") %>% 
  mutate(MV_fisher = ifelse(is.na(MV), 0, 1))


#NORMALIZING DATA
Normalized_Data_Metc_AIL <- NORmalization(Data_readin_AIL) %>%
  arrange(-Normalized_response) %>%
  group_by(Compound.Name, Sex) %>%
  mutate(Nrow = row_number()) %>% ungroup() %>%
  mutate(Replicate_number = Nrow,
         Compound.Name = as.factor(Compound.Name))
Normalized_Data_Metc_ARIST <- NORmalization(Data_readin_ARIST) %>%
  arrange(-Normalized_response) %>%
  group_by(Compound.Name, Sex) %>%
  mutate(Nrow = row_number()) %>% ungroup() %>%
  mutate(Replicate_number = Nrow,
         Compound.Name = as.factor(Compound.Name))
Normalized_Data_Metc_TAG <- NORmalization(Data_readin_TAG) %>%
  arrange(-Normalized_response) %>%
  group_by(Compound.Name, Sex) %>%
  mutate(Nrow = row_number()) %>% ungroup() %>%
  mutate(Replicate_number = Nrow,
         Compound.Name = as.factor(Compound.Name))

All_metcalfa_data_norm <-
  bind_rows(Normalized_Data_Metc_AIL,
            Normalized_Data_Metc_ARIST,
            Normalized_Data_Metc_TAG) %>%
  mutate(Compound.Name = tolower(Compound.Name))


#NORMALIZING NA_DATA
NAs_Normalized_Data_Metc_AIL <- NORmalization(NAs_EAD_readIN_AIL) %>%
  arrange(-Normalized_response) %>%
  group_by(Compound.Name, Sex) %>%
  mutate(Nrow = row_number()) %>% ungroup() %>%
  mutate(Replicate_number = Nrow,
         Compound.Name = as.factor(Compound.Name))
NAs_Normalized_Data_Metc_ARIST <- NORmalization(NAs_EAD_readIN_ARIST) %>%
  arrange(-Normalized_response) %>%
  group_by(Compound.Name, Sex) %>%
  mutate(Nrow = row_number()) %>% ungroup() %>%
  mutate(Replicate_number = Nrow,
         Compound.Name = as.factor(Compound.Name))
NAs_Normalized_Data_Metc_TAG <- NORmalization(NAs_EAD_readIN_TAG) %>%
  arrange(-Normalized_response) %>%
  group_by(Compound.Name, Sex) %>%
  mutate(Nrow = row_number()) %>% ungroup() %>%
  mutate(Replicate_number = Nrow,
         Compound.Name = as.factor(Compound.Name))

NAs_All_metcalfa_data_norm <-
  bind_rows(NAs_Normalized_Data_Metc_AIL,
            NAs_Normalized_Data_Metc_ARIST,
            NAs_Normalized_Data_Metc_TAG) %>%
  mutate(Compound.Name = tolower(Compound.Name))


#READIN MS DATA
MS_Tagetes <- rMSin("Data/metcalfa_tagetes_Seb.xlsx") %>% 
  filter(!is.na(Compound.Name)) %>% 
  mutate(CAS = as.character(CAS),
         Chemical.group = as.character(Chemical.group)) %>% 
  dplyr::select(Compound.Name, CAS, Area.from.MS,"Sample descriptor", "Published.kovats(closest)")
MS_Ailanthus <-rMSin("Data/metcalfa_ailanthus_SebféleTemplate.xlsx") %>% 
  filter(!is.na(Compound.Name)) %>% 
  mutate(CAS = as.character(CAS),
         Chemical.group = as.character(Chemical.group))%>% 
  dplyr::select(Compound.Name, CAS, Area.from.MS,"Sample descriptor", "Published.kovats(closest)")
MS_Aristolochia <- rMSin("Data//metcalfa_aristolochia_SebféleTemplate.xlsx") %>% 
  filter(!is.na(Compound.Name)) %>% 
  mutate(CAS = as.character(CAS),
         Chemical.group = as.character(Chemical.group))%>% 
  dplyr::select(Compound.Name, CAS, Area.from.MS,"Sample descriptor", "Published.kovats(closest)")

All_MS <- MS_Tagetes %>%
  full_join(MS_Ailanthus) %>%
  full_join(MS_Aristolochia)

