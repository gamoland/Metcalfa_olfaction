library(tidyverse)
library(openxlsx)
library(stringi)
source("Scripts/Function_normalization_v1.R")
source("Scripts/Function_read_in_EAD_v1.R")
source("Scripts/Function_read_in_MS_v1.R")
source("Scripts/Heatmap_functions_v2.R")
source("Scripts/PCA_function_v1.R")

Tagetes <- rEADin(File_Dir_In = "../OlfactoMix_platform/Data/Metcalfa/metcalfa_tagetes_SebféleTemplate.xlsx") %>% 
mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "random")
  
Ailanthus <- rEADin(File_Dir_In = "../OlfactoMix_platform/Data/Metcalfa/metcalfa_ailanthus_SebféleTemplate.xlsx") %>% 
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "random")

Normalized_Data_t <- NORmalization(Read_Data_In = Tagetes) 

Top_3_data <- Normalized_Data_t %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "random") %>%
  mutate(Randoms = stri_extract_first_regex(str = Compound.Name, pattern = "[:alpha:]+")) %>%
  filter(Randoms != "unknown") %>%
  arrange(-Normalized_response) %>%
  group_by(Compound.Name, Sex) %>%
  mutate(Nrow = row_number()) %>% 
  ungroup() %>% 
  filter(Nrow <= 3) %>%
  mutate(Replicate_number = Nrow,
         Compound.Name = as.factor(Compound.Name))
  
   slice_max(order_by = Normalized_response, n = 3) %>%
   ungroup()

Order <-Top_3_data %>%
  group_by(Compound.Name) %>%
  summarise(Sum_resp= sum(Normalized_response)) %>%
  ungroup() %>%
  arrange(-Sum_resp)


mutate(Compound.Name =fct_inorder(Compound.Name))
Top_3_data$Compound.Name <- Top_3_data$Compound.Name %>% fct_relevel(levels(Order$Compound.Name))

Chosing_Heatmap(Top_3_data, Heatmap_type_In = "Ave" )
PCA(Top_3_data, Group = "Sex")

