library(tidyverse)
library(openxlsx)

# function to read in MS data
rMSin <- function(File_Dir_In, Volatilome_name = "Volatilome"){ # frunction rMSin takes 2 arguments
                                                                # file_Dir_In and volatilome_name
  
  # File_Dir_In <- "Data/Template for annotating_Metcalfa.xlsx"
  # Volatilome_name = "Volatilome"
  
Volatilome <- read.xlsx(File_Dir_In, # reading data from excel 
                        sheet = Volatilome_name, 
                        startRow = 4) %>% 
  filter(!is.na(No.peak)) # removing empty entries. reading data without empty entries

Volatilome_sample <- read.xlsx(xlsxFile = File_Dir_In, # reading data from excel 
                        sheet = Volatilome_name, 
                        rows = c(2), cols = c(3:4), colNames = F) %>% 
  spread(key = X1, value =  X2) %>% 
  mutate(Common = "a")

New_tibble_joint <- Volatilome %>% 
  mutate(Common = "a") %>% 
  #mutate(CAS = as.character(CAS)) %>% 
  left_join(Volatilome_sample) %>% 
  dplyr::select(-Common)

New_tibble_joint %>% return() # returning  object volatilome
}
  
