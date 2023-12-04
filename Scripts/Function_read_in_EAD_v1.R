library(tidyverse)
library(openxlsx)
#library(xlsx) Should be removed

rEADin <- function(File_Dir_In, Volatilome_name = "Volatilome"){ #function name, input dataset
  
  File_dir <- File_Dir_In #name for the input dataset
  
  SheetNames <- getSheetNames(File_dir) #get the sheets' names
  # Var_1 <- sym(CAS_code) 
  Read_In_sheets <- function(Sheet_Names_In  , CAS_code = Var_1) { #where do i give the input sheet?
    
    Selected_sheet <- Sheet_Names_In #
    
    M_females <- read.xlsx(xlsxFile = File_dir, sheet = Selected_sheet, startRow = 7,check.names = T, skipEmptyRows = T) # read in data from the chosen sheet
    Key <- read.xlsx(xlsxFile = File_dir, sheet = Selected_sheet, startRow = 1, rows = 1:5, cols = c(4,6), colNames = FALSE) #read in the "key infos" frrom the same sheet
    
    Spread_Key <- Key %>% 
      spread(key = X1, value = X2)  #ungathering
    
    Lenght_of_data <- length(M_females) #counts the coulumns in tha data, how can i count the rows? RowNumber or n to count the observation
    First_data <- M_females %>% dplyr::select(1:17) #getting the "constants"
    Replicates_number <- (Lenght_of_data - length(First_data))/3 # counting how many replicates the dataset has
    
    New_tibble <- tibble() # new tibble
    i_a <- 1 #where should the index start for the loop, it's always 1
    for (i_a in 1:Replicates_number) { #how many loop should be done
      i_b <- 18:20 + 3*(i_a-1) #where should we find the next replicates
      CAS_code <- sym(colnames(M_females)[4])
      Loop_data <- M_females %>% 
        dplyr::select(c(1:17, i_b)) %>% # select the given coulumns. is i_b a list?
        dplyr::rename(Rep = 18, Time = 19, MV = 20) %>% 
       # mutate(CAS = as.character(CAS)) %>% 
        mutate(Replicate_number = i_a) %>% 
        mutate(Replicate_number = i_a, Rep = as.character(Rep),
               Time = as.numeric(Time),
               MV = as.numeric(MV),
               Comments = as.character(Comments),
               GC.EAD.KOVATS = as.numeric(GC.EAD.KOVATS),
               Chemical.group = as.character(Chemical.group),
               Compound.Name = as.character(Compound.Name),
               Closest.published.kovats = as.numeric(Closest.published.kovats),
               No.peak = as.numeric(No.peak),
               GC.MS.RT = as.numeric(GC.MS.RT),
               GC.EAD.RT = as.numeric(GC.EAD.RT),
               CAS.INCHI = if_else(is.na(as.character(!!CAS_code)), as.character(Compound.Name), 
                                   as.character(!!CAS_code))) %>%  #new variable to count the replicate
        dplyr::select(-n.1,-N.1)
      
      New_tibble <- New_tibble %>% bind_rows(Loop_data) #binding the tibble with the the loop dataset
      i_a <- i_a +1 # end of the loop
    }
    
    New_tibble <- New_tibble %>%  
      filter(!is.na(MV)) %>% #filtering out the empty MV rows
      mutate(Sheet = Selected_sheet) # to be able to join
    
    Full_data_sheet <- New_tibble %>% mutate(A = "1") %>% 
      left_join(Spread_Key %>% mutate(A = "1")) #joining the two dataset
    
    Full_data_sheet %>% return()
  }
  
  Sheets_names <- SheetNames %>% 
    enframe() %>% #make it tablellike
    filter(!value %in% c(Volatilome_name, "Kovats n-alkanes", "Volatilome_Key")) 
  
  Names <- Sheets_names$value 
  i_a <- 1
  New_tibble_for_everything <- tibble()
  for(i_a in 1:length(Names)) {
    Resulting_from_function <- Read_In_sheets(Names[i_a]) #this is where i give the sheet with the loop? so upto this point, the function was just an emty shell, and now it can run?
    New_tibble_for_everything <- New_tibble_for_everything %>% 
      bind_rows(Resulting_from_function)
    i_a <- i_a + 1
  }
  
  Volatilome <- read.xlsx(xlsxFile = File_dir, # reading data from excel 
                          sheet = Volatilome_name, 
                          rows = c(2), cols = c(3:4), colNames = F)
  
  New_tibble_joint <-  switch(if_else(length(Volatilome) > 1, "A", "B"),
                              "A" = New_tibble_for_everything %>% 
                                mutate(Common = "a") %>% 
                                left_join(Volatilome %>% 
                                            spread(key = X1, value =  X2) %>% 
                                            mutate(Common = "a") %>% 
                                            rename(Sample_name = `Sample descriptor`)) %>% 
                                dplyr::select(-Common),
                              "B" = New_tibble_for_everything
  )
  
  New_tibble_joint %>% return()
}

