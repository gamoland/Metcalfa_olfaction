############## normalization function v1 ###############

# This function nornalizes GC minivolt readings

# ist normalization step
NORmalization <-function(Read_Data_In){ # function to read in data
  
Normalization_Data <- Read_Data_In %>% # create object for the normalized data
  mutate(Abs_mv = abs(MV), # add a variable Abs_mv to the data: the absoulute values of minivolt
         Log_abs_mv = if_else(Abs_mv == 0, 0, log(Abs_mv))) %>%  # taking the log of absollute value of minivolt
  group_by(Species, Sex, Recording_point, Organ, Rep, Variable_1) %>% # grouping
  mutate(Average_normalized_response = exp(mean(Log_abs_mv))) %>% # adding a new variable to the data:the 
                                                                  # i.e taking the exponent of the mean of the log of absoulute vallue of minivolt
  ungroup() # ungrouping the grouped data
# 2nd normalization step: normalizing the normalized_data from step 1
Normalization_Data_1 <- Normalization_Data %>% # creating new object from the normalized data
  mutate(Normalized_for_step4_1 = Abs_mv/Average_normalized_response) %>% # add new variable(Normalized_for_step4_1) which 
                                                                       # equals dividing absolute minivolt by
                                                                       # average normalised response from step 1
  group_by(Species, Sex, Recording_point, Organ, Rep, Variable_1) %>%  #grouping the data
  mutate(Average_for_step4_2 = mean(Normalized_for_step4_1)) %>%  # taking the mean  of Normalized_for_step4_1           
  ungroup() # ungrouping data
Normalization_Data_1 <- Normalization_Data %>% 
  mutate(Normalized_for_step4_1 = Abs_mv/Average_normalized_response) %>% 
  group_by(Species, Sex, Recording_point, Organ, Rep, Variable_1) %>% 
  mutate(Average_for_step4_2 = mean(Normalized_for_step4_1)) %>% 
  ungroup()

Normalization_Data_2 <- Normalization_Data_1 %>% 
  mutate(Normalized_response = (Abs_mv/Average_normalized_response)/Average_for_step4_2) %>% 
  dplyr::select(-Abs_mv, -Log_abs_mv, -Average_normalized_response, -Normalized_for_step4_1, -Average_for_step4_2)
Normalization_Data_2 %>% return()
# 3rd normalization step
Normalization_Data_2 <- Normalization_Data_1 %>% # create new object from the 2nd normalization step
  mutate(Normalized_response = (Abs_mv/Average_normalized_response)/Average_for_step4_2) %>% # add a new variable which is equals to                                                                            # dividing MV by average_normalized_response of step 1 and dividing                                                                            # that by the normalized response of step 2 (Average_for_step4_3)
  dplyr::select(-Abs_mv, -Log_abs_mv, -Average_normalized_response, -Normalized_for_step4_1, -Average_for_step4_2) 
Normalization_Data_2 %>% return() # return normalization_Data_2 without the variable from the line above. (-) sign
}



