Data_for_mannwhitney <- NAs_Normalized_Data_Metc_TAG %>%
  # mutate(Normalized_response = ifelse(is.na(Normalized_response), 0, Normalized_response),
  mutate(Sex = as.factor(Sex),
         Variable_1 = as.factor(Variable_1))

file.create("Output/Metcalfa/MW_output.txt") # creating a file for the output

Table_for_p_correction_MW_STAGE <- tibble()
Table_for_p_correction_MW_SEX <- tibble()

b <- 1
for (b in 1:nrow(Compounds_TAG)) {
  #
  # Min_response_adult <- Data_for_mannwhitney %>%
  #   filter(Compound.Name == Compounds_TAG$Compound.Name[b],
  #          Variable_1 == "adult") %>%
  #   drop_na(Normalized_response) %>%
  #   summarise_at(vars(Normalized_response),
  #                list(min = min))
  #
  # Min_response_nimph <- Data_for_mannwhitney %>%
  #   filter(Compound.Name == Compounds_TAG$Compound.Name[b],
  #          Variable_1 == "nimph") %>%
  #   drop_na(Normalized_response) %>%
  #   summarise_at(vars(Normalized_response),
  #                list(min = min))
  #
  # Min_response_male <- Data_for_mannwhitney %>%
  #   filter(Compound.Name == Compounds_TAG$Compound.Name[b],
  #          Sex == "male") %>%
  #   drop_na(Normalized_response) %>%
  #   summarise_at(vars(Normalized_response),
  #                list(min = min))
  # Min_response_female <- Data_for_mannwhitney %>%
  #   filter(Compound.Name == Compounds_TAG$Compound.Name[b],
  #          Sex == "female") %>%
  #   drop_na(Normalized_response) %>%
  #   summarise_at(vars(Normalized_response),
  #                list(min = min))
  
  MW_data_SEX <- Data_for_mannwhitney %>%
    filter(Compound.Name == Compounds_TAG$Compound.Name[b],
           Sex != "nimph") %>%
    drop_na(Normalized_response)
  # mutate(Response = ifelse(is.na(Normalized_response),
  #                          ifelse(Sex == "male",
  #                                 Min_response_male,
  #                                 Min_response_female),
  #                          Normalized_response))
  
  MW_data_STAGE <- Data_for_mannwhitney %>%
    filter(Compound.Name == Compounds_TAG$Compound.Name[b]) %>%
    drop_na(Normalized_response)
  # mutate(Response = ifelse(is.na(Normalized_response),
  #                          ifelse(Variable_1 == "adult",
  #                                 Min_response_adult,
  #                                 Min_response_nimph),
  #                          Normalized_response))
  
  
  MW_test_SEX <- MW_data_SEX %>%
    rstatix::wilcox_test(Normalized_response ~ Sex,
                         exact = T,
                         detailed = T) %>%
    add_significance() %>%
    add_column(.before = "estimate", compound = Compounds_TAG$Compound.Name[b])
  Table_for_p_correction_MW_SEX <- Table_for_p_correction_MW_SEX %>%
    bind_rows(MW_test_SEX)
  
  MW_test_STAGE <- MW_data_STAGE %>%
    rstatix::wilcox_test(Normalized_response ~ Variable_1,
                         exact = T,
                         detailed = T) %>%
    add_significance() %>%
    add_column(.before = "estimate", compound = Compounds_TAG$Compound.Name[b])
  Table_for_p_correction_MW_STAGE <-
    Table_for_p_correction_MW_STAGE %>%
    bind_rows(MW_test_STAGE)
  
  Effect_size_SEX <- MW_data_SEX %>%
    wilcox_effsize(Normalized_response ~ Sex)
  Effect_size_STAGE <- MW_data_STAGE %>%
    wilcox_effsize(Normalized_response ~ Variable_1)
  
  sink(file = "Output/Metcalfa/MW_output.txt", append = T)
  print(Compounds_TAG$Compound.Name[b])
  print(MW_test_SEX)
  print(Effect_size_SEX)
  print(MW_test_STAGE)
  print(Effect_size_STAGE)
  cat("\n")
  cat("\n")
  sink()
  print(b)
  
}


Table_for_p_correction_MW_STAGE <- Table_for_p_correction_MW_STAGE %>%
  filter(compound %in% c("1-Hexanol",
                         "3-Hexen-1-ol (Z)-",
                         "beta-cis-Ocimene",
                         "Caryophyllene",
                         "DMNT",
                         "Methyl salicylate",
                         "Myroxide",
                         "Piperitone")) %>%
  adjust_pvalue(p.col = "p", method = "fdr") %>%
  add_significance()

Table_for_p_correction_MW_SEX <- Table_for_p_correction_MW_SEX %>%
  filter(compound %in% c("1-Hexanol",
                         "3-Hexen-1-ol (Z)-",
                         "beta-cis-Ocimene",
                         "Caryophyllene",
                         "DMNT",
                         "Methyl salicylate",
                         "Myroxide",
                         "Piperitone")) %>%
  adjust_pvalue(p.col = "p", method = "holm") %>%
  add_significance()

write.xlsx(Table_for_p_correction_MW_SEX, file = "Output/Metcalfa/MW_SEX.xlsx")
write.xlsx(Table_for_p_correction_MW_STAGE, file = "Output/Metcalfa/MW_STAGE.xlsx")
