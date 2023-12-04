source("Scripts/Metcalfa/permutation_SEX.R")
source("Scripts/Metcalfa/permutation_STAGE.R")
source("Scripts/Metcalfa/permutation_SEX_NR.R")
source("Scripts/Metcalfa/permutation_STAGE_NR.R")
NAs_Normalized_Data_Metc_TAG_perm <- NAs_Normalized_Data_Metc_TAG %>%
  mutate(
    Sex = as.factor(Sex),
    Variable_1 = as.factor(Variable_1),
    MV_abs = abs(MV),
    Rep = as.factor(Rep)
  ) %>%
  filter(MV_abs != "N/A")

Compounds_TAG_filtered <- c(
  "1-Hexanol",
  "3-Hexen-1-ol (Z)-",
  "beta-cis-Ocimene",
  "Caryophyllene",
  "DMNT",
  "Methyl salicylate",
  "Myroxide",
  "Piperitone"
)

perm_sex_all_MV <- list()
perm_stage_all_MV <- list()
perm_sex_all_NR <- list()
perm_stage_all_NR <- list()
i <- 1
for (i in 1:length(Compounds_TAG_filtered)) {
  
  Selected_compound = Compounds_TAG_filtered[[i]]
  Input_data <- NAs_Normalized_Data_Metc_TAG_perm %>%
    filter(Compound.Name == Selected_compound)
  
  perm_sex_all_MV[[Selected_compound]] <- permutation_SEX(Input_data)
  perm_stage_all_MV[[Selected_compound]] <- permutation_STAGE(Input_data)
  perm_sex_all_NR[[Selected_compound]] <- permutation_SEX_NR(Input_data)
  perm_stage_all_NR[[Selected_compound]] <- permutation_STAGE_NR(Input_data)
}

Result_NR_SEX <- matrix(ncol = 2, nrow = 8) 
Result_NR_STAGE <- matrix(ncol = 2, nrow = 8)
Result_MV_abs_SEX <- matrix(ncol = 2, nrow = 8)
Result_MV_abs_STAGE <- matrix(ncol = 2, nrow = 8)
a <- 1
for (a in 1:length(Compounds_TAG_filtered)) {
  Selected_compound = Compounds_TAG_filtered[[a]]
  
  Result_NR_SEX[a,1] <- Selected_compound
  Result_NR_SEX[a,2] <- perm_sex_all_NR[[Selected_compound]][[3]]
  
  Result_MV_abs_SEX[a,1] <- Selected_compound
  Result_MV_abs_SEX[a,2] <- perm_sex_all_MV[[Selected_compound]][[3]]
  
  
  Result_MV_abs_STAGE[a,1] <- Selected_compound
  Result_MV_abs_STAGE[a,2] <- perm_stage_all_MV[[Selected_compound]][[3]]
  
  
  Result_NR_STAGE[a,1] <- Selected_compound
  Result_NR_STAGE[a,2] <- perm_stage_all_NR[[Selected_compound]][[3]]
  
}


NAs_Normalized_Data_Metc_TAG_filtered <- NAs_Normalized_Data_Metc_TAG %>% 
  filter(Compound.Name %in% c(
    "1-Hexanol",
    "3-Hexen-1-ol (Z)-",
    "beta-cis-Ocimene",
    "Caryophyllene",
    "DMNT",
    "Methyl salicylate",
    "Myroxide",
    "Piperitone"
  ))

NAs_Normalized_Data_Metc_TAG_filtered_SEX <- NAs_Normalized_Data_Metc_TAG_mm %>% 
  filter(Sex != "nimph")

NAs_Normalized_Data_Metc_TAG_filtered <- as.matrix(NAs_Normalized_Data_Metc_TAG_filtered)
NAs_Normalized_Data_Metc_TAG_filtered_SEX <- as.tibble(NAs_Normalized_Data_Metc_TAG_filtered_SEX)


Normalized_response_m <- as.list(NAs_Normalized_Data_Metc_TAG_filtered$Normalized_response[NAs_Normalized_Data_Metc_TAG_filtered$Sex == "male"])


Perm_SEX_all_MV_COIN <- independence_test(data = NAs_Normalized_Data_Metc_TAG_filtered_SEX, 
                                          formula = Normalizes_response ~ Sex)