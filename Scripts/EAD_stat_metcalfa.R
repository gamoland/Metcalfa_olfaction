library(stringi)
library(tidyverse)
library(tidyr)
library(xlsx)
library(chemodiv)
library(ggplot2)
library(ggplus)
library(ggforce)
library(gmodels)
library(DescTools)
library(qqplotr)
library(rstatix)
library(coin)
library(ggpubr)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
library(fitdistrplus)
library(logspline)
library("lme4")
library("MASS")
library("car")
library(emmeans)
library(DHARMa)
library("grafify")

#glm source
source("Scripts/glm.R")
#PERMUTATION-TEST
source("Scripts/permutation_test.R")
#MANN-WHITNEY TEST, WILCOXON SUM RANK TEST
source("Scripts/mannwhitney_test.R")



#READING EVERYTHING IN
source("Scripts/Reading_in_data.R")


#CREATING TABLE FOR THE COMPOUNDS FROM TAGETES
Compounds_TAG <- Normalized_Data_Metc_TAG %>%
  group_by(Compound.Name) %>%
  summarise() %>%
  as.tibble()


#DESCRIPTIVE STAT FOR FEMALE-MALE AND ADULT-NIMPH COMPARISION
Descriptives_all_compound_dev_stage <- tibble()
Descriptives_all_compound_sex <- tibble()

a <- 1
for (a in 1:nrow(Compounds_TAG)) {
  Data_per_compound <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == Compounds_TAG$Compound.Name[a]) %>%
    group_by(Compound.Name, Sex, Variable_1, Normalized_response) %>%
    summarise()
  
  Descriptives_dev_stage <- Data_per_compound %>%
    group_by(Variable_1) %>%
    summarise(
      n = n(),
      mean = mean(Normalized_response, na.rm = TRUE),
      sd = sd(Normalized_response, na.rm = TRUE),
      stderr = sd / sqrt(n),
      LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
      UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
      median = median(Normalized_response, na.rm = TRUE),
      min = min(Normalized_response, na.rm = TRUE),
      max = max(Normalized_response, na.rm = TRUE),
      IQR = IQR(Normalized_response, na.rm = TRUE),
      LCLmed = MedianCI(Normalized_response, na.rm = TRUE)[2],
      UCLmed = MedianCI(Normalized_response, na.rm = TRUE)[3]
    ) %>%
    mutate(Compound = Compounds_TAG$Compound.Name[a])
  
  Descriptives_all_compound_dev_stage <-
    Descriptives_all_compound_dev_stage %>%
    bind_rows(Descriptives_dev_stage)
  
  Descriptives_sex <- Data_per_compound %>%
    group_by(Sex) %>%
    summarise(
      n = n(),
      mean = mean(Normalized_response, na.rm = TRUE),
      sd = sd(Normalized_response, na.rm = TRUE),
      stderr = sd / sqrt(n),
      LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
      UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
      median = median(Normalized_response, na.rm = TRUE),
      min = min(Normalized_response, na.rm = TRUE),
      max = max(Normalized_response, na.rm = TRUE),
      IQR = IQR(Normalized_response, na.rm = TRUE),
      LCLmed = MedianCI(Normalized_response, na.rm = TRUE)[2],
      UCLmed = MedianCI(Normalized_response, na.rm = TRUE)[3]
    ) %>%
    mutate(Compound = Compounds_TAG$Compound.Name[a])
  
  Descriptives_all_compound_sex <-
    Descriptives_all_compound_sex %>%
    bind_rows(Descriptives_sex) %>%
    filter(Sex != "nimph")
  
  Boxplot_dev_stage <- Data_per_compound %>%
    ggplot(aes(x = Variable_1, y = Normalized_response, fill = Variable_1)) +
    stat_boxplot(geom = "errorbar", width = 0.5) +
    geom_boxplot(fill = "light blue") +
    geom_jitter(width = .1, size = 2) +
    stat_summary(
      fun.y = mean,
      geom = "point",
      shape = 10,
      size = 3.5,
      color = "black"
    ) +
    ggtitle(paste("Boxplot of", Compounds_TAG$Compound.Name[a])) +
    theme_bw() + theme(legend.position = "none")
  print(Boxplot_dev_stage)
  
  Boxplot_sex <- Data_per_compound %>%
    filter(Sex != "nimph") %>%
    ggplot(aes(x = Sex, y = Normalized_response, fill = Sex)) +
    stat_boxplot(geom = "errorbar", width = 0.5) +
    geom_boxplot(fill = "light blue") +
    geom_jitter(width = .1, size = 2) +
    stat_summary(
      fun.y = mean,
      geom = "point",
      shape = 10,
      size = 3.5,
      color = "black"
    ) +
    ggtitle(paste("Boxplot of", Compounds_TAG$Compound.Name[a])) +
    theme_bw() + theme(legend.position = "none")
  print(Boxplot_sex)
  
  QQplot_dev_stage <- Data_per_compound %>%
    ggplot(mapping = aes(
      sample = Normalized_response,
      color = Variable_1,
      fill = Variable_1
    )) +
    stat_qq_band(
      alpha = 0.5,
      conf = 0.95,
      qtype = 1,
      bandType = "boot"
    ) +
    stat_qq_line(identity = TRUE) +
    stat_qq_point(col = "black") +
    facet_wrap( ~ Variable_1, scales = "free") +
    stat_qq_point(col = "black") +
    facet_wrap( ~ Sex, scales = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_bw() +
    ggtitle(paste("QQplot of", Compounds_TAG$Compound.Name[a]))
  
  QQplot_sex <- Data_per_compound %>%
    filter(Sex != "nimph") %>%
    ggplot(mapping = aes(
      sample = Normalized_response,
      color = Sex,
      fill = Sex
    )) +
    stat_qq_band(
      alpha = 0.5,
      conf = 0.95,
      qtype = 1,
      bandType = "boot"
    ) +
    stat_qq_line(identity = TRUE) +
    stat_qq_point(col = "black") +
    facet_wrap( ~ Variable_1, scales = "free") +
    stat_qq_point(col = "black") +
    facet_wrap( ~ Sex, scales = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_bw() +
    ggtitle(paste("QQplot of", Compounds_TAG$Compound.Name[a]))
  
  print(Compounds_TAG$Compound.Name[a])
}


#FISHER EXACT TEST
Data_for_fisherTest <- NAs_Normalized_Data_Metc_TAG %>%
  group_by(Compound.Name, Sex, Variable_1, MV_fisher) %>%
  summarise(n = n())


Table_for_p_correction_FISHER_STAGE <- tibble()
Table_for_p_correction_FISHER_SEX <- tibble()

Table_for_choosing_compounds_SEX <- tibble()
Table_for_choosing_compounds_STAGE <- tibble()


Table_for_choosing_compounds_A <- tibble()

file.create("Output/fisher_output.txt") # creating a file for the output

i <- 1
for (i in 1:nrow(Compounds_TAG)) {
  #Creating the matrix for fisher-test
  Data_per_compound_fisher_STAGE_for_A <- Data_for_fisherTest %>%
    filter(Compound.Name == Compounds_TAG$Compound.Name[17]) %>%
    spread(key = MV_fisher, value = n, sep = "")  %>%
    mutate(
      MV_fisher0 = ifelse("MV_fisher0" %in% names(.), ifelse(is.na(MV_fisher0), 0, MV_fisher0), 0),
      MV_fisher1 = ifelse("MV_fisher1" %in% names(.), ifelse(is.na(MV_fisher1), 0, MV_fisher1), 0)
    )
  
  
  Table_for_choosing_compounds_A <-
    Table_for_choosing_compounds_A %>%
    bind_rows(
      Data_per_compound_fisher_STAGE_for_A %>%
        pivot_wider(
          names_from = Sex,
          values_from = MV_fisher1,
          id_cols = Compound.Name
        )
    )
  
  Data_per_compound_fisher_STAGE <-
    Data_per_compound_fisher_STAGE_for_A %>%
    group_by(Variable_1) %>%
    summarise(
      individuals_yes = sum(MV_fisher1),
      individuals_non = sum(MV_fisher0)) %>%
    column_to_rownames("Variable_1")
  
  
  Data_per_compound_fisher_SEX <- Data_for_fisherTest %>%
    filter(Compound.Name == Compounds_TAG$Compound.Name[17],
           Sex != "nimph") %>%
    spread(key = MV_fisher, value = n, sep = "")  %>%
    mutate(
      MV_fisher0 = ifelse("MV_fisher0" %in% names(.), ifelse(is.na(MV_fisher0), 0, MV_fisher0), 0),
      MV_fisher1 = ifelse("MV_fisher1" %in% names(.), ifelse(is.na(MV_fisher1), 0, MV_fisher1), 0)
    ) %>%
    group_by(Sex) %>%
    summarise(
    individuals_yes = sum(MV_fisher1),
    individuals_non = sum(MV_fisher0)
    ) %>%
    column_to_rownames("Sex")
  
  #fisher-test
  
  Fisher_tests_STAGE <- fisher_test(Data_per_compound_fisher_STAGE,
                                    alternative = "two.sided",
                                    detailed = T) %>%
    add_significance() %>%
    add_column(.before = "n", compound = Compounds_TAG$Compound.Name[i])
  Table_for_p_correction_FISHER_STAGE <-
    Table_for_p_correction_FISHER_STAGE %>%
    bind_rows(Fisher_tests_STAGE)
  
  Fisher_tests_SEX <- fisher_test(Data_per_compound_fisher_SEX,
                                  alternative = "two.sided",
                                  detailed = T) %>%
    add_significance() %>%
    add_column(.before = "n", compound = Compounds_TAG$Compound.Name[i])
  Table_for_p_correction_FISHER_SEX <-
    Table_for_p_correction_FISHER_SEX %>%
    bind_rows(Fisher_tests_SEX)
  
  
  #sinking the required outputs into a .txt file; not a nice solution, but it works...sry Seb
  sink(file = "Output/fisher_output.txt", append = T)
  print(Compounds_TAG$Compound.Name[i])
  print(Data_per_compound_fisher_STAGE)
  print(Fisher_tests_STAGE)
  print(Data_per_compound_fisher_SEX)
  print(Fisher_tests_SEX)
  cat("\n")
  cat("\n")
  sink(file = NULL)
  
  Data_bind_SEX <- Data_per_compound_fisher_SEX %>%
    summarise(across(c(individuals_non:individuals_yes), sum)) %>%
    add_column(.before = "individuals_non", compound = Compounds_TAG$Compound.Name[i])
  Table_for_choosing_compounds_SEX <-
    Table_for_choosing_compounds_SEX %>%
    bind_rows(Data_bind_SEX)
  
  Data_bind_STAGE <- Data_per_compound_fisher_STAGE %>%
    summarise(across(c(individuals_non:individuals_yes), sum)) %>%
    add_column(.before = "individuals_non", compound = Compounds_TAG$Compound.Name[i])
  Table_for_choosing_compounds_STAGE <-
    Table_for_choosing_compounds_STAGE %>%
    bind_rows(Data_bind_STAGE)
  
  print(i)
}

Table_for_p_correction_FISHER_STAGE <- Table_for_p_correction_FISHER_STAGE %>%
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

Table_for_p_correction_FISHER_SEX <- Table_for_p_correction_FISHER_SEX %>%
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

write.xlsx(Table_for_p_correction_FISHER_SEX, file = "Output/FISHER_SEX.xlsx", )
write.xlsx(Table_for_p_correction_FISHER_STAGE, file = "Output/FISHER_STAGE.xlsx")

++#MIXED MODEL
#ANOVA

Normalized_Data_Metc_TAG_lme <- Normalized_Data_Metc_TAG %>%
  mutate(
    Sex = as.factor(Sex),
    Variable_1 = as.factor(Variable_1),
    MV_abs = abs(MV),
    Rep = as.factor(Rep),
    Normalized_response_log = log(Normalized_response)
  ) %>% 
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

Normalized_Data_Metc_TAG_lme_SEX <- Normalized_Data_Metc_TAG_lme 

Rep_number_lmer <- Normalized_Data_Metc_TAG_lme_SEX %>% 
  group_by(Compound.Name, Sex) %>% 
  summarise(n = n())

write.xlsx(x = Rep_number_lmer, file = "Output/Rep_number.xlsx")

qqp(Normalized_Data_Metc_TAG_lme_SEX$Normalized_response_log, "norm")
qqp(Normalized_Data_Metc_TAG_lme_SEX$Normalized_response, "lnorm")
poisson <- fitdistr(Normalized_Data_Metc_TAG_lme_SEX$MV_abs, "Poisson")
qqp(Normalized_Data_Metc_TAG_lme_SEX$MV_abs, "pois", lambda = poisson$estimate,)
gamma <- fitdistr(Normalized_Data_Metc_TAG_lme_SEX$Normalized_response, "gamma")
qqp(Normalized_Data_Metc_TAG_lme_SEX$MV_abs, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])


#assumptions-outlier
Outliers_anova_TAG_SEX <- Normalized_Data_Metc_TAG_lme_SEX %>%
  group_by(Compound.Name) %>% 
  identify_outliers(Normalized_response_log)

#normality of residuals&homogen
LME_SEX_Compound<- lmer(data = Normalized_Data_Metc_TAG_lme_SEX, 
                  formula = Normalized_response_log ~  Compound.Name + (1 | Rep))
simulationOutput_lme <- simulateResiduals(fittedModel = LME_SEX_Compound)
plot(simulationOutput_lme)
anova(LME_SEX_Compound)
LMER_results <- summary(LME_SEX_Compound)

ls1 = lsmeans(LME_SEX_Compound, pairwise~Compound.Name)
ls1$contrasts
A_confint_84 <- confint(ls1, level=0.84)[1]

LME_TEST <- lme(data = Normalized_Data_Metc_TAG_lme_SEX, 
                fixed = Normalized_response_log ~  Compound.Name, 
                random = ~ 1 | Rep, 
                na.action = na.omit)
summary(LME_TEST)
anova(LME_TEST)


posthoc2 <- emmeans(LME_SEX_Compound, "Compound.Name")
summary(posthoc2, adjust = "fdr", infer = c(TRUE, FALSE), level = .84)

PostHoc_results <- posthoc_Pairwise(Model = LME_SEX_Compound, Fixed_Factor = "Compound.Name", P_Adj = "fdr", level = .84)
summary(PostHoc_results, level = .84)
# write.xlsx(x = PostHoc_results$contrasts, file = "Output/posthoc.xlsx")

A_confint_84_posthoc <- confint(PostHoc_results, level=0.84)[1]

# CIs <- intervals(LME_TEST, level = .84, which = "fixed")

# CIs_ <- as.tibble(CIs[["fixed"]]) %>% 
  # mutate(Compound.Name = c("1-hexanol",
  #                          "(Z)-3-Hexen-1-ol",
  #                          "(Z)-beta-Ocimene",
  #                          "ß-caryophyllene",
  #                          "DMNT",
  #                          "MeSa",
  #                          "myroxide",
  #                          "piperitone"))
# write.xlsx(CIs_, file = "Output/CIs.xlsx")

ConfIntVals <- read.xlsx(xlsxFile = "Output/CIs.xlsx", sheet = "readin", rows = 1:9, cols = 1:4 )

CI_boxplot <- Normalized_Data_Metc_TAG_lme_SEX %>% 
  dplyr::select(Compound.Name, Sex, Variable_1, Rep, Normalized_response_log, Normalized_response) %>% 
  group_by(Compound.Name) %>% 
  mutate(MEAN = mean(Normalized_response_log)) %>% 
  # ungroup() %>% 
  # group_by(Compound.Name, MEAN) %>% 
  # summarise() %>% 
  # mutate(Compound.Name = factor(Compound.Name, levels = Compound.Name[order(MEAN, decreasing = T)])) %>% 
  # # mutate(y0 = quantile(Normalized_response_log, 0.05, na.rm = T), 
  #        y16 = quantile(Normalized_response_log, 0.42, na.rm = T),
  #        y50 = median(Normalized_response_log, na.rm = T), 
  #        y84 = quantile(Normalized_response_log, 0.42, na.rm = T), 
  #        y100 = quantile(Normalized_response_log, 0.95, na.rm = T)) %>% 
  # ungroup() %>% 
  left_join(ConfIntVals) %>% 
  ungroup()
  

results_boxplot <- CI_boxplot %>% 
  ggplot(aes(x=fct_rev(fct_reorder(Compound.Name, MEAN)), y = Normalized_response_log)) +
  # geom_boxplot(aes(ymin = y0, lower = y16, middle = y50, upper = y84, ymax = y100), stat = "identity") +
  geom_jitter(width = .15, size = 1, color = "#343a40") +
  # geom_errorbar(aes(ymin = y16, ymax = y84), width = 0.2) +
  geom_linerange(aes(ymin = lower, ymax = upper), color = "#a4133c", size = .8) +
  stat_summary(color = "#a4133c") +
  labs(title = "Difference between responses", x = "compound", y = "response strenght (mV)") + 
  annotate("text", x=1, y=-2, label= "A") +
  annotate("text", x=2, y=-2, label= "A") +
  annotate("text", x=3, y=-2, label= "A") +
  annotate("text", x=4, y=-2, label= "B") +
  annotate("text", x=5, y=-2, label= "B") +
  annotate("text", x=6, y=-2, label= "B") +
  annotate("text", x=7, y=-2, label= "B") +
  annotate("text", x=8, y=-2, label= "C") +
  # ggtitle("Title") +
  # (main = "blabla") +
  # geom_linerange(aes(ymin = lower, ymax = upper), color = "black", size = .8) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 15))
results_boxplot

  ggsave("Output/lmer_plot_all.pdf",width = 300, height = 150, units = "mm" )

#Gergőnek
LME_TEST_gergo <- lme(data = NAs_Normalized_Data_Metc_TAG_lme, 
                fixed = Normalized_response_log ~  Sex * Compound.Name, 
                random = ~ 1 | Rep, 
                na.action = na.omit)
LME_SEX_gergo<- lmer(data = NAs_Normalized_Data_Metc_TAG_lme_SEX, 
                        formula = Normalized_response_log ~  Compound.Name * Sex + (1 | Rep))
simulationOutput_lme_gergo <- simulateResiduals(fittedModel = LME_SEX_gergo)
plot(simulationOutput_lme_gergo)
anova(LME_SEX_gergo)
summary(LME_TEST_gergo)
anova(LME_TEST_gergo)
