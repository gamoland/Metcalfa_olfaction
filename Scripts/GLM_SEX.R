glm_sex <- function() {
  GLM_data_SEX <- tibble()
  #1HEXANOL
  Responses_1HEXOL_SEX <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "1-Hexanol",
           Sex != "nimph") %>%
    group_by(Sex, Normalized_response, Nrow) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_1HEXOL_SEX <- fitDist(
    data = Responses_1HEXOL_SEX,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_1HEXOL_SEX)
  gamlsstest_1HEXOL_SEX <-
    gamlss(
      data = na.omit(Responses_1HEXOL_SEX),
      formula = Normalized_response ~ Sex,
      sigma.formula = ~ Sex,
      family = PE
    )
  summary(gamlsstest_1HEXOL_SEX)
  
  plot <-
    ggplot(Responses_1HEXOL_SEX, aes(x = Sex, y = Normalized_response)) +
    geom_point()
  plot
  
  glmtest_1HEXOL_SEX <- glm(
    data = Responses_1HEXOL_SEX,
    formula = log(Normalized_response) ~ Sex,
    family = gaussian(link = "identity")
  )
  summary(glmtest_1HEXOL_SEX)
  par(mfrow = c(2, 2))
  plot(glmtest_1HEXOL_SEX)
  Data_to_bind <-
    as.data.frame(t(coef(summary(
      glmtest_1HEXOL_SEX
    ))[2, ])) %>%
    add_column(.before = "Estimate", compound = "1-Hexanol")
  GLM_data_SEX <- GLM_data_SEX %>%
    bind_rows(Data_to_bind)
  
  #3-Hexen-1-ol (Z)-
  Responses_3HEX1OL_SEX <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "3-Hexen-1-ol (Z)-",
           Sex != "nimph") %>%
    group_by(Sex, Normalized_response, Nrow) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_3HEX1OL_SEX <- fitDist(
    data = Responses_3HEX1OL_SEX,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_3HEX1OL_SEX)
  
  plot <-
    ggplot(Responses_3HEX1OL_SEX, aes(x = Sex, y = Normalized_response)) +
    geom_point()
  
  glmtest_3HEX1OL_SEX <- glm(
    data = Responses_3HEX1OL_SEX,
    formula = Normalized_response ~ Sex,
    family = gaussian(link = "identity")
  )
  summary(glmtest_3HEX1OL_SEX)
  par(mfrow = c(2, 2))
  plot(glmtest_3HEX1OL_SEX)
  Data_to_bind <-
    as.data.frame(t(coef(summary(
      glmtest_3HEX1OL_SEX
    ))[2, ])) %>%
    add_column(.before = "Estimate", compound = "3-Hexen-1-ol (Z)-")
  GLM_data_SEX <- GLM_data_SEX %>%
    bind_rows(Data_to_bind)
  
  # beta-cis-Ocimene
  Responses_BcisO_SEX <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "beta-cis-Ocimene",
           Sex != "nimph") %>%
    group_by(Sex, Normalized_response, Nrow) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_BcisO_SEX <- fitDist(
    data = Responses_BcisO_SEX,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_BcisO_SEX)
  
  plot <-
    ggplot(Responses_BcisO_SEX, aes(x = Sex, y = Normalized_response)) +
    geom_point()
  plot
  
  glmtest_BcisO_SEX <- glm(data = Responses_BcisO_SEX,
                           formula = Normalized_response ~ Sex,
                           inverse.gaussian(link = "1/mu^2"))
  summary(glmtest_BcisO_SEX)
  par(mfrow = c(2, 2))
  plot(glmtest_BcisO_SEX)
  Data_to_bind <-
    as.data.frame(t(coef(summary(
      glmtest_BcisO_SEX
    ))[2, ])) %>%
    add_column(.before = "Estimate", compound = "beta-cis-Ocimene")
  GLM_data_SEX <- GLM_data_SEX %>%
    bind_rows(Data_to_bind)
  
  
  #beta-trans-Ocimene
  Responses_BtransO_SEX <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "beta-trans-Ocimene",
           Sex != "nimph") %>%
    group_by(Sex, Normalized_response, Nrow) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_BtransO_SEX <- fitDist(
    data = Responses_BtransO_SEX,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_BtransO_SEX)
  
  plot <-
    ggplot(Responses_BtransO_SEX, aes(x = Sex, y = Normalized_response)) +
    geom_point()
  plot
  
  glmtest_BtransO_SEX <- glm(data = Responses_BtransO_SEX,
                             formula = Normalized_response ~ Sex,
                             gaussian(link = "identity"))
  summary(glmtest_BtransO_SEX)
  par(mfrow = c(2, 2))
  plot(glmtest_BtransO_SEX)
  Data_to_bind <-
    as.data.frame(t(coef(summary(
      glmtest_BtransO_SEX
    ))[2, ])) %>%
    add_column(.before = "Estimate", compound = "beta-trans-Ocimene")
  GLM_data_SEX <- GLM_data_SEX %>%
    bind_rows(Data_to_bind)
  
  
  #Caryophyllene
  Responses_CARY_SEX <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "Caryophyllene",
           Sex != "nimph") %>%
    group_by(Sex, Normalized_response, Nrow) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_CARY_SEX <- fitDist(
    data = Responses_CARY_SEX,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_CARY_SEX)
  
  plot <-
    ggplot(Responses_CARY_SEX, aes(x = Sex, y = Normalized_response)) +
    geom_point()
  plot
  
  glmtest_CARY_SEX <- glm(
    data = Responses_CARY_SEX,
    formula = log(Normalized_response) ~ Sex,
    gaussian(link = "identity")
  )
  summary(glmtest_CARY_SEX)
  par(mfrow = c(2, 2))
  plot(glmtest_CARY_SEX)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_CARY_SEX))[2, ])) %>%
    add_column(.before = "Estimate", compound = "Caryophyllene")
  GLM_data_SEX <- GLM_data_SEX %>%
    bind_rows(Data_to_bind)
  
  
  #DMNT
  Responses_DMNT_SEX <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "DMNT",
           Sex != "nimph") %>%
    group_by(Sex, Normalized_response, Nrow) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_DMNT_SEX <- fitDist(
    data = Responses_DMNT_SEX,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_DMNT_SEX)
  
  plot <-
    ggplot(Responses_DMNT_SEX, aes(x = Sex, y = Normalized_response)) +
    geom_point()
  plot
  
  glmtest_DMNT_SEX <- glm(
    data = Responses_DMNT_SEX,
    formula = log(Normalized_response) ~ Sex,
    gaussian(link = "identity")
  )
  summary(glmtest_DMNT_SEX)
  par(mfrow = c(2, 2))
  plot(glmtest_DMNT_SEX)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_DMNT_SEX))[2, ])) %>%
    add_column(.before = "Estimate", compound = "DMNT")
  GLM_data_SEX <- GLM_data_SEX %>%
    bind_rows(Data_to_bind)
  
  
  # #Germacrene D
  # Responses_GERD_SEX <- NAs_Normalized_Data_Metc_TAG %>%
  #   filter(Compound.Name == "Germacrene D",
  #          Sex != "nimph") %>%
  #   group_by(Sex, Normalized_response, Nrow) %>%
  #   summarise() %>%
  #   mutate(Normalized_response = as.vector(Normalized_response))
  # Dist_for_data_GERD_SEX <- fitDist(data = Responses_GERD_SEX,
  #                               y = Normalized_response,
  #                               k = 2, type = "realAll",
  #                               trace = FALSE,
  #                               try.gamlss = TRUE)
  # summary(Dist_for_data_GERD_SEX)
  #
  # plot <- ggplot(Responses_GERD_SEX, aes(x = Sex, y = Normalized_response)) +
  #   geom_point()
  # plot
  #
  # glmtest_GERD_SEX <- glm(data = Responses_GERD_SEX,
  #                     formula = Normalized_response ~ Sex,
  #                     inverse.gaussian(link = "1/mu^2"))
  # summary(glmtest_GERD_SEX)
  # par(mfrow = c(2, 2))
  # plot(glmtest_GERD_SEX)
  # Data_to_bind <- as.data.frame(t(coef(summary(glmtest_GERD_SEX))[2,])) %>%
  #   add_column(.before = "Estimate", compound = "Germacrene D")
  # GLM_data_SEX <- GLM_data_SEX %>%
  #   bind_rows(Data_to_bind)
  #
  #
  # #linalool
  # Responses_linalool_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
  #   filter(Compound.Name == "linalool") %>%
  #   group_by(Variable_1, Normalized_response) %>%
  #   summarise() %>%
  #   mutate(Normalized_response = as.vector(Normalized_response))
  # Dist_for_data_linalool <- fitDist(data = Responses_linalool_STAGE,
  #                                   y = Normalized_response,
  #                                   k = 2, type = "realAll",
  #                                   trace = FALSE,
  #                                   try.gamlss = TRUE)
  # summary(Dist_for_data_linalool)
  #
  # plot <- ggplot(Responses_linalool_STAGE, aes(x = Variable_1, y = Normalized_response)) +
  #   geom_point()
  # plot
  #
  # glmtest_linalool <- glm(data = Responses_linalool_STAGE,
  #                         formula = Normalized_response ~ Variable_1,
  #                         quasi(link=power(1/3)))
  # summary(glmtest_linalool)
  # par(mfrow = c(2, 2))
  # plot(glmtest_linalool)
  # Data_to_bind <- as.data.frame(t(coef(summary(glmtest_linalool))[2,])) %>%
  #   add_column(.before = "Estimate", compound = "linalool")
  # GLM_data_STAGE <- GLM_data_STAGE %>%
  #   bind_rows(Data_to_bind)
  #
  #
  #MeSa
  Responses_MeSa_SEX <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "Methyl salicylate",
           Sex != "nimph") %>%
    group_by(Sex, Normalized_response, Nrow) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_MeSa_SEX <- fitDist(
    data = Responses_MeSa_SEX,
    y = Normalized_response,
    
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_MeSa_SEX)
  
  plot <-
    ggplot(Responses_MeSa_SEX, aes(x = Sex, y = Normalized_response)) +
    geom_point()
  plot
  
  glmtest_MeSa_SEX <- glm(
    data = Responses_MeSa_SEX,
    formula = log(Normalized_response) ~ Sex,
    gaussian(link = "identity")
  )
  summary(glmtest_MeSa_SEX)
  par(mfrow = c(2, 2))
  plot(glmtest_MeSa_SEX)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_MeSa_SEX))[2, ])) %>%
    add_column(.before = "Estimate", compound = "Methyl salicylate")
  GLM_data_SEX <- GLM_data_SEX %>%
    bind_rows(Data_to_bind)
  
  
  #Myroxide
  Responses_MYR_SEX <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "Myroxide",
           Sex != "nimph") %>%
    group_by(Sex, Normalized_response, Nrow) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_MYR_SEX <- fitDist(
    data = Responses_MYR_SEX,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_MYR_SEX)
  
  plot <-
    ggplot(Responses_MYR_SEX, aes(x = Sex, y = Normalized_response)) +
    geom_point()
  plot
  
  glmtest_MYR_SEX <- glm(data = Responses_MYR_SEX,
                         formula = Normalized_response ~ Sex,
                         gaussian(link = "identity"))
  summary(glmtest_MYR_SEX)
  par(mfrow = c(2, 2))
  plot(glmtest_MYR_SEX)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_MYR_SEX))[2, ])) %>%
    add_column(.before = "Estimate", compound = "Myroxide")
  GLM_data_SEX <- GLM_data_SEX %>%
    bind_rows(Data_to_bind)
  
  
  #Piperitone
  Responses_PIP_SEX <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "Piperitone",
           Sex != "nimph") %>%
    group_by(Sex, Normalized_response, Nrow) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_PIP_SEX <- fitDist(
    data = Responses_PIP_SEX,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_PIP_SEX)
  
  plot <-
    ggplot(Responses_PIP_SEX, aes(x = Sex, y = Normalized_response)) +
    geom_point()
  plot
  
  glmtest_PIP_SEX <- glm(data = Responses_PIP_SEX,
                         formula = Normalized_response ~ Sex,
                         gaussian(link = "identity"))
  summary(glmtest_PIP_SEX)
  par(mfrow = c(2, 2))
  plot(glmtest_PIP_SEX)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_PIP_SEX))[2, ])) %>%
    add_column(.before = "Estimate", compound = "Piperitone")
  GLM_data_SEX <- GLM_data_SEX %>%
    bind_rows(Data_to_bind)
  
  GLM_data_SEX_FINAL <- GLM_data_SEX %>%
    adjust_pvalue(p.col = "Pr(>|t|)", method = "fdr") %>%
    add_significance(add_significance("Pr(>|t|).adj"))
  
  return(GLM_data_SEX_FINAL)
}