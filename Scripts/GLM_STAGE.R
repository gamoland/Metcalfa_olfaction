library("onewaytests")

glm_stage <- function() {
  GLM_data_STAGE <-  tibble()
  #1HEXANOL
  Responses_1HEXOL_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "1-Hexanol") %>%
    group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_1HEXOL <- fitDist(
    data = Responses_1HEXOL_STAGE,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_1HEXOL)
  gamlsstest_1HEXOL <- gamlss(
    data = na.omit(Responses_1HEXOL_STAGE),
    formula = Normalized_response ~ Variable_1,
    sigma.formula = ~ Variable_1,
    family = IG
  )
  summary(gamlsstest_1HEXOL)
  
  plot <-
    ggplot(Responses_1HEXOL_STAGE,
           aes(x = Variable_1, y = Normalized_response)) +
    geom_point()
  plot
  
  bf.test(Normalized_response ~ Variable_1, data = Responses_1HEXOL_STAGE)
  
  glmtest_1HEXOL <- glm(
    data = Responses_1HEXOL_STAGE,
    formula = Normalized_response ~ Variable_1,
    family = inverse.gaussian("1/mu^2")
  )
  summary(glmtest_1HEXOL)
  par(mfrow = c(2, 2))
  plot(glmtest_1HEXOL)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_1HEXOL))[2, ])) %>%
    add_column(.before = "Estimate", compound = "1-Hexanol")
  GLM_data_STAGE <- GLM_data_STAGE %>%
    bind_rows(Data_to_bind)
  
  #3-Hexen-1-ol (Z)-
  Responses_3HEX1OL_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "3-Hexen-1-ol (Z)-") %>%
    group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_3HEX1OL <- fitDist(
    data = Responses_3HEX1OL_STAGE,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_3HEX1OL)
  
  plot <-
    ggplot(Responses_3HEX1OL_STAGE,
           aes(x = Variable_1, y = Normalized_response)) +
    geom_point()
  
  
  bf.test(Normalized_response ~ Variable_1, data = Responses_3HEX1OL_STAGE)
  
  glmtest_3HEX1OL <- glm(data = Responses_3HEX1OL_STAGE,
                         formula = Normalized_response ~ Variable_1,
                         Gamma(link = "log"))
  summary(glmtest_3HEX1OL)
  par(mfrow = c(2, 2))
  plot(glmtest_3HEX1OL)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_3HEX1OL))[2, ])) %>%
    add_column(.before = "Estimate", compound = "3-Hexen-1-ol (Z)-")
  GLM_data_STAGE <- GLM_data_STAGE %>%
    bind_rows(Data_to_bind)
  
  # beta-cis-Ocimene
  Responses_BcisO_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "beta-cis-Ocimene") %>%
    group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_BcisO <- fitDist(
    data = Responses_BcisO_STAGE,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_BcisO)
  
  plot <-
    ggplot(Responses_BcisO_STAGE,
           aes(x = Variable_1, y = Normalized_response)) +
    geom_point()
  
  bf.test(Normalized_response ~ Variable_1, data = Responses_BcisO_STAGE)
  
  
  glmtest_BcisO <- glm(data = Responses_BcisO_STAGE,
                       formula = Normalized_response ~ Variable_1)
  summary(glmtest_BcisO)
  par(mfrow = c(2, 2))
  plot(glmtest_BcisO)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_BcisO))[2, ])) %>%
    add_column(.before = "Estimate", compound = "beta-cis-Ocimene")
  GLM_data_STAGE <- GLM_data_STAGE %>%
    bind_rows(Data_to_bind)
  
  
  #beta-trans-Ocimene
  Responses_BtransO_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "beta-trans-Ocimene") %>%
    group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response),
           Sex = as.factor(Sex),
           Variable_1 = as.factor(Variable_1)) %>% 
    na.omit
  Dist_for_data_BtransO <- fitDist(
    data = Responses_BtransO_STAGE,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_BtransO)
  
  plot <-
    ggplot(Responses_BtransO_STAGE,
           aes(x = Variable_1, y = Normalized_response)) +
    geom_point()
  plot
  
  bf.test(Normalized_response ~ Variable_1, data = Responses_BtransO_STAGE)
  
  
  glmtest_BtransO <- glm(data = Responses_BtransO_STAGE,
                         formula = Normalized_response ~ Variable_1,
                         family = inverse.gaussian(link = "identity"))
  summary(glmtest_BtransO)
  par(mfrow = c(2, 2))
  plot(glmtest_BtransO)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_BtransO))[2, ])) %>%
    add_column(.before = "Estimate", compound = "beta-trans-Ocimene")
  GLM_data_STAGE <- GLM_data_STAGE %>%
    bind_rows(Data_to_bind)
  
  
  #Caryophyllene
  Responses_CARY_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "Caryophyllene") %>%
    group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_CARY <- fitDist(
    data = Responses_CARY_STAGE,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_CARY)
  
  plot <-
    ggplot(Responses_CARY_STAGE,
           aes(x = Variable_1, y = Normalized_response)) +
    geom_point()
  plot
  
  bf.test(Normalized_response ~ Variable_1, data = Responses_CARY_STAGE)
  
  glmtest_CARY <- glm(
    data = Responses_CARY_STAGE,
    formula = log(Normalized_response) ~ Variable_1,
    gaussian(link = "identity")
  )
  summary(glmtest_CARY)
  par(mfrow = c(2, 2))
  plot(glmtest_CARY)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_CARY))[2, ])) %>%
    add_column(.before = "Estimate", compound = "Caryophyllene")
  GLM_data_STAGE <- GLM_data_STAGE %>%
    bind_rows(Data_to_bind)
  
  
  #DMNT
  Responses_DMNT_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "DMNT") %>%
    group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_DMNT <- fitDist(
    data = Responses_DMNT_STAGE,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_DMNT)
  
  bf.test(Normalized_response ~ Variable_1, data = Responses_DMNT_STAGE)
  
  plot <-
    ggplot(Responses_DMNT_STAGE,
           aes(x = Variable_1, y = Normalized_response)) +
    geom_point()
  plot
  
  glmtest_DMNT <- glm(
    data = Responses_DMNT_STAGE,
    formula = Normalized_response ~ Variable_1,
    quasi(link = "log", variance = "constant")
  )
  summary(glmtest_DMNT)
  par(mfrow = c(2, 2))
  plot(glmtest_DMNT)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_DMNT))[2, ])) %>%
    add_column(.before = "Estimate", compound = "DMNT")
  GLM_data_STAGE <- GLM_data_STAGE %>%
    bind_rows(Data_to_bind)
  
  
  # #Germacrene D
  # Responses_GERD_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
  #   filter(Compound.Name == "Germacrene D") %>%
  #   group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
  #   summarise() %>%
  #   mutate(Normalized_response = as.vector(Normalized_response))
  # Dist_for_data_GERD <- fitDist(data = Responses_GERD_STAGE,
  #                               y = Normalized_response,
  #                               k = 2, type = "realAll",
  #                               trace = FALSE,
  #                               try.gamlss = TRUE)
  # summary(Dist_for_data_GERD)
  #
  # plot <- ggplot(Responses_GERD_STAGE, aes(x = Variable_1, y = Normalized_response)) +
  #   geom_point()
  # plot
  #
  # glmtest_GERD <- glm(data = Responses_GERD_STAGE,
  #                     formula = Normalized_response ~ Variable_1,
  #                     gaussian(link = "log"))
  # summary(glmtest_GERD)
  # par(mfrow = c(2, 2))
  # plot(glmtest_GERD)
  # Data_to_bind <- as.data.frame(t(coef(summary(glmtest_GERD))[2,])) %>%
  #   add_column(.before = "Estimate", compound = "Germacrene D")
  # GLM_data_STAGE <- GLM_data_STAGE %>%
  #   bind_rows(Data_to_bind)
  
  
  #linalool
  # Responses_linalool_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
  #   filter(Compound.Name == "linalool") %>%
  #   group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
  #   summarise() %>%
  #   mutate(Normalized_response = as.vector(Normalized_response))
  # Dist_for_data_linalool <- fitDist(data = Responses_linalool_STAGE,
  #                               y = Normalized_response,
  #                               k = 2, type = "realAll",
  #                               trace = FALSE,
  #                               try.gamlss = TRUE)
  # summary(Dist_for_data_linalool)
  #
  # plot <- ggplot(Responses_linalool_STAGE, aes(x = Variable_1, y = Normalized_response)) +
  #   geom_point()
  # plot
  #
  # glmtest_linalool <- glm(data = Responses_linalool_STAGE,
  #                     formula = Normalized_response ~ Variable_1,
  #                     quasi(link=power(1/3)))
  # summary(glmtest_linalool)
  # par(mfrow = c(2, 2))
  # plot(glmtest_linalool)
  # Data_to_bind <- as.data.frame(t(coef(summary(glmtest_linalool))[2,])) %>%
  #   add_column(.before = "Estimate", compound = "linalool")
  # GLM_data_STAGE <- GLM_data_STAGE %>%
  #   bind_rows(Data_to_bind)
  
  
  #MeSa
  Responses_MeSa_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "Methyl salicylate") %>%
    group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_MeSa <- fitDist(
    data = Responses_MeSa_STAGE,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_MeSa)
  
  plot <-
    ggplot(Responses_MeSa_STAGE,
           aes(x = Variable_1, y = Normalized_response)) +
    geom_point()
  plot
  
  bf.test(Normalized_response ~ Variable_1, data = Responses_MeSa_STAGE)
  
  
  glmtest_MeSa <- glm(data = Responses_MeSa_STAGE,
                      formula = Normalized_response ~ Variable_1,
                      gaussian(link = "log"))
  summary(glmtest_MeSa)
  par(mfrow = c(2, 2))
  plot(glmtest_MeSa)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_MeSa))[2, ])) %>%
    add_column(.before = "Estimate", compound = "Methyl salicylate")
  GLM_data_STAGE <- GLM_data_STAGE %>%
    bind_rows(Data_to_bind)
  
  
  #Myroxide
  Responses_MYR_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "Myroxide") %>%
    group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_MYR <- fitDist(
    data = Responses_MYR_STAGE,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_MYR)
  
  plot <-
    ggplot(Responses_MYR_STAGE,
           aes(x = Variable_1, y = Normalized_response)) +
    geom_point()
  plot
  
  bf.test(Normalized_response ~ Variable_1, data = Responses_MYR_STAGE)
  
  
  glmtest_MYR <- glm(
    data = Responses_MYR_STAGE,
    formula = Normalized_response ~ Variable_1,
    quasi(link = "identity", variance = "mu^2")
  )
  summary(glmtest_MYR)
  par(mfrow = c(2, 2))
  plot(glmtest_MYR)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_MYR))[2, ])) %>%
    add_column(.before = "Estimate", compound = "Myroxide")
  GLM_data_STAGE <- GLM_data_STAGE %>%
    bind_rows(Data_to_bind)
  
  
  #Piperitone
  Responses_PIP_STAGE <- NAs_Normalized_Data_Metc_TAG %>%
    filter(Compound.Name == "Piperitone") %>%
    group_by(Variable_1, Normalized_response, Nrow, Sex) %>%
    summarise() %>%
    mutate(Normalized_response = as.vector(Normalized_response))
  Dist_for_data_PIP <- fitDist(
    data = Responses_PIP_STAGE,
    y = Normalized_response,
    k = 2,
    type = "realAll",
    trace = FALSE,
    try.gamlss = TRUE
  )
  summary(Dist_for_data_PIP)
  
  plot <-
    ggplot(Responses_PIP_STAGE,
           aes(x = Variable_1, y = Normalized_response)) +
    geom_point()
  plot
  
  bf.test(Normalized_response ~ Variable_1, data = Responses_PIP_STAGE)
  
  
  glmtest_PIP <- glm(data = Responses_PIP_STAGE,
                     formula = Normalized_response ~ Variable_1,
                     gaussian(link = "log"))
  summary(glmtest_PIP)
  par(mfrow = c(2, 2))
  plot(glmtest_PIP)
  Data_to_bind <-
    as.data.frame(t(coef(summary(glmtest_PIP))[2, ])) %>%
    add_column(.before = "Estimate", compound = "Piperitone")
  
  GLM_data_STAGE <- GLM_data_STAGE %>%
    bind_rows(Data_to_bind)
  
  GLM_data_STAGE_final <- GLM_data_STAGE %>%
    adjust_pvalue(p.col = "Pr(>|t|)", method = "fdr") %>%
    add_significance(add_significance("Pr(>|t|).adj"))
  
  return(GLM_data_STAGE_final)
}
