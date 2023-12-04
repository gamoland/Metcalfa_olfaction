permutation_STAGE_NR <-
  function(Input_data) {
    
    mean_diff_MV_STAGE <-
      abs(mean(Input_data$Normalized_response[Input_data$Variable_1 ==
                                   "adult"]) -
            mean(Input_data$Normalized_response[Input_data$Variable_1 ==
                                     "nimph"]))
    median_diff_MV_STAGE <-
      abs(median(Input_data$Normalized_response[Input_data$Variable_1 ==
                                     "adult"]) -
            median(Input_data$Normalized_response[Input_data$Variable_1 ==
                                       "nimph"]))
    
    set.seed(2020)
    N <- length(Input_data$Variable_1)
    P <- 100000
    responses <- Input_data$Normalized_response
    PermSamples <- matrix(0, nrow = N, ncol = P)
    
    for (i in 1:P) {
      PermSamples[, i] <- sample(responses,
                                 size = N,
                                 replace = FALSE)
    }
    PermSamples[, 1:5]
    
    Perm_diff_mean_MV_STAGE <- rep(0, P)
    Perm_diff_median_MV_STAGE <- rep(0, P)
    for (i in 1:P) {
      Perm_diff_mean_MV_STAGE[i] <-
        abs(mean(PermSamples[Input_data$Variable_1 == "adult", i]) -
              mean(PermSamples[Input_data$Variable_1 == "nimph", i]))
      
      Perm_diff_median_MV_STAGE[i] <-
        abs(median(PermSamples[Input_data$Variable_1 == "adult", i]) -
              median(PermSamples[Input_data$Variable_1 == "nimph" , i]))
    }
    round(Perm_diff_mean_MV_STAGE[1:15], 1)
    
    p_mean_MV_STAGE <-
      mean(Perm_diff_mean_MV_STAGE >= mean_diff_MV_STAGE)
    p_median_MV_STAGE <-
      mean(Perm_diff_median_MV_STAGE >= median_diff_MV_STAGE)
    
    # Perm_diff_mean_MV_STAGE_df <- data.frame(difs = Perm_diff_mean_MV_STAGE)
    # Perm_diff_median_MV_STAGE_df <- data.frame(difs = Perm_diff_median_MV_STAGE)
    #
    # mean_boxplot <- Perm_diff_mean_MV_STAGE_df %>%
    #   ggplot(aes(x = difs)) +
    #   geom_histogram(color = "black",
    #                  fill = "pink",
    #                  alpha = .4) +
    #   geom_vline(
    #     color = "navy",
    #     lwd = 1,
    #     lty = 2,
    #     xintercept = 0.1650
    #   ) +
    #   # xlim(-1, 1)
    #   theme_classic() +
    #   ggtitle("Mean Differences from 10000 Permutations of Raw Data_STAGE")
    # print(mean_boxplot)
    #
    # median_boxplot <- Perm_diff_median_MV_STAGE_df %>%
    #   ggplot(aes(x = difs)) +
    #   geom_histogram(color = "black",
    #                  fill = "pink",
    #                  alpha = .4) +
    #   geom_vline(
    #     color = "navy",
    #     lwd = 1,
    #     lty = 2,
    #     xintercept = 0.1525
    #   ) +
    #   # xlim(-1, 1)
    #   theme_classic() +
    #   ggtitle("Median Differences from 10000 Permutations of Raw Data_STAGE")
    # print(median_boxplot)
    
    return(
      list(
        "median_diff_MV_STAGE" = median_diff_MV_STAGE,
        "mean_diff_MV_STAGE" = mean_diff_MV_STAGE,
        "p_mean_MV_STAGE" = p_mean_MV_STAGE,
        "p_median_MV_STAGE" = p_median_MV_STAGE
      )
    )
  }