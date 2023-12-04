permutation_SEX_NR <-
  function(Input_data) {
    mean_diff_MV_SEX <-
      abs(mean(Input_data$Normalized_response[Input_data$Sex ==
                                   "female"]) -
            mean(Input_data$Normalized_response[Input_data$Sex ==
                                     "male"]))
    median_diff_MV_SEX <-
      abs(median(Input_data$Normalized_response[Input_data$Sex ==
                                     "female"]) -
            median(Input_data$Normalized_response[Input_data$Sex ==
                                       "male"]))
    
    set.seed(2020)
    N <- length(Input_data$Sex)
    P <- 100000
    responses <- Input_data$Normalized_response
    PermSamples <- matrix(0, nrow = N, ncol = P)
    
    for (i in 1:P) {
      PermSamples[, i] <- sample(responses,
                                 size = N,
                                 replace = FALSE)
    }
    PermSamples[, 1:5]
    
    Perm_diff_mean_MV_SEX <- rep(0, P)
    Perm_diff_median_MV_SEX <- rep(0, P)
    for (i in 1:P) {
      Perm_diff_mean_MV_SEX[i] <-
        abs(mean(PermSamples[Input_data$Sex == "female", i]) -
              mean(PermSamples[Input_data$Sex == "male", i]))
      
      Perm_diff_median_MV_SEX[i] <-
        abs(median(PermSamples[Input_data$Sex == "female", i]) -
              median(PermSamples[Input_data$Sex == "male" , i]))
    }
    round(Perm_diff_mean_MV_SEX[1:15], 1)
    
    p_mean_MV_SEX <- mean(Perm_diff_mean_MV_SEX >= mean_diff_MV_SEX)
    p_median_MV_SEX <-
      mean(Perm_diff_median_MV_SEX >= median_diff_MV_SEX)
    #
    # Perm_diff_mean_MV_SEX_df <- data.frame(difs = Perm_diff_mean_MV_SEX)
    # Perm_diff_median_MV_SEX_df <- data.frame(difs = Perm_diff_median_MV_SEX)
    #
    # mean_boxplot <- Perm_diff_mean_MV_SEX_df %>%
    #   ggplot(aes(x = difs)) +
    #   geom_histogram(color = "black",
    #                  fill = "pink",
    #                  alpha = .4) +
    #   geom_vline(
    #     color = "navy",
    #     lwd = 1,
    #     lty = 2,
    #     xintercept = 0.0404
    #   ) +
    #   # xlim(-1, 1)
    #   theme_classic() +
    #   ggtitle("Mean Differences from 10000 Permutations of Raw Data")
    # print(mean_boxplot)
    #
    # median_boxplot <- Perm_diff_median_MV_SEX_df %>%
    #   ggplot(aes(x = difs)) +
    #   geom_histogram(color = "black",
    #                  fill = "pink",
    #                  alpha = .4) +
    #   geom_vline(
    #     color = "navy",
    #     lwd = 1,
    #     lty = 2,
    #     xintercept = 0.1049
    #   ) +
    #   # xlim(-1, 1)
    #   theme_classic() +
    #   ggtitle("Median Differences from 10000 Permutations of Raw Data")
    # print(median_boxplot)
    
    return(
      list(
        "median_diff_MV_SEX" = median_diff_MV_SEX,
        "mean_diff_MV_SEX" = mean_diff_MV_SEX,
        "p_mean_MV_SEX" = p_mean_MV_SEX,
        "p_median_MV_SEX" = p_median_MV_SEX
      )
    )
  }