Npclassifier_API <- function(Input) {
  Input <- Input
  Return_list <- tibble()
  i <- 1
  for (i in 1:nrow(Input)) {
    Selected <- Input[i,] 
    print(i)
    Npclass_raw <- GET(paste0("https://npclassifier.ucsd.edu/classify?smiles=",  URLencode(Selected$CanonicalSMILES, reserved = TRUE)))
    #URLencode(Selected$CanonicalSMILES, reserved = TRUE)
    Npclass_content <-content(Npclass_raw,
                              "text", encoding = "UTF-8")
    Np_class_from_JSON <- fromJSON(Npclass_content)
    Test <- if_else(is.null(Np_class_from_JSON$pathway_results[1][[1]]),
                    "Not resolved", "Resolved")
    Np_pathway1 <- Np_class_from_JSON$pathway_results[1]
    Np_pathway2 <- Np_class_from_JSON$pathway_results[2]
    Np_pathway3 <- Np_class_from_JSON$pathway_results[3] 
    Superclass_results1 <- if_else(is.list(Np_class_from_JSON$superclass_results[1]),
                                   "", as.character(Np_class_from_JSON$superclass_results[1]))
    
    Superclass_results2 <- if_else(is.list(Np_class_from_JSON$superclass_results[2]),
                                   "", as.character(Np_class_from_JSON$superclass_results[2]))
    
    Class_results1 <- if_else(is.list(Np_class_from_JSON$class_results[1]),
                              "", as.character(Np_class_from_JSON$class_results[1]))
    Class_results2 <- if_else(is.list(Np_class_from_JSON$class_results[2]),
                              "", as.character(Np_class_from_JSON$class_results[2]))
    
    To_add <- switch (Test,
                      "Resolved" = tibble(Np_pathway1, Np_pathway2, Superclass_results1,Superclass_results2,
                                          Class_results1, Class_results2
                      ),
                      "Not resolved" = tibble(Np_pathway1 = "Not resolvable"))
    Resolved <- Selected %>% bind_cols(To_add) 
    Return_list <- Return_list %>% bind_rows(Resolved)
    i <- i + 1
  }
  return(Return_list)
}