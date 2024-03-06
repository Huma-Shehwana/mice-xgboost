rm(list = ls())

#source("Simulation_function.R")
#library(mice)     # for imputation and amputation
library(purrr)    # for functional programming
library(furrr)    # for functional futures
library(mvtnorm)  # for multivariate normal data
library(magrittr) # for pipes
library(dplyr)    # for data manipulation
library(tibble)   # for tibbles
library(mixgb)
library(microbenchmark)


seed <- 123
set.seed(seed) 
Num_ds <- 10
prop_values <- c(20, 30, 40, 50)
                 #, 60, 70,80, 90)
m=5
maxit=5


#################################################################################
########################     Step 1: Data generation    #########################
#################################################################################

sigma <- matrix(data = c(1, 0.7, 0.7, 1), 
                                ncol = 2)

simdata <- replicate(n = Num_ds, 
                      expr = mvtnorm::rmvnorm(n = 1000, 
                                             mean = c(8, 3), 
                                             sigma = sigma) %>% 
                      as_tibble() %>% # make into a tibble
                      rename(x = V1, z = V2) %>% # rename columns
                      mutate(y = 6 * x + 3 * z + rnorm(1000)), # add y
                      simplify = FALSE) # keep as list of generated sets

true_vals <- c(0,6,3)

#Regression
simdata %>% 
    map(~.x %$% # for every simulated set in simdata....
        lm(y ~ x + z) %>% # fit linear model
        coefficients) %>% # extract coefficients
  Reduce("+", .) / length(simdata) # add all and divide by length (= average)



#################################################################################
#############          2. Missing data 
#################################################################################


apply_ampute <- function(simdata, prop_value, seed) {
  simdata %>%
    furrr::future_map(function(x) {
      x %>%
        ampute(prop = prop_value / 100, mech = "MAR", type = "RIGHT") %>%
        .$amp 
    }, .options = furrr_options(seed = seed))
}

missing_MAR_list <- map(prop_values, ~ apply_ampute(simdata, .x, seed))
names(missing_MAR_list) <- prop_values


NAs_in_data <- map(missing_MAR_list, ~ map_dbl(.x, ~ sum(is.na(.x))))
print(NAs_in_data)



#################################################################################
#############          3. Imputation
#################################################################################


############################     1 - Default    ##################################

t1 <- Sys.time()

impute_MAR <- missing_MAR_list %>%
  furrr::future_map(function(mat_list) {
          lapply(mat_list, function(mat) {
                mice::mice(mat, 
                          m = m, 
                          maxit = maxit,
                          print = FALSE)
                })
    }, .options = furrr_options(seed = 123))
totalTime_default <- Sys.time() - t1


eval_default <- impute_MAR %>% 
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z)) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds })


save(eval_default, file = "MI_default.RData")


############################     2 - RF    ##################################
#rm(list = c("t1", "totalTime_default"))

t1_rf <- Sys.time()

impute_MAR <- missing_MAR_list %>%
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      mice::mice(mat, 
                 m = m, method = "rf",
                 maxit = maxit,
                 print = FALSE)
    })
  }, .options = furrr_options(seed = 123))
totalTime_RF <- Sys.time() - t1_rf


eval_RF <- impute_MAR %>% 
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z)) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(eval_RF, file = "MI_RF.RData")




############################     3 - CART   ##################################


t1_cart <- Sys.time()

impute_MAR <- missing_MAR_list %>%
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      mice::mice(mat, 
                 m = m, method = "cart",
                 maxit = maxit,
                 print = FALSE)
    })
  }, .options = furrr_options(seed = 123))

totalTime_CART <-  Sys.time() - t1_cart


eval_CART <- impute_MAR %>% 
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z)) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(eval_CART, file = "MI_CART.RData")


########################    4 - XGBoost ####################################

t1_xgb <- Sys.time()

impute_MAR <- missing_MAR_list %>%
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      mice::mice(mat, 
                 m = m, method = "xgb",
                 maxit = maxit,xgb.params=NULL,
                 print = FALSE)
    })
  }, .options = furrr_options(seed = 123))

totalTime_xgb <- Sys.time() - t1_xgb

eval_xgb <- impute_MAR %>% 
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z)) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(eval_xgb, file = "MI_xgb.RData")


########################################################################

imputed_data <- lapply(impute_MAR, function(inner_list) lapply(inner_list, '[[','imp'))
var_names <-  colnames(simdata[[1]])

# Calculation of testing Error
testError_list = list()

for(k in var_names){
            print(k)
            missing_index <- missing_MAR_list %>% 
                            map(~map(.x, ~ which(is.na(.x[, k]))))  #finds index of missing data
             actual_data <- lapply(missing_index, function(x) {
                                        lapply(1:length(x), function(i) {simdata[[i]][x[[i]],k]})}) # retrieves actual data which was amputed as missing
  
             imputed_var <- lapply(imputed_data, function(inner_list) #retrieves imputed data for each missing variable
                                        lapply(inner_list, '[[',k))
  
             testError_list[[k]] <- lapply(seq_along(imputed_var), function(i) {
                                              lapply(seq_along(imputed_var[[i]]), function(j) { 
                                                       apply(imputed_var[[i]][[j]],2, 
                                                             function(col) {as.matrix(colMeans(actual_data[[i]][[j]] - col))})})
                                                        })
          }
 
 

########################################################################
##########             Training and test arrangement       #############
########################################################################


train_error <- lapply(impute_MAR, function(inner_list) lapply(inner_list, '[[','trainError'))
train_error_f1 <- lapply(train_error, function(x){ lapply(seq_along(x), function(y) {x[[y]][!sapply(x[[y]], is.null)]})}) # removes NULL element from training data
p = list()

for(i in 1:length(var_names)){

  train_error_var <- lapply(train_error_f1, function(inner_list) lapply(inner_list, '[[',var_names[i]))
  trainingErrorArranged <- lapply(seq_along(train_error_var), function(i) {
              lapply(seq_along(train_error_var[[i]]), function(j) {
                      df <- train_error_var[[i]][[j]]
                      df$Percent <- prop_values[i]
                      df$Dataset <- j
                      df
                    })
})

trainingError_df <- do.call(rbind, do.call(c, trainingErrorArranged))
trainingError_df$variable <- var_names[i]

transformed_train <- trainingError_df %>%
          pivot_longer(cols = c("V1","V2","V3","V4","V5"), names_to = "Imputation", values_to = "Error") %>%
          mutate(ErrorType = "TrainingError")

transformed_train$Imputation <- gsub("^V", "", transformed_train$Imputation )


# test error arrangement

index <- match(var_names[i], names(testError_list))
testError_list_var <- testError_list[[index]]

testingErrorArranged <- lapply(seq_along(testError_list_var), function(i) {
  lapply(seq_along(testError_list_var[[i]]), function(j) {
    df <- as.data.frame(matrix(testError_list_var[[i]][[j]],nrow = 1))
    df$Percent <- as.numeric(i)
    df$Dataset <- as.numeric(j)
    df
  })
})

testingError_df <- do.call(rbind, do.call(c, testingErrorArranged))
testingError_df$variable <- var_names[i]

transformed_test <- testingError_df %>%
  pivot_longer(cols = c("V1","V2","V3","V4","V5"), names_to = "Imputation", values_to = "Error") %>%
  mutate(ErrorType = "TestingError")

transformed_test$Imputation <- gsub("^V", "", transformed_test$Imputation )
transformed_test$Percent<-prop_values[match(transformed_test$Percent , 1:length(prop_values))]



combined_Error<-rbind(transformed_train, transformed_test)


combined_Error$Percent<-paste(combined_Error$Percent, "% missingness", sep = "")
combined_Error$Dataset <- factor(combined_Error$Dataset, levels = unique(combined_Error$Dataset))


color_palette <- c("TrainingError" = "blue", "TestingError" = "red")

# Plot using ggplot2
p[[i]] <- ggplot(combined_Error, aes(x = Dataset, y = Error, color = ErrorType, group = interaction(ErrorType, Imputation))) +
  geom_line() +
  facet_wrap(~ Percent) +
  scale_color_manual(values = color_palette) +
  labs(x = "Dataset", y = "Error", color = "Error Type") + ggtitle(paste("Variable: ", unique(combined_Error$variable))) +
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5)) 
}

pdf("Training_test_XGBDefault_subsample70.pdf")
p
dev.off()





########################################################################################################################
#######################################           XGB - No subsample.     ##############################################
########################################################################################################################


########################    4 - XGBoost ####################################

t1_xgbnosub <- Sys.time()

impute_MAR_nosub <- missing_MAR_list %>%
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      mice::mice(mat, 
                 m = m, method = "xgb",
                 maxit = maxit,xgb.params=NULL,
                 print = FALSE)
    })
  }, .options = furrr_options(seed = 123))

totalTime_xgbnosub <- Sys.time() - t1_xgbnosub

eval_xgb_nosub <- impute_MAR_nosub %>% 
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z)) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(eval_xgb_nosub,file =  "MI_xgb_nosub.RData")


########################################################################

imputed_data <- lapply(impute_MAR_nosub, function(inner_list) lapply(inner_list, '[[','imp'))
var_names <-  colnames(simdata[[1]])

# Calculation of testing Error
testError_list = list()

for(k in var_names){
  print(k)
  missing_index <- missing_MAR_list %>% 
    map(~map(.x, ~ which(is.na(.x[, k]))))  #finds index of missing data
  actual_data <- lapply(missing_index, function(x) {
    lapply(1:length(x), function(i) {simdata[[i]][x[[i]],k]})}) # retrieves actual data which was amputed as missing
  
  imputed_var <- lapply(imputed_data, function(inner_list) #retrieves imputed data for each missing variable
    lapply(inner_list, '[[',k))
  
  testError_list[[k]] <- lapply(seq_along(imputed_var), function(i) {
    lapply(seq_along(imputed_var[[i]]), function(j) { 
      apply(imputed_var[[i]][[j]],2, 
            function(col) {as.matrix(colMeans(actual_data[[i]][[j]] - col))})})
  })
}



########################################################################
##########             Training and test arrangement       #############
########################################################################


train_error <- lapply(impute_MAR_nosub, function(inner_list) lapply(inner_list, '[[','trainError'))
train_error_f1 <- lapply(train_error, function(x){ lapply(seq_along(x), function(y) {x[[y]][!sapply(x[[y]], is.null)]})}) # removes NULL element from training data
p = list()

for(i in 1:length(var_names)){
  
  train_error_var <- lapply(train_error_f1, function(inner_list) lapply(inner_list, '[[',var_names[i]))
  trainingErrorArranged <- lapply(seq_along(train_error_var), function(i) {
    lapply(seq_along(train_error_var[[i]]), function(j) {
      df <- train_error_var[[i]][[j]]
      df$Percent <- prop_values[i]
      df$Dataset <- j
      df
    })
  })
  
  trainingError_df <- do.call(rbind, do.call(c, trainingErrorArranged))
  trainingError_df$variable <- var_names[i]
  
  transformed_train <- trainingError_df %>%
    pivot_longer(cols = c("V1","V2","V3","V4","V5"), names_to = "Imputation", values_to = "Error") %>%
    mutate(ErrorType = "TrainingError")
  
  transformed_train$Imputation <- gsub("^V", "", transformed_train$Imputation )
  
  
  # test error arrangement
  
  index <- match(var_names[i], names(testError_list))
  testError_list_var <- testError_list[[index]]
  
  testingErrorArranged <- lapply(seq_along(testError_list_var), function(i) {
    lapply(seq_along(testError_list_var[[i]]), function(j) {
      df <- as.data.frame(matrix(testError_list_var[[i]][[j]],nrow = 1))
      df$Percent <- as.numeric(i)
      df$Dataset <- as.numeric(j)
      df
    })
  })
  
  testingError_df <- do.call(rbind, do.call(c, testingErrorArranged))
  testingError_df$variable <- var_names[i]
  
  transformed_test <- testingError_df %>%
    pivot_longer(cols = c("V1","V2","V3","V4","V5"), names_to = "Imputation", values_to = "Error") %>%
    mutate(ErrorType = "TestingError")
  
  transformed_test$Imputation <- gsub("^V", "", transformed_test$Imputation )
  transformed_test$Percent<-prop_values[match(transformed_test$Percent , 1:length(prop_values))]
  
  
  
  combined_Error<-rbind(transformed_train, transformed_test)
  
  
  combined_Error$Percent<-paste(combined_Error$Percent, "% missingness", sep = "")
  combined_Error$Dataset <- factor(combined_Error$Dataset, levels = unique(combined_Error$Dataset))
  
  
  color_palette <- c("TrainingError" = "blue", "TestingError" = "red")
  
  # Plot using ggplot2
  p[[i]] <- ggplot(combined_Error, aes(x = Dataset, y = Error, color = ErrorType, group = interaction(ErrorType, Imputation))) +
    geom_line() +
    facet_wrap(~ Percent) +
    scale_color_manual(values = color_palette) +
    labs(x = "Dataset", y = "Error", color = "Error Type") + ggtitle(paste("Variable: ", unique(combined_Error$variable))) +
    theme_minimal()+theme(plot.title = element_text(hjust = 0.5)) 
}

pdf("Training_test_XGBDefault_nosub.pdf")
p
dev.off()

############################################################################################################################


########################    4 - XGBoost - Random parameter ####################################

t1_xgb_rp <- Sys.time()


random_param_set <- missing_MAR_list %>%
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      xgb_param_calc(mat,response = NULL, select_features=NULL, num_cores = 6)
    })
  }, .options = furrr_options(seed = 123))


mice_results_xgb_param_random_maxit5 <- map2(missing_MAR_list, random_param_set, function(data_inner, params_inner) {
  map2(data_inner, params_inner, function(data_single, params_single) {
  mice(data_single, m = m, method = "xgb", maxit = maxit,xgb.params =  params_single$parameter, print = FALSE)
  })})
totalTime_xgbnosub <- Sys.time() - t1_xgb_rp


eval_xgb_param_random_maxit5 <- mice_results_xgb_param_random_maxit5 %>% 
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z)) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})
  
save(eval_xgb_param_random_maxit5,file =  "MI_xgb_param_random_maxit5.RData")

  
########################################################################

imputed_data <- lapply(mice_results_xgb_param_random_maxit5, function(inner_list) lapply(inner_list, '[[','imp'))
var_names <-  colnames(simdata[[1]])

# Calculation of testing Error
testError_list = list()

for(k in var_names){
  print(k)
  missing_index <- missing_MAR_list %>% 
    map(~map(.x, ~ which(is.na(.x[, k]))))  #finds index of missing data
  actual_data <- lapply(missing_index, function(x) {
    lapply(1:length(x), function(i) {simdata[[i]][x[[i]],k]})}) # retrieves actual data which was amputed as missing
  
  imputed_var <- lapply(imputed_data, function(inner_list) #retrieves imputed data for each missing variable
    lapply(inner_list, '[[',k))
  
  testError_list[[k]] <- lapply(seq_along(imputed_var), function(i) {
    lapply(seq_along(imputed_var[[i]]), function(j) { 
      apply(imputed_var[[i]][[j]],2, 
            function(col) {as.matrix(colMeans(actual_data[[i]][[j]] - col))})})
  })
}



########################################################################
##########             Training and test arrangement       #############
########################################################################


train_error <- lapply(mice_results_xgb_param_random_maxit5, function(inner_list) lapply(inner_list, '[[','trainError'))
train_error_f1 <- lapply(train_error, function(x){ lapply(seq_along(x), function(y) {x[[y]][!sapply(x[[y]], is.null)]})}) # removes NULL element from training data
p = list()

for(i in 1:length(var_names)){
  
  train_error_var <- lapply(train_error_f1, function(inner_list) lapply(inner_list, '[[',var_names[i]))
  trainingErrorArranged <- lapply(seq_along(train_error_var), function(i) {
    lapply(seq_along(train_error_var[[i]]), function(j) {
      df <- train_error_var[[i]][[j]]
      df$Percent <- prop_values[i]
      df$Dataset <- j
      df
    })
  })
  
  trainingError_df <- do.call(rbind, do.call(c, trainingErrorArranged))
  trainingError_df$variable <- var_names[i]
  
  transformed_train <- trainingError_df %>%
    pivot_longer(cols = c("V1","V2","V3","V4","V5"), names_to = "Imputation", values_to = "Error") %>%
    mutate(ErrorType = "TrainingError")
  
  transformed_train$Imputation <- gsub("^V", "", transformed_train$Imputation )
  
  
  # test error arrangement
  
  index <- match(var_names[i], names(testError_list))
  testError_list_var <- testError_list[[index]]
  
  testingErrorArranged <- lapply(seq_along(testError_list_var), function(i) {
    lapply(seq_along(testError_list_var[[i]]), function(j) {
      df <- as.data.frame(matrix(testError_list_var[[i]][[j]],nrow = 1))
      df$Percent <- as.numeric(i)
      df$Dataset <- as.numeric(j)
      df
    })
  })
  
  testingError_df <- do.call(rbind, do.call(c, testingErrorArranged))
  testingError_df$variable <- var_names[i]
  
  transformed_test <- testingError_df %>%
    pivot_longer(cols = c("V1","V2","V3","V4","V5"), names_to = "Imputation", values_to = "Error") %>%
    mutate(ErrorType = "TestingError")
  
  transformed_test$Imputation <- gsub("^V", "", transformed_test$Imputation )
  transformed_test$Percent<-prop_values[match(transformed_test$Percent , 1:length(prop_values))]
  
  
  
  combined_Error<-rbind(transformed_train, transformed_test)
  
  
  combined_Error$Percent<-paste(combined_Error$Percent, "% missingness", sep = "")
  combined_Error$Dataset <- factor(combined_Error$Dataset, levels = unique(combined_Error$Dataset))
  
  
  color_palette <- c("TrainingError" = "blue", "TestingError" = "red")
  
  # Plot using ggplot2
  p[[i]] <- ggplot(combined_Error, aes(x = Dataset, y = Error, color = ErrorType, group = interaction(ErrorType, Imputation))) +
    geom_line() +
    facet_wrap(~ Percent) +
    scale_color_manual(values = color_palette) +
    labs(x = "Dataset", y = "Error", color = "Error Type") + ggtitle(paste("Variable: ", unique(combined_Error$variable))) +
    theme_minimal()+theme(plot.title = element_text(hjust = 0.5)) 
}

pdf("Training_test_XGB_param_random_maxit5.pdf")
p
dev.off()


########################    4 - XGBoost - All parameter ####################################

t1_xgb_ap <- Sys.time()


all_param_set <- missing_MAR_list %>%
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      xgb_param_calc(mat,response = "all", select_features=NULL, num_cores = 6)
    })
  }, .options = furrr_options(seed = 123))
time_param <- Sys.time() - t1_xgb_ap


t1_xgb_ap_mice <- Sys.time()
mice_results_xgb_param_all_maxit5 <- map2(missing_MAR_list, all_param_set, function(data_inner, params_inner) {
  map2(data_inner, params_inner, function(data_single, params_single) {
    mice(data_single, m = m, method = "xgb", maxit = maxit,xgb.params =  params_single$parameter, print = FALSE)
  })})
totalTime_xgb_ap_mice <- Sys.time() - t1_xgb_ap_mice


save(mice_results_xgb_param_all_maxit5,file =  "MI_xgb_param_all_maxit5.RData")



eval_xgb_param_all_maxit5 <- mice_results_xgb_param_all_maxit5 %>% 
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z)) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})


########################################################################

imputed_data <- lapply(mice_results_xgb_param_all_maxit5, function(inner_list) lapply(inner_list, '[[','imp'))
var_names <-  colnames(simdata[[1]])

# Calculation of testing Error
testError_list = list()

for(k in var_names){
  print(k)
  missing_index <- missing_MAR_list %>% 
    map(~map(.x, ~ which(is.na(.x[, k]))))  #finds index of missing data
  actual_data <- lapply(missing_index, function(x) {
    lapply(1:length(x), function(i) {simdata[[i]][x[[i]],k]})}) # retrieves actual data which was amputed as missing
  
  imputed_var <- lapply(imputed_data, function(inner_list) #retrieves imputed data for each missing variable
    lapply(inner_list, '[[',k))
  
  testError_list[[k]] <- lapply(seq_along(imputed_var), function(i) {
    lapply(seq_along(imputed_var[[i]]), function(j) { 
      apply(imputed_var[[i]][[j]],2, 
            function(col) {as.matrix(colMeans(actual_data[[i]][[j]] - col))})})
  })
}



########################################################################
##########             Training and test arrangement       #############
########################################################################


train_error <- lapply(mice_results_xgb_param_random_maxit5, function(inner_list) lapply(inner_list, '[[','trainError'))
train_error_f1 <- lapply(train_error, function(x){ lapply(seq_along(x), function(y) {x[[y]][!sapply(x[[y]], is.null)]})}) # removes NULL element from training data
p = list()

for(i in 1:length(var_names)){
  
  train_error_var <- lapply(train_error_f1, function(inner_list) lapply(inner_list, '[[',var_names[i]))
  trainingErrorArranged <- lapply(seq_along(train_error_var), function(i) {
    lapply(seq_along(train_error_var[[i]]), function(j) {
      df <- train_error_var[[i]][[j]]
      df$Percent <- prop_values[i]
      df$Dataset <- j
      df
    })
  })
  
  trainingError_df <- do.call(rbind, do.call(c, trainingErrorArranged))
  trainingError_df$variable <- var_names[i]
  
  transformed_train <- trainingError_df %>%
    pivot_longer(cols = c("V1","V2","V3","V4","V5"), names_to = "Imputation", values_to = "Error") %>%
    mutate(ErrorType = "TrainingError")
  
  transformed_train$Imputation <- gsub("^V", "", transformed_train$Imputation )
  
  
  # test error arrangement
  
  index <- match(var_names[i], names(testError_list))
  testError_list_var <- testError_list[[index]]
  
  testingErrorArranged <- lapply(seq_along(testError_list_var), function(i) {
    lapply(seq_along(testError_list_var[[i]]), function(j) {
      df <- as.data.frame(matrix(testError_list_var[[i]][[j]],nrow = 1))
      df$Percent <- as.numeric(i)
      df$Dataset <- as.numeric(j)
      df
    })
  })
  
  testingError_df <- do.call(rbind, do.call(c, testingErrorArranged))
  testingError_df$variable <- var_names[i]
  
  transformed_test <- testingError_df %>%
    pivot_longer(cols = c("V1","V2","V3","V4","V5"), names_to = "Imputation", values_to = "Error") %>%
    mutate(ErrorType = "TestingError")
  
  transformed_test$Imputation <- gsub("^V", "", transformed_test$Imputation )
  transformed_test$Percent<-prop_values[match(transformed_test$Percent , 1:length(prop_values))]
  
  
  
  combined_Error<-rbind(transformed_train, transformed_test)
  
  
  combined_Error$Percent<-paste(combined_Error$Percent, "% missingness", sep = "")
  combined_Error$Dataset <- factor(combined_Error$Dataset, levels = unique(combined_Error$Dataset))
  
  
  color_palette <- c("TrainingError" = "blue", "TestingError" = "red")
  
  # Plot using ggplot2
  p[[i]] <- ggplot(combined_Error, aes(x = Dataset, y = Error, color = ErrorType, group = interaction(ErrorType, Imputation))) +
    geom_line() +
    facet_wrap(~ Percent) +
    scale_color_manual(values = color_palette) +
    labs(x = "Dataset", y = "Error", color = "Error Type") + ggtitle(paste("Variable: ", unique(combined_Error$variable))) +
    theme_minimal()+theme(plot.title = element_text(hjust = 0.5)) 
}

pdf("Training_test_XGB_param_all_maxit5.pdf")
p
dev.off()


###################################################################################################################################

t1_mixgb <- Sys.time()
impute_MAR_mixgb <- missing_MAR_list %>%
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {
      mixgb(mat, 
                 m = m,
                 maxit = maxit,
                 print = FALSE)
    })
  }, .options = furrr_options(seed = 123))

totalTime_mixgb <-  Sys.time() - t1_cart


eval_mixgb <- impute_MAR_mixgb %>% 
  furrr::future_map(function(mat_list) {
    lapply(mat_list, function(mat) {mat %>%
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z))  %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(eval_mixgb, file = "MI_mixgb.RData")
















cor_list <- map(simdata, ~ cor(.x)) %>% Reduce("+", .) / length(simdata)

cors_all <- cbind( 
  cor_withY = cor_list[c("y","x","z"),"y"],
         cor_withX =cor_list[c("y","x","z"),"x"],
         cor_withZ = cor_list[c("y","x","z"),"z"]
         )


saveRDS(cors_all, "Correlations_XYZ.rds")

