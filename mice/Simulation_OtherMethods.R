rm(list = ls())

source("Simulation_function.R")
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
Num_ds <- 1
prop_list = c(20,30, 40,50)
methods<- c("norm","rf","cart","xgb","mixgb","xgbParam")
m=5
maxit=5

##### Step 1: Data generation


sigma <- matrix(data = c(1, 0.7, 0.7, 1), 
                ncol = 2)

# Defining a data matrix
simdata <- replicate(n = Num_ds, 
                     expr = mvtnorm::rmvnorm(n = 1000, 
                                             mean = c(8, 3), 
                                             sigma = sigma) %>% 
                       as_tibble() %>% # make into a tibble
                       rename(x = V1, z = V2) %>% # rename columns
                       mutate(y = 6 * x + 3 * z + rnorm(1000)), # add y
                     simplify = FALSE) # keep as list of generated sets

# regression to get estimates
simdata %>% 
  map(~.x %$% # for every simulated set in simdata....
        lm(y ~ x + z) %>% # fit linear model
        coefficients) %>% # extract coefficients
  Reduce("+", .) / length(simdata) # add all and divide by length (= average)


######## Step 2: Simulations    


OutputFile <- "Simulations_1000_donors"

prop_results <- lapply(prop_list, function(prop_value) {
  
  
        mbased_MAR <- simdata %>%
                      furrr::future_map(function(x) {
                      x %>% 
                      ampute(prop = prop_value/100, 
                      mech = "MAR", type = "RIGHT") %>% .$amp 
                      }, .options = furrr_options(seed = seed))
  
  



        NAs_in_data <- as_tibble(flatten_dbl(mbased_MAR %>% map(function(x){sum(is.na(x))})))
        print("Number of NAs in datasets")
        print(NAs_in_data)
        print(sprintf("Starting simulation for %s %% missing data", prop_value))


        
        eval <- lapply(methods, function(method) {
  
                      t1 <- Sys.time()
                      sim_eval_result <- sim_eval(mbased_MAR, m = m, 
                                                  maxit = maxit, method = method,
                                                  param = paramEstimation(mbased_MAR, method))
                      
                      totalTime <- t1 - Sys.time()
                      results <- list(method = method, eval_result = sim_eval_result, 
                                      time_taken = totalTime, missing_prop = prop_value)
                      
                      save(results, file = sprintf("sim_%s_%s_%s.RData", Num_ds,prop_value, method))
      })


})

    


cor_list <- map(simdata, ~ cor(.x)) %>% Reduce("+", .) / length(simdata)

cors_all <- cbind( 
  cor_withY = cor_list[c("y","x","z"),"y"],
         cor_withX =cor_list[c("y","x","z"),"x"],
         cor_withZ = cor_list[c("y","x","z"),"z"]
         )


saveRDS(cors_all, "Correlations_XYZ.rds")


start_points <- seq(21, 500, by = 25)
end_points <- seq(25, 500, by = 25)
indices <- sapply(start_points, function(x) {seq(x,length.out = 5)})

load("XGBParam_TrainingPredictionBias_x.RData")
Training_prediction_x<- apply(indices, 2, function(x) {mean(existing_data[x,])})
load("XGBParam_TrainingPredictionBias_y.RData")
Training_prediction_y<- apply(indices, 2, function(x) {mean(existing_data[x,])})
load("XGBParam_TrainingPredictionBias_z.RData")
Training_prediction_z<- apply(indices, 2, function(x) {mean(existing_data[x,])})



Training_prediction <- as.data.frame(cbind(Training_prediction_x, Training_prediction_y, Training_prediction_z))
Training_prediction$Datasets <- rep(1:5, 4)
Training_prediction$PropNA <- rep(c(20,30,40,50),each = 5)
save(Training_prediction, file = "TrainingPredictions.RData")



start_points <- seq(21, 500, by = 25)
end_points <- seq(25, 500, by = 25)
indices <- sapply(start_points, function(x) {seq(x,length.out = 5)})

load("TrainingPredictionBias_x.RData")
existing_data <- existing_data[975:1475,,drop=FALSE]
Training_prediction_x2<- apply(indices, 2, function(x) {mean(existing_data[x,])})


load("TrainingPredictionBias_y.RData")
existing_data <- existing_data[975:1475,,drop=FALSE]
Training_prediction_y2<- apply(indices, 2, function(x) {mean(existing_data[x,])})


load("TrainingPredictionBias_z.RData")
existing_data <- existing_data[975:1475,,drop=FALSE]
Training_prediction_z2<- apply(indices, 2, function(x) {mean(existing_data[x,])})
Training_prediction_2 <- as.data.frame(cbind(Training_prediction_x2, Training_prediction_y2, Training_prediction_z2))


colnames(Training_prediction_2) <- paste("DefaultParam", colnames(Training_prediction_2), sep = "_")


TrainPredict_combined <- cbind(Training_prediction_2,Training_prediction)
save(TrainPredict_combined, file = "AllTrainingPredictions.RData")





 