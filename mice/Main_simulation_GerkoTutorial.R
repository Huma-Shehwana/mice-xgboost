rm(list = ls())


#library(mice)     # for imputation and amputation
library(purrr)    # for functional programming
library(furrr)    # for functional futures
library(mvtnorm)  # for multivariate normal data
library(magrittr) # for pipes
library(dplyr)    # for data manipulation
library(tibble)   # for tibbles
library(mixgb)
library(microbenchmark)


seed = 123
set.seed(seed)     




#######################################################################################
####################    Parameter Estimation function    ######################
#######################################################################################



paramEstimation<-function(mbased_MAR, method){
  if(method !="xgbParam"){
    param <- replicate(length(mbased_MAR), NULL)
  } else if(method=="xgbParam"){
    param <- 
      mbased_MAR %>%
      furrr::future_map(function(x) {
        x %>% 
          xgb_bayes_opt_completeCases_oneVariable(response = NULL, select_features=NULL)
      }, .options = furrr_options(seed = seed))
  }
  return(param)
}



#####################################################################################################
###################  Simulation and evaluation function  #####################################
#####################################################################################################
sim_eval<-function(data_mar, m, maxit, method,param){
  
  
  one_data_run <- function(x, param, m, maxit, method ) {
    
    if(is.null(param)) {
      x %>%
        mice(m = m,
             maxit = maxit,
             method = method,
             print = F)
    } else {
      x %>%
        mice(m = m,
             maxit = maxit,
             method = method,
             print = F, xgb.params = param)
    }
  }
  
  print(paste("Imputation using", method))
  
  if(method!="mixgb"){
  result_list <- furrr::future_map2(data_mar, param, ~one_data_run(.x, .y, m, maxit, method), .options = furrr_options(seed = seed))
  
  eval <- result_list %>% 
    map(~.x %>% # for every simulated multiple imputation....
          complete("all") %>% # create a list of completed data sets
          map(~.x %$% # for every completed data set....
                lm(y ~ x + z)) %>% # fit linear model
          pool() %>%  # pool coefficients
          summary(conf.int = TRUE) %>% # summary of coefficients
          mutate(true = c(0, 6, 3), # add true
                 cov = conf.low < true & true < conf.high, # coverage
                 bias = estimate - true,
                 width = conf.high - conf.low) %>% # bias
          column_to_rownames("term")) %>% # `term` as rownames
    Reduce("+", .) / length(data_mar) # add all and divide by length 
  } else if (method=="mixgb") {
    
    mixgb_results<-mixgb_results<-data_mar %>%
                  furrr::future_map(function(x) { x %>%
                  mixgb(m=m,maxit=1)}, .options = furrr_options(seed = seed))
    
    eval <- mixgb_results %>% 
      map(~.x %>% # for every simulated multiple imputation....
            map(~.x %$% lm(y ~ x + z)) %>% # fit linear model
            pool() %>%  # pool coefficients
            summary(conf.int = TRUE) %>% # summary of coefficients
            mutate(true = c(0, 6, 3), # add true
                   cov = conf.low < true & true < conf.high, # coverage
                   bias = estimate - true,
                   width = conf.high - conf.low) %>% # bias
            column_to_rownames("term")) %>% # `term` as rownames
      Reduce("+", .) / length(data_mar) # add all and divide by length 
}
  
  return(eval)
  
}



##### Step 1: Data generation


sigma <- matrix(data = c(1, 0.7, 0.7, 1), 
                ncol = 2)
simdata <- replicate(n = 2, 
                     expr = mvtnorm::rmvnorm(n = 1000, 
                                             mean = c(8, 3), 
                                             sigma = sigma) %>% 
                       as_tibble() %>% # make into a tibble
                       rename(x = V1, z = V2) %>% # rename columns
                       mutate(y = 6 * x + 3 * z + rnorm(50)), # add y
                     simplify = FALSE) # keep as list of generated sets


simdata %>% 
  map(~.x %$% # for every simulated set in simdata....
        lm(y ~ x + z) %>% # fit linear model
        coefficients) %>% # extract coefficients
  Reduce("+", .) / length(simdata) # add all and divide by length (= average)


######## Step 2: Simulations    

#prop = 0.5
prop = c(0.2,0.3)
         #0.4,0.5)

#methods<- c("norm")
methods<- c("norm","rf","cart","xgb","xgbParam", "mixgb")
m=5
maxit=5
OutputFile <- "Simulations_1000_donors"

sim_res <- data.frame(matrix(ncol = length(methods)*11, nrow = 0))





prop_results <- lapply(prop_list, function(prop_value) {
  mbased_MAR <- simdata %>%
    furrr::future_map(function(x) {
      x %>% 
        ampute(prop = prop_value, 
               mech = "MAR", type = "RIGHT") %>% .$amp 
    }, .options = furrr_options(seed = seed))
  
  

totalTime <- vector()

for(j in 1:length(prop)){

mbased_MAR <- 
  simdata %>%
  furrr::future_map(function(x) {
    x %>% 
      ampute(prop = prop[j], 
             mech = "MAR", type = "RIGHT") %>% .$amp 
  }, .options = furrr_options(seed = seed))

NAs_in_data <- as_tibble(flatten_dbl(mbased_MAR %>% map(function(x){sum(is.na(x))})))
print("Number of NAs in datasets")
print(NAs_in_data)


eval = list()

print(sprintf("Starting simulation for %s %% missing data", prop[j]))

for(i in 1:length(methods)){
  time_taken <- system.time( eval[[i]] <- sim_eval(mbased_MAR, m = m, 
                                       maxit = maxit, method = methods[i],
                                       param=paramEstimation(mbased_MAR, methods[i]))
                             )
  totalTime <- c(totalTime, time_taken[3])
}

eval_cbind<-bind_cols(eval)

if(j==1) {
names_cbind<-c(sapply(1:length(eval), function(x) {paste(methods[x],colnames(eval[[x]]), sep='_')}))
colnames(sim_res) <- names_cbind
}

colnames(eval_cbind) <- colnames(sim_res)
sim_res <- rbind(sim_res, eval_cbind)
write.table(sim_res, sprintf("%s_prop%s.tsv", OutputFile, prop[j]), sep = '\t')
}      


cor_list <- map(simdata, ~ cor(.x)) %>% Reduce("+", .) / length(simdata)

sim_res <- sim_res %>% 
  mutate(cor_withY = as.vector(replicate(length(prop), cor_list[c("y","x","z"),"y"])),
         cor_withX = as.vector(replicate(length(prop), cor_list[c("y","x","z"),"x"])),
         cor_withZ = as.vector(replicate(length(prop), cor_list[c("y","x","z"),"z"])),
         Missing_prop = as.vector(rep(prop,each = 3))
         )


write.table(sim_res, sprintf("%s.tsv",OutputFile), sep='\t')








 