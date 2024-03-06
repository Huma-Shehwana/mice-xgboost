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
    
    file_path <- "Param.RData"
    
    if (file.exists(file_path)) {
      load(file_path)
      param_data <- rbind(param_data, param)
    } else {
      param_data <- param
    }
    save(param_data, file = file_path)
  }
  return(param)
}



#####################################################################################################
###################     Simulation and evaluation function      #####################################
#####################################################################################################



sim_eval<-function(data_mar, m, maxit, method,param){
  
  
  one_data_run <- function(x, param, m, maxit, method ) {
    
      x %>%
        mice(m = m,
             maxit = maxit,
             method = method,
             print = F, xgb.params = param)
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

