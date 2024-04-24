#' Description of your function.
#'
#'@aliases mice.impute.xgb
#' @param y: Vector of values of response variable
#' @param ry: Vector of true/false representing observed values in target variable. TRUE representing respective value was observed and FALSE represents missing value.
#' @param wy: Vector of true/false representing missing values in target variable. TRUE representing respective value was missing and FALSE represents observed value.
#' @param x: dataframe. Data of all variables except target variable.
#' @param xgb.params: list of XGBoost parameters obtained from xgb_param_calc function. User can also manually provide a list of parameters.
#' @param nthreads: Number of threads. By default, function uses all available - 1 cores.
#' @param match.type: Possible value can be "predicted.observed", "original.observed" or "predicted"
#' @return Imputed values for target values and training error (actual - predicted)
#' 
#' XGBoost uses a set of hyperparameters for tuning the algorithm. User has three options to select these parameters
#' 1) Using default parameters (adapted from Deng et al., 2023).
#' 2) Run xgb_param_calc function to let function randomly select one variable as a target and tune the parameters using only one target variable.
#' 3) Run xgb_param_calc function to iteratively select best set of hyperparameters for all variables.
#' 
#'mice.impute.xgb
#' @examples
#' mice(nhanes, method = "xgb", param = NULL)
#'
#' @export


mice.impute.xgb<- function(y, ry, x, wy=NULL,xgb.params = list(), nthreads = -1,match.type = "predicted",...) {
  
  ###########################################################################################################
  ############################          Parameters
  ###########################################################################################################
  
  
  # hardcoding general parameters
  nrounds = 100 # number of trees
  nthread <- nthreads
  early_stopping_rounds = 30 # if algorithm does not improve evaluation metric for 30 cycles, it can stop
  print_every_n = 50 # to make console less crowded but to also check if algorithm is also working.
  verbose = 0
  
  #if no hyperparameters are provided, algorithm uses default set of parameters adapted from Deng etal, 2023
  # If only one set of hyperparameters are provided it is sued to impute all variables, otherwise, hyperparameters specific to each variable is sued to train the model
  
  if(is.null(xgb.params)){
    params <- list(tree_method = "hist", eta = 0.3, gamma = 0, max_depth = 3, 
                   min_child_weight = 1, subsample = 1, sampling_method = "uniform", 
                   colsample_bytree = 1, colsample_bylevel = 1, colsample_bynode = 1, 
                   lambda = 1, alpha = 0,max_bin = 256, num_parallel_tree = 1, nthread = -1) 
  } else if(all(names(xgb.params) %in% colnames(x))) {
    
    param_oneVariable <- xgb.params[[setdiff(names(xgb.params), colnames(x))]]
    params <- c(param_oneVariable, list(tree_method = "hist", num_parallel_tree = 8))
    
  } else {
    
    params <- c(xgb.params, list(device = "cpu",tree_method = "hist", num_parallel_tree = 1, 
                                 nthread = nthreads))
  }
  
  ############################################################################################################
  
  var.type <- switch(class(y)[1],
                     numeric = "numeric",
                     integer = "numeric",
                     factor = ifelse(nlevels(y) == 2, "binary", "multiclass"),
                     ordered = ifelse(nlevels(y) == 2, "binary", "multiclass"),
                     stop(paste0("Unsupported data type: ", class(y)))
  )
  
  obj.type <- switch(var.type,
                     numeric = "reg:squarederror",
                     binary =  "predicted","binary:logistic",
                     multiclass = ifelse(match.type=="predicted", "multi:softmax", "multi:softprob"),
                     stop(paste0("Unsupported data type: ", class(y)))
  )
  
  eval.metric <- switch(var.type,
                        numeric = "MAPE",
                        binary = "logloss",
                        multiclass = "mlogloss",
                        stop(paste0("Unsupported data type: ", class(y)))
  )
  
  ############################################################################################################
  
  if (is.null(wy)) 
    wy <- !ry
  
  nmis <- sum(wy)                 # number of missing entries for target variable
  xobs <- x[ry, , drop = FALSE]   # All variable data (-target variable), where target variable was observed
  xmis <- x[wy, , drop = FALSE]   # All variable data (-target variable), where target variable was missing
  yobs <- y[ry]                   #observed target values
  
  
  if(var.type != "numeric"){
    obs.y <- as.integer(obs.y) - 1
    bin.t <- sort(table(obs.y))
    
    
    if (is.na(bin.t[2])) {
      imp <- levels(obs.y)[as.integer(names(bin.t[1])) + 1]
      msg <- paste("The binary variable", var, "in the data only have single class. Imputation models can't be built.")
      stop(msg)
    }
  }
  
  ######################################################################################
  
  dobs <- xgboost::xgb.DMatrix(data = xobs, label = yobs, nthread = nthreads) #Defining xgboost object
  dmis <- xgboost::xgb.DMatrix(data = xmis, nthread = nthreads) #Defining xgboost object
  
  watchlist <- list(train = dobs)
  
  if(var.type!="multiclasss"){
    xgb.fit <- xgboost::xgb.train(
      data = dobs, objective = obj.type, watchlist = watchlist,
      params = params, nrounds = nrounds, early_stopping_rounds = early_stopping_rounds, 
      print_every_n = print_every_n, verbose = verbose
    ) # training the model on observed data using parameters set earlier on this page
  } else if(var.type=="multiclass") {
    N.class <- length(levels(y))
    xgb.fit <- xgboost::xgb.train(
      data = dobs,num_class = N.class, objective = obj.type, watchlist = watchlist,
      params = params, nrounds = nrounds, early_stopping_rounds = early_stopping_rounds, 
      print_every_n = print_every_n, verbose = verbose
    ) # training the model on observed data using parameters set earlier on this page
  }
  
  yhatmis <- predict(xgb.fit, dmis) #obtaining the prediction for missing y variable
  yhatobs <- predict(xgb.fit, dobs) # obtaining the prediction for observed y variable
  
  
  if(var.type == "numeric"){
    if(match.type == "predicted.observed"){
      idx <- matchindex(d = yhatobs, t = yhatmis, k = k)
      imp <-yobs[idx]
    } else if(match.type == "predicted" | is.null(match.type)) {
      imp <- yhatmis
    } else if(match.type == "original.observed"){
      idx <- matchindex(d = yobs, t = yhatmis, k = k)
      imp <-yobs[idx]
    }
    
  }
  
  if(var.type == "binary"){
    if(match.type == "predicted"){
      yhatmis <- ifelse(yhatmis >= 0.5, 1, 0)
      yhatmis <- levels(y)[yhatmis + 1]
      imp <- yhatmis
    } else if(match.type == "predicted.observed"){
      idx <- matchindex(d = yhatobs, t = yhatmis, k = k)
      imp <-yobs[idx]
    } else if(match.type == "original.observed"){
      idx <- matchindex(d = yobs, t = yhatmis, k = k)
      imp <-yobs[idx]
    }
    
    
    
    
    
    
    
  } else if(var.type=="multiclass"){
    imp <- levels(sorted.dt[[var]])[yhatmis + 1]
  }
  bias <- mean(yobs - yhatobs)
  
  error = bias
  
  return(list(imp, error))
  
  
  
  
  
}



