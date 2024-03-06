#' Description of your function.
#'
#'@aliases mice.impute.xgb
#' @param arg1 Description of arg1.
#' @param arg2 Description of arg2.
#' @return Description of the return value.
#'
#' @examples
#' # Example usage of your_function
#' your_function(arg1_value, arg2_value)
#'
#' @export


mice.impute.xgb<- function(y, ry, x, wy = NULL,xgb.params = list(), ...) {
  
  
  
############################################################################################################
############################################################################################################
  if(is.null(xgb.params)){
    params <- list(device = "cpu",tree_method = "hist", 
                       eta = 0.3, gamma = 0, max_depth = 3, min_child_weight = 1, 
                       subsample = 1, sampling_method = "uniform", 
                      colsample_bytree = 1, colsample_bylevel = 1, 
                      colsample_bynode = 1, lambda = 1, alpha = 0,
                      max_bin = 256, num_parallel_tree = 1, nthread = -1) 
  } else if(all(names(xgb.params) %in% colnames(x))) {
    browser()
    param_oneVariable <- xgb.params[[setdiff(names(xgb.params), colnames(x))]]
    params <- c(param_oneVariable, list(device = "cpu",tree_method = "hist", num_parallel_tree = 1, nthread = -1))
  } else {
    params <- c(xgb.params, list(device = "cpu",tree_method = "hist", num_parallel_tree = 1, nthread = -1))
  }
  

############################################################################################################
############################################################################################################
  
  
  nrounds = 200
  nthread <- params$nthread
  early_stopping_rounds = 20
  print_every_n = 30
  verbose = 0
  
############################################################################################################
############################################################################################################

  if (is.null(wy)) wy <- !ry
  
  nmis <- sum(wy)
  xobs <- x[ry, , drop = FALSE]
  xmis <- x[wy, , drop = FALSE]
  yobs <- y[ry]

  
  dobs <- xgboost::xgb.DMatrix(data = xobs, label = yobs, nthread = nthread)
  dmis <- xgboost::xgb.DMatrix(data = xmis, nthread = nthread)
  
  watchlist <- list(train = dobs)
  obj.type <- "reg:squarederror"
  
  xgb.fit <- xgboost::xgb.train(
                      data = dobs, objective = obj.type, watchlist = watchlist,
                      params = params, nrounds = nrounds, early_stopping_rounds = early_stopping_rounds, 
                      print_every_n = print_every_n, verbose = verbose
                      )
 
   yhatmis <- predict(xgb.fit, dmis)
   yhatobs <- predict(xgb.fit, dobs)
   idx <- matchindex(d = yhatobs, t = yhatmis, k = 5)
   imp <-yobs[idx]
   
   bias <- mean(yhatobs - yobs)
   
   error = bias
   
   return(list(imp, error))

   ############################################################################################################
   ############################################################################################################
   ############################################################################################################
   
   
  
  # 
  # file_path <- sprintf("TrainingPredictionBias_%s.RData",missing_elements)
  # 
  # if (file.exists(file_path)) {
  #   load(file_path)
  #   existing_data <- rbind(existing_data, bias)
  # } else {
  #   existing_data <- bias
  # }
  # save(existing_data, file = file_path)
  # 
  

  

  
}



