mice.impute.xgbParam<- function(y, ry, x, wy = NULL,xgb.params_call = list(),dm,...) {
  
  if (is.null(wy)) wy <- !ry
  
  
  nmis <- sum(wy)
  xobs <- x[ry, , drop = FALSE]
  xmis <- x[wy, , drop = FALSE]
  yobs <- y[ry]
  
  
  xgb.params<-xgb.params_call[[1]]
  nrounds<-xgb.params_call[[2]]
  

  dobs <- xgboost::xgb.DMatrix(data = xobs, label = yobs)
  dmis <- xgboost::xgb.DMatrix(data = xmis)
  
  watchlist <- list(train = dobs)
  

  xgb.fit <- xgboost::xgb.train(
    data = dobs, 
    params = xgb.params, nrounds = nrounds, early_stopping_rounds = NULL, watchlist=watchlist, verbose = 0, ...
  )
  
  #yhatmis <- predict(xgb.fit, dmis)
  #yhatmis
  
  
  
  allCols<-c("x", "y", "z")
  missing_elements <- allCols[!(allCols %in% colnames(x))]
  yhatobs <- predict(xgb.fit, xobs)
  bias <- mean(yobs-yhatobs)
  names(bias) <- missing_elements
  
  file_path <- sprintf("XGBParam_TrainingPredictionBias_%s.RData",missing_elements)
  
  if (file.exists(file_path)) {
    load(file_path)
    existing_data <- rbind(existing_data, bias)
  } else {
    existing_data <- bias
  }
  save(existing_data, file = file_path)
  
browser()
  yhatmis <- predict(xgb.fit, rbind(xobs,xmis), predleaf = TRUE)
  nodes_obs <- yhatmis[1:nrow(xobs), , drop = FALSE]
  nodes_mis <- yhatmis[(nrow(xobs) + 1):nrow(yhatmis), , drop = FALSE]
  
  select_donors <- function(i) {
    # Function to extract all eligible donors for each missing value
    print(i)
    donors <- split(yobs, nodes_obs[, i])
    donors[as.character(nodes_mis[, i])]
  }
  
  donors_selected<-sapply(seq_len(nrounds), FUN = select_donors)
  
  matched_sample<-apply(donors_selected, MARGIN = 1, FUN = function(s) sample(unlist(s), 1))
  
  matched_sample
}

