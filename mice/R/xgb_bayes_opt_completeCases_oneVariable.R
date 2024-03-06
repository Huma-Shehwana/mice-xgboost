
#' Title
#'
#' @param data 
#' @param response 
#' @param select_features 
#'
#' @return
#' @export
#'
#' @examples
#' 
xgb_bayes_opt_completeCases_oneVariable<-function(data,response = NULL, select_features=NULL, num_cores = 6){
  

  num.cc <- sum(complete.cases(data))
  
  if (num.cc == 0) {
    stop("No complete cases in this dataset.")
  }
  
  
  if (num.cc > 0 & num.cc < 10) {
    warnings("Less than 10 complete cases. Results may not be reliable.")
  }
  
 
  cc.data <- data[complete.cases(data), ]
  Names <- colnames(cc.data)
  Types <- as.vector(apply(data,2,class))
  
  if (any(Types == "character")) {
    stop("This datset contains variables with character type. Please set stringsAsFactors = TRUE")
  }
  
  na.col <- which(colSums(is.na(data)) != 0)
  
  if (is.null(response)) {
    # r.idx <- sample(1:ncol(cc.data), size = 1)
    r.idx <- sample(na.col, size = 1)
    response <- Names[r.idx]
    response
  } else if (!is.character(response)) {
    response <- Names[response]
  }
  
  
  if (is.null(select_features)) {
    select_features <- setdiff(Names, response)
  } else if (!is.character(select_features)) {
    select_features <- Names[select_features]
  }
  
  p <- length(select_features) + 1
  if (p == 2) {
    obs.data <- Matrix::sparse.model.matrix(reformulate(select_features, response), data = cc.data)
  } else {
    obs.data <- Matrix::sparse.model.matrix(reformulate(select_features, response), data = cc.data)[, -1]
  }
  response_data <- as.vector(unlist(cc.data[,response]))
  
  
  
  nfold = 5
  
  cv_folds = rBayesianOptimization::KFold(response_data, 
                                          nfolds= nfold,
                                          seed= 0)
  
  
scoring_function <- function(eta, gamma, max_depth, min_child_weight, subsample,lambda = lambda,alpha = alpha, nfold = 5) {
    
    dtrain <- xgboost::xgb.DMatrix(obs.data, label = response_data, missing = NA)
    pars <- list(
      eta = eta, gamma = gamma, max_depth = max_depth, min_child_weight = min_child_weight,
      lambda=lambda, alpha = alpha, booster = "gbtree",
      objective = "reg:squarederror",
      eval_metric = "mape", #Kullback–Leibler divergence # change evaluation metric
      verbosity = 0, nfold = nfold
    )
    set.seed(12345)
    
    xgbcv <- xgboost::xgb.cv(
      params = pars,
      data = dtrain,
      folds = cv_folds,
      nrounds = 200,
      prediction = TRUE,
      showsd = TRUE,
      early_stopping_rounds = 20,
      maximize = FALSE
      #stratified = TRUE
    )
    

    return(
      list(
        Score = -min(xgbcv$evaluation_log$test_mape_mean),
        nrounds = xgbcv$best_iteration
      )
    )
  }
  

bounds <- list(eta = c(0.001, 1),
               gamma=c(0,10),
               max_depth = c(2L, 20L),
               min_child_weight = c(0, 50),
               lambda = c(0, 100),
               alpha = c(0, 100)
)
  
  
 # cl <- makeCluster(num_cores)
#  registerDoParallel(cl)
 # clusterExport(cl, c("obs.data", "response_data", "nfold", "bounds", "scoring_function"), envir = environment())
  
  # Load required libraries on each cluster
  # clusterEvalQ(cl, {
  #   library(ParBayesianOptimization)
  #   library(Matrix)
  #   library(xgboost)
  # })

     bayes_out <- ParBayesianOptimization::bayesOpt(FUN = scoring_function, 
                           bounds = bounds, 
                           initPoints = 10, 
                           iters.n = 10,
                           cv_folds = cv_folds,  # Pass cv_folds to the function
                           parallel = FALSE, verbose = FALSE)
  
    stopCluster(cl)
     
  
  opt_params <- append(list(booster = "gbtree", 
                            objective = "reg:squarederror", 
                            eval_metric = "mape"), 
                       ParBayesianOptimization::getBestPars(bayes_out))
  browser()
  
  folds <- list(fold1 = as.integer(seq(1, nrow(obs.data), by = 5)),
                fold2 = as.integer(seq(2, nrow(obs.data), by = 5)),
                fold3 = as.integer(seq(3, nrow(obs.data), by = 5)),
                fold4 = as.integer(seq(4, nrow(obs.data), by = 5)),
                fold5 = as.integer(seq(5, nrow(obs.data), by = 5)))
  
  # Run cross validation 
  xgbcv <- xgboost::xgb.cv(params = opt_params,
                  data = obs.data,
                  label = response_data,
                  nround = 200,
                  folds = cv_folds,
                  prediction = TRUE,
                  early_stopping_rounds = 20,
                  verbose = 0,
                  maximize = F)
  
  nrounds = xgbcv$best_iteration
  #opt_params$nrounds <- nrounds
  
  return(list(opt_params, nrounds))
  
}

















