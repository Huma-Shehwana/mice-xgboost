



xgb_param_calc<-function(data,response = NULL, select_features=NULL, num_cores = 6){
  
  
  library(xgboost)
  library(insuranceData) # example dataset https://cran.r-project.org/web/packages/insuranceData/insuranceData.pdf
  library(tidyverse)
  library(mlrMBO)
  library(rBayesianOptimization) 
  
  
  ##############################################################################
  
  par_opt <- function(cc.data, select_features, response, nfold = 5, seedval = 12345, early_stopping_rounds = 20, nround = 200){
    
    p <- length(select_features) + 1
    if (p == 2) {
      obs.data <- Matrix::sparse.model.matrix(reformulate(select_features, response), data = cc.data)
    } else {
      obs.data <- Matrix::sparse.model.matrix(reformulate(select_features, response), data = cc.data)[, -1]
    }
    response_data <- as.vector(unlist(cc.data[,response]))
    
    
    dtrain <- xgboost::xgb.DMatrix(obs.data, label = response_data, missing = NA)
    
    browser()
    cv_folds = rBayesianOptimization::KFold(response_data, 
                                            nfolds= nfold,
                                            seed= seedval)
    
    obj.fun  <- smoof::makeSingleObjectiveFunction(
      name = "xgb_cv_bayes",
      fn =   function(x){
        set.seed(seedval)
        cv <- xgb.cv(params = list(
          booster = "gbtree",
          eta = x["eta"],
          max_depth = x["max_depth"],
          min_child_weight = x["min_child_weight"],
          gamma = x["gamma"],
          alpha = x["alpha"],
          lambda = x["lambda"],
          objective = 'reg:squarederror', 
          eval_metric     = "mape"),
          data = dtrain,
          nround = nround,
          early_stopping_rounds = early_stopping_rounds,
          folds = cv_folds,
          prediction = TRUE,
          maximize = FALSE,
          showsd = TRUE,
          verbose = 0)
        
        min(cv$evaluation_log[, 'test_mape_mean'])
      },
      par.set = makeParamSet(
        makeNumericParam("eta",lower = 0.001, upper = 1),
        makeNumericParam("alpha",lower = 0, upper = 100),
        makeNumericParam("lambda",lower = 0, upper = 100),
        makeNumericParam("gamma",lower = 0,upper = 10),
        makeIntegerParam("max_depth",lower= 2,upper = 20),
        makeIntegerParam("min_child_weight", lower= 0,upper = 50)
      ),
      minimize = TRUE
    )
    
    
    set.seed(seedval)
    control = makeMBOControl()
    des = generateDesign(n=10,
                         par.set = getParamSet(obj.fun), 
                         fun = lhs::randomLHS) 
    
    control = setMBOControlTermination(control, iters = 10)
    run = mbo(fun = obj.fun, 
              control = control, 
              design = des)
    
    
    best_solution <- run$opt.path$env$path[which.min(run$opt.path$env$path$y),]
    best_MAPE <- best_solution$y
    best_parameters <- run$x
    
    best_parameters
  }
  
  
  ##############################################################################
  
  
  num.cc <- sum(complete.cases(data))
  
  if (num.cc == 0) {
    stop("No complete cases in this dataset.")
  }
  
  
  if (num.cc > 0 & num.cc < 50) {
    warnings("Less than 10 complete cases. Results may not be reliable.")
  }
  
  
  cc.data <- data[complete.cases(data), ]
  Names <- colnames(cc.data)
  Types <- as.vector(apply(data,2,class))
  
  if (any(Types == "character")) {
    stop("This datset contains variables with character type. Please set stringsAsFactors = TRUE")
  }
  
  na.col <- which(colSums(is.na(data))!=0)
  
  
  if (is.numeric(select_features)) {
    select_features <- Names[select_features]
  } else if (is.null(select_features) & response !="all") {
    select_features <- setdiff(Names, response)
  } else if (response =="all" & is.null(select_features)){
    select_features = Names
  }
  
  if (is.null(response)) {
    # r.idx <- sample(1:ncol(cc.data), size = 1)
    r.idx <- sample(na.col, size = 1)
    response <- Names[r.idx]
  } else if (is.numeric(response)) {
    response <- Names[response]
  } else if(response =="all"){
    response <- Names
  }
  


  
  
  params = list()
  
  
  if(length(response)>1){
    for(i in response){
      response_1 = i
      select_features <- setdiff(Names, response_1)
      browser()
      params[[i]] <- par_opt(cc.data, select_features, response_1)
      
    }
  } else {
    params = par_opt(cc.data, select_features, response)
  }
  
 
  return(params)
  
}






