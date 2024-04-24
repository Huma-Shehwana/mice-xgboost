

xgb_param_calc<-function(data,response = NULL, select_features=NULL, num_cores = -1, iter=50){
  

  
  par_opt <- function(cc.data, training_features, response_var, nfold = 5, seedval = seed, early_stopping_rounds = 20, nround = 200){
    
    p <- length(training_features) + 1
    
    if (p == 2) {
      obs.data <- Matrix::sparse.model.matrix(reformulate(training_features, response_var), data = cc.data)
    } else {
      obs.data <- Matrix::sparse.model.matrix(reformulate(training_features, response_var), data = cc.data)[, -1]
    }
    response_data <- as.vector(unlist(cc.data[,response_var]))
    
    
    dtrain <- xgboost::xgb.DMatrix(obs.data, label = response_data, missing = NA)
    
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
                                                        eval_metric = "mape"),
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
    
    control = setMBOControlTermination(control, iters = iter)
    run = mbo(fun = obj.fun, 
              control = control, 
              design = des)
    

    plot_fig<-run$opt.path$env$path  %>% 
              mutate(Round = row_number()) %>% ggplot(aes(x= Round, y= y)) + 
              geom_point() +
              labs(title = sprintf("Response Variable: %s , Features: %s", response_var, paste(training_features,collapse = ",")))+
              ylab("MAPE") + theme(plot.title = element_text(hjust = 0.5))
    

    best_solution <- run$opt.path$env$path[which.min(run$opt.path$env$path$y),]
    best_MAPE <- best_solution$y
    best_parameters <- run$x
    
    list("parameter" = best_parameters, "fig" = plot_fig)
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
  
  #cc.data <- data
  
  Names <- colnames(cc.data)
  Types <- as.vector(apply(cc.data,2,class))
  
  if (any(Types == "character")) {
    stop("This datset contains variables with character type. Please set stringsAsFactors = TRUE")
  }
  
  na.col <- which(colSums(is.na(data))!=0)
  
  
  
  
  if(!is.null(select_features) &!identical(select_features, "all") & !all(select_features %in% seq_along(Names))){
    stop("Please check the format of select_features. It can be either NULL, all or index of columns that you want to include as a feature")
  }
  
  if (is.null(response)) { 
              if(is.null(select_features) | identical(select_features, "all")) {
                        r.idx <- sample(na.col, size = 1)
                        response_var <- Names[r.idx]
                        training_features <- setdiff(Names, response_var)
              }else if(is.numeric(select_features)){
                        training_features <- Names[select_features]
                        response_var <- sample(setdiff(Names[na.col], training_features), size = 1)
              } 

    } else if (is.numeric(response)) {
                        if(!all(response %in% seq_along(Names))){
                                  stop("Please check the format of response variable It can be either NULL, all or index of columns that you want to use as a target variable")
                          }
                  if(length(response)>1){ 
                        message("You have selected more than one variables for parameter optimization. Please note that only the first index will be used for parameter estimation")
                  } 
                  response_var <- Names[response[1]]
                  if(is.null(select_features) | identical(select_features, "all")) {
                          training_features <- setdiff(Names, response_var)
                  } else if(is.numeric(select_features)){
                    if(response %in% select_features){
                      message("Response and feature variables should not be same. All variables other than response variable will be used as features")
                      training_features <- setdiff(Names, response_var)
                      } else {
                        training_features <- Names[select_features]
                      }
                  } 
    } else if(response =="all"){
    response_var <- Names
    training_features <- NULL
    message("You have chosen all variables to be iteratively used as a target variable for parameter optimization. Please note that select_features will be automatically chosen for this process.")
    } else {
      stop("Please check the format of response variable It can be either NULL, \"all\" or index of column that you want to use for parameter optimization")
      }
  
  params = list()
  
  if(length(response_var)>1){
    message("Performing bayesian optimization iteratively using each vairable as response variable")
    params <- future_map(setNames(response_var, response_var), ~par_opt(cc.data, setdiff(Names, .x), .x))
    
  } else {
    message("Performing bayesian optimization using ", response_var, " as a response variable and ", training_features, " as features")
    params = par_opt(cc.data, training_features, response_var)
  }
  
 
  return(params)
  
}






