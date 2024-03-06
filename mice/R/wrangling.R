methods<- c("norm","rf","cart","xgb","mixgb","xgbParam")
pattern = "sim_5_20_%s.RData"

norm_res <- load(sprintf(pattern, methods[[1]]))
norm_resEval <- results$eval_result[,c("true","cov","bias","width")]
colnames(norm_resEval) = paste(colnames(norm_resEval), "norm", sep = "_")


rf_res <- load(sprintf(pattern, methods[[2]]))
rf_resEval <- results$eval_result[,c("cov","bias","width")]
colnames(rf_resEval) = paste(colnames(rf_resEval), "RF", sep = "_")


cart_res <- load(sprintf(pattern, methods[[3]]))
cart_resEval <- results$eval_result[,c("cov","bias","width")]
colnames(cart_resEval) = paste(colnames(cart_resEval), "CART", sep = "_")


xgb_res <- load(sprintf(pattern, methods[[4]]))
xgb_resEval <- results$eval_result[,c("cov","bias","width")]
colnames(xgb_resEval) = paste(colnames(xgb_resEval), "XGB", sep = "_")


mixgb_res <- load(sprintf(pattern, methods[[5]]))
mixgb_resEval <- results$eval_result[,c("cov","bias","width")]
colnames(mixgb_resEval) = paste(colnames(mixgb_resEval), "mixgb", sep = "_")

xgbParam_res <- load(sprintf(pattern, methods[[6]]))
xgbParam_resEval <- results$eval_result[,c("cov","bias","width")]
colnames(xgbParam_resEval) = paste(colnames(xgbParam_resEval), "xgbParam", sep = "_")


allres <- cbind(norm_resEval,rf_resEval, cart_resEval, xgb_resEval,mixgb_resEval,xgbParam_resEval )

allres<-allres[,c(1,
            grep("cov_*", colnames(allres)), grep("bias_*", colnames(allres)), grep("width_*", colnames(allres)))]

save(allres,file ="SimResults_20%.RData")

