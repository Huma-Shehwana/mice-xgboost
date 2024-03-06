data <- read.table("Simulations_1000.tsv", header = TRUE)

estimates <- c("bias", "cov", "width", "prop")

columns_bias <- data[, grep("bias", names(data))]
columns_cov <- data[, grep("cov", names(data))]
columns_width <- data[, grep("width", names(data))]
columns_prop <- data[, grep("prop", names(data))]

result <- cbind(columns_bias, columns_cov, columns_width)

write.table(result, "Formatted_sim_1000.txt",sep = "\t")
