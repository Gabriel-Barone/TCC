series_generation_function <- function (A, Tf, predictors, rho) {
  
  Tf = Tf+1
  
  # Parameter to generate our AR processes
  sd.rho <- sqrt(1-rho^2)
  
  #Creating empty dataframe
  data <- as.data.frame(matrix(NA, A*Tf, 2*length(predictors) + 3))
  colnames(data) <- c("unit", "t", paste0("H", 1:length(predictors)), paste0("H*", 1:length(predictors)), "V")
  
  # Create the unit identifiers.
  for (a in 1:A) {
    data[(1 + (a-1)*Tf):(Tf + (a-1)*Tf), 1] <- rep(a, Tf)
  }
  
  # Create the time identifiers.
  data[, 2] <- rep(seq(from = 0, to  = (Tf-1), by = 1), A)
  
  #CREATING TRUE PREDICTORS-----------------
  # Create the common time trend
  delta <- matrix(rnorm(Tf * length(predictors), mean = 0, sd = 1), Tf, length(predictors))
  
  # Create the stationary components
  lambda <- matrix(NA, Tf, 0)
  for (g in 1:length(predictors)) {
    for (i in 1:predictors[g]) {
      col_name <- paste0("predictor_", g, "_group_", i)
      new_columns <- matrix(NA, Tf, 1)
      new_columns <- matrix(arima.sim(list(ar = rho), n =Tf, sd=sd.rho))
      colnames(new_columns) <- col_name
      lambda <- cbind(lambda, new_columns)
    }
  }
  
  #Create groups
  groups <- matrix(NA, A, 0)
  for (g in 1:length((predictors))) {

    col_name <- paste0("predictor_", g)
    new_columns <- matrix(sample(1:A,A), A, 1)
    colnames(new_columns) <- col_name
    groups <- cbind(groups, new_columns)
  }
  
  
  # Create variables
  for (a in 1:A) {
  
    count_lambda <- 1
    
    for (g in 1:length(predictors)) {

      lambda_coef <-  matrix(0, predictors[g], 1)
      lambda_coef[ceiling(groups[a,g]/(A/predictors[g])),1] <- 1
      
      # Error term
      epsilon <- matrix(rnorm(Tf, mean = 0, sd = 0.1), Tf, 1)
      
      # Serie
      data[(1 + (a-1)*Tf):(Tf + (a-1)*Tf), g + 2] <-
        delta[,g] + lambda[,count_lambda :(count_lambda + predictors[g] - 1)] %*% lambda_coef +
         + epsilon
      count_lambda <- count_lambda + predictors[g]
    }
  }
  
  #CREATING FALSE PREDICTORS---------------
  # Create the common time trend
  delta <- matrix(rnorm(Tf * length(predictors), mean = 0, sd = 1), Tf, length(predictors))
  
  # Create the stationary components
  lambda <- matrix(NA, Tf, 0)
  for (g in 1:length(predictors)) {
    for (i in 1:predictors[g]) {
      col_name <- paste0("predictor_", g, "_group_", i)
      new_columns <- matrix(NA, Tf, 1)
      new_columns <- matrix(arima.sim(list(ar = rho), n =Tf, sd=sd.rho))
      colnames(new_columns) <- col_name
      lambda <- cbind(lambda, new_columns)
    }
  }
  
  #Create groups
  groups <- matrix(NA, A, 0)
  for (g in 1:length((predictors))) {
    
    col_name <- paste0("predictor_", g)
    new_columns <- matrix(sample(1:A,A), A, 1)
    colnames(new_columns) <- col_name
    groups <- cbind(groups, new_columns)
  }
  
  # Create variables
  for (a in 1:A) {
    
    count_lambda <- 1
    
    for (g in 1:length(predictors)) {
      
      lambda_coef <-  matrix(0, predictors[g], 1)
      lambda_coef[ceiling(groups[a,g]/(A/predictors[g])),1] <- 1
      
      # Error term
      epsilon <- matrix(rnorm(Tf, mean = 0, sd = 0.1), Tf, 1)
      
      # Serie
      data[(1 + (a-1)*Tf):(Tf + (a-1)*Tf), length(predictors)+ g + 2] <-
        delta[,g] + lambda[,count_lambda :(count_lambda + predictors[g] - 1)] %*% lambda_coef +
        + epsilon
      count_lambda <- count_lambda + predictors[g]
    }
  }
  
  #CREATE V-----------
  
  phi_v <- runif(1, -1, 1)
  phi <- rnorm(length(predictors), 0, 1)
  epsilon <- matrix(rnorm(A*Tf, mean = 0, sd = 0.1), A*Tf, 1)
  for (a in 1:A) {
    data[(2 + (a-1)*Tf), 2*length(predictors) + 3] <- sum(phi * data[(1 + (a-1)*Tf), 3:(length(predictors) + 2)]) + epsilon[(2 + (a-1)*Tf),]
    for (i in (3 + (a-1)*Tf):(Tf + (a-1)*Tf)) {
      data[i, 2*length(predictors) + 3] <- sum(phi * data[i-1, 3:(length(predictors) + 2)]) + phi_v * data[i-1, 2*length(predictors)+3] + epsilon[i,]
      }
  }
  
  #OMIT NA --------------
  
  data <- na.omit(data)
  row.names(data) <- 1:nrow(data)
  
  return(data)
}

#times_observed
Tf = 60

#list_of_predictors_and_how_many_groups_of_variables_share_them
predictors<- c(2,4,5,10)

#number_of_units
A = 20

#AR(1)_param
rho = 0.5

#CREATE FIRST 1000 SERIES ------------

# Create an empty data frame to store the results
result_df <- data.frame()

# Set the seed for reproducibility
set.seed(1)
# Generate 10 series
for (i in 1:1000) {
  # Generate series
  data <- series_generation_function(A, Tf, predictors, rho)
  
  # Add series_id column
  data$series_id <- i
  
  # Append the data to the result_df
  result_df <- rbind(result_df, data)
  
  cat("Processing iteration", i, "out of 1000\n")
}

# Write the result to a CSV file
setwd("C:/Users/gpime/OneDrive/Área de Trabalho/TCC")
write.csv(result_df, "tcc_output_data_1.csv", row.names = FALSE)




#CREATE SECOND 1000 SERIES ------------

# Create an empty data frame to store the results
result_df <- data.frame()

# Set the seed for reproducibility
set.seed(2)
# Generate 10 series
for (i in 1:1000) {
  # Generate series
  data <- series_generation_function(A, Tf, predictors, rho)
  
  # Add series_id column
  data$series_id <- i
  
  # Append the data to the result_df
  result_df <- rbind(result_df, data)
  
  cat("Processing iteration", i, "out of 1000\n")
}

# Write the result to a CSV file
setwd("C:/Users/gpime/OneDrive/Área de Trabalho/TCC")
write.csv(result_df, "tcc_output_data_2.csv", row.names = FALSE)




#CREATE THIRD 1000 SERIES ------------

# Create an empty data frame to store the results
result_df <- data.frame()

# Set the seed for reproducibility
set.seed(3)
# Generate 10 series
for (i in 1:1000) {
  # Generate series
  data <- series_generation_function(A, Tf, predictors, rho)
  
  # Add series_id column
  data$series_id <- i
  
  # Append the data to the result_df
  result_df <- rbind(result_df, data)
  
  cat("Processing iteration", i, "out of 1000\n")
}

# Write the result to a CSV file
setwd("C:/Users/gpime/OneDrive/Área de Trabalho/TCC")
write.csv(result_df, "tcc_output_data_3.csv", row.names = FALSE)




#CREATE FOURTH 1000 SERIES ------------

# Create an empty data frame to store the results
result_df <- data.frame()

# Set the seed for reproducibility
set.seed(4)
# Generate 10 series
for (i in 1:1000) {
  # Generate series
  data <- series_generation_function(A, Tf, predictors, rho)
  
  # Add series_id column
  data$series_id <- i
  
  # Append the data to the result_df
  result_df <- rbind(result_df, data)
  
  cat("Processing iteration", i, "out of 1000\n")
}

# Write the result to a CSV file
setwd("C:/Users/gpime/OneDrive/Área de Trabalho/TCC")
write.csv(result_df, "tcc_output_data_4.csv", row.names = FALSE)




#CREATE FITH 1000 SERIES ------------

# Create an empty data frame to store the results
result_df <- data.frame()

# Set the seed for reproducibility
set.seed(5)
# Generate 10 series
for (i in 1:1000) {
  # Generate series
  data <- series_generation_function(A, Tf, predictors, rho)
  
  # Add series_id column
  data$series_id <- i
  
  # Append the data to the result_df
  result_df <- rbind(result_df, data)
  
  cat("Processing iteration", i, "out of 1000\n")
}

# Write the result to a CSV file
setwd("C:/Users/gpime/OneDrive/Área de Trabalho/TCC")
write.csv(result_df, "tcc_output_data_5.csv", row.names = FALSE)





