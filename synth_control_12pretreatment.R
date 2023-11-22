require(tidysynth)
require(dplyr)

setwd("C:/Users/gpime/OneDrive/Área de Trabalho/TCC")

# Read the CSV file into a data frame
df <- read.csv("tcc_output_data_1.csv")

#Creating empty dataframe
diagnostic <- as.data.frame(matrix(NA, 4000, 5))
colnames(diagnostic) <- c("series_id", "synth_control_type", "p_value_<=0.05", "p_value_<=0.1", "R^2")

# Create the unit identifiers.
for (r in 1:(nrow(diagnostic)/4)) {
  diagnostic[(4*r-3):(4*r), 1] <- rep(r, 4)
}

diagnostic <- diagnostic %>%
  group_by(series_id) %>%
  mutate(synth_control_type = c('T_predictors', 'F_predictors', 'ALL_predictors', 'NO_predictors')[row_number()]) %>%
  ungroup()

diagnostic['p_value_<=0.05'] <- 0
diagnostic['p_value_<=0.1'] <- 0

calculate_R2 <- function(data) {
  
  # Extract real_y and synth_y from the filtered data
  real_y_filtered <- data$real_y
  synth_y_filtered <- data$synth_y
  # Calculate mean of real_y
  mean_real_y <- mean(real_y_filtered)
  # Calculate R^2
  numerator <- sum((real_y_filtered - synth_y_filtered)^2)
  denominator <- sum((real_y_filtered - mean_real_y)^2)
  r_squared <- 1 - (numerator / denominator)
  
}

#FIRST SERIES-----------

iteration_count <- 0

for(r in 1:(nrow(diagnostic)/4)) {
  print(r)
  
  series_used <- df[df$series_id == r,]
  
  #TRUE PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                        # V = mean(V, na.rm = T),
                         H1 = mean(`H1`, na.rm = T),
                         H2 = mean(`H2`, na.rm = T),
                         H3 = mean(`H3`, na.rm = T),
                         H4 = mean(`H4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
      V_38 = V) %>%
      generate_predictor(time_window = 43,
      V_43 = V) %>%
      generate_predictor(time_window = 49,
      V_49 = V) %>%
      
    
    # Generate the fitted weights for the synthetic control
    generate_weights(optimization_window = 38:50, # time to use in the optimization task
                     ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-3)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-3)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-3)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-3)] <- calculate_R2(synth_control)
    
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-3)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-3)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-3)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
    
  #FALSE PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H.1`, na.rm = T),
                         H2 = mean(`H.2`, na.rm = T),
                         H3 = mean(`H.3`, na.rm = T),
                         H4 = mean(`H.4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-2)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-2)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-2)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-2)] <- calculate_R2(synth_control)
    
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-2)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-2)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-2)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
      
  #ALL PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H1`, na.rm = T),
                         H2 = mean(`H2`, na.rm = T),
                         H3 = mean(`H3`, na.rm = T),
                         H4 = mean(`H4`, na.rm = T),
                         `H*1` = mean(`H.1`, na.rm = T),
                         `H*2` = mean(`H.2`, na.rm = T),
                         `H*3` = mean(`H.3`, na.rm = T),
                         `H*4` = mean(`H.4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-1)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-1)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-1)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-1)] <- calculate_R2(synth_control)
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-1)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-1)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-1)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
      
  #NO PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r)] = 1
      diagnostic$`p_value_<=0.1`[(4*r)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r)] <- calculate_R2(synth_control)
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4 * r)] <- NA
    diagnostic$`p_value_<=0.05`[(4 * r)] <- NA
    diagnostic$`p_value_<=0.1`[(4 * r)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  # Check if the current iteration is a multiple of 50
  if (r %% 50 == 0) {
    # Increment the iteration count
    iteration_count <- iteration_count + 1
    
    # Construct the filename based on the current range of r
    filename <- paste0("tcc_diagnosticc_1_to_", (iteration_count * 50), ".csv")
    
    # Save the diagnostic data frame to CSV
    write.csv(diagnostic, filename, row.names = FALSE)
  }
}

#SECOND SERIES-----------

# Read the CSV file into a data frame
df <- read.csv("tcc_output_data_2.csv")

iteration_count <- 0

for(r in 1:(nrow(diagnostic)/4)) {
  print(r)
  
  series_used <- df[df$series_id == r,]
  
  #TRUE PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H1`, na.rm = T),
                         H2 = mean(`H2`, na.rm = T),
                         H3 = mean(`H3`, na.rm = T),
                         H4 = mean(`H4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-3)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-3)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-3)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-3)] <- calculate_R2(synth_control)
    
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-3)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-3)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-3)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #FALSE PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H.1`, na.rm = T),
                         H2 = mean(`H.2`, na.rm = T),
                         H3 = mean(`H.3`, na.rm = T),
                         H4 = mean(`H.4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-2)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-2)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-2)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-2)] <- calculate_R2(synth_control)
    
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-2)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-2)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-2)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #ALL PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H1`, na.rm = T),
                         H2 = mean(`H2`, na.rm = T),
                         H3 = mean(`H3`, na.rm = T),
                         H4 = mean(`H4`, na.rm = T),
                         `H*1` = mean(`H.1`, na.rm = T),
                         `H*2` = mean(`H.2`, na.rm = T),
                         `H*3` = mean(`H.3`, na.rm = T),
                         `H*4` = mean(`H.4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-1)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-1)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-1)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-1)] <- calculate_R2(synth_control)
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-1)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-1)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-1)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #NO PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r)] = 1
      diagnostic$`p_value_<=0.1`[(4*r)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r)] <- calculate_R2(synth_control)
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4 * r)] <- NA
    diagnostic$`p_value_<=0.05`[(4 * r)] <- NA
    diagnostic$`p_value_<=0.1`[(4 * r)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  # Check if the current iteration is a multiple of 50
  if (r %% 50 == 0) {
    # Increment the iteration count
    iteration_count <- iteration_count + 1
    
    # Construct the filename based on the current range of r
    filename <- paste0("tcc_diagnosticc_1001_to_", (iteration_count * 50 + 1001), ".csv")
    
    # Save the diagnostic data frame to CSV
    write.csv(diagnostic, filename, row.names = FALSE)
  }
}


#THIRD SERIES-----------

# Read the CSV file into a data frame
df <- read.csv("tcc_output_data_3.csv")

iteration_count <- 0

for(r in 1:(nrow(diagnostic)/4)) {
  print(r)
  
  series_used <- df[df$series_id == r,]
  
  #TRUE PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H1`, na.rm = T),
                         H2 = mean(`H2`, na.rm = T),
                         H3 = mean(`H3`, na.rm = T),
                         H4 = mean(`H4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-3)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-3)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-3)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-3)] <- calculate_R2(synth_control)
    
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-3)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-3)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-3)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #FALSE PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H.1`, na.rm = T),
                         H2 = mean(`H.2`, na.rm = T),
                         H3 = mean(`H.3`, na.rm = T),
                         H4 = mean(`H.4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-2)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-2)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-2)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-2)] <- calculate_R2(synth_control)
    
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-2)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-2)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-2)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #ALL PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H1`, na.rm = T),
                         H2 = mean(`H2`, na.rm = T),
                         H3 = mean(`H3`, na.rm = T),
                         H4 = mean(`H4`, na.rm = T),
                         `H*1` = mean(`H.1`, na.rm = T),
                         `H*2` = mean(`H.2`, na.rm = T),
                         `H*3` = mean(`H.3`, na.rm = T),
                         `H*4` = mean(`H.4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-1)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-1)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-1)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-1)] <- calculate_R2(synth_control)
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-1)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-1)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-1)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #NO PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r)] = 1
      diagnostic$`p_value_<=0.1`[(4*r)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r)] <- calculate_R2(synth_control)
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4 * r)] <- NA
    diagnostic$`p_value_<=0.05`[(4 * r)] <- NA
    diagnostic$`p_value_<=0.1`[(4 * r)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  # Check if the current iteration is a multiple of 50
  if (r %% 50 == 0) {
    # Increment the iteration count
    iteration_count <- iteration_count + 1
    
    # Construct the filename based on the current range of r
    filename <- paste0("tcc_diagnosticc_2001_to_", (iteration_count * 50 + 2001), ".csv")
    
    # Save the diagnostic data frame to CSV
    write.csv(diagnostic, filename, row.names = FALSE)
  }
}

#FOURTH SERIES-----------

# Read the CSV file into a data frame
df <- read.csv("tcc_output_data_4.csv")

iteration_count <- 0

for(r in 1:(nrow(diagnostic)/4)) {
  print(r)
  
  series_used <- df[df$series_id == r,]
  
  #TRUE PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H1`, na.rm = T),
                         H2 = mean(`H2`, na.rm = T),
                         H3 = mean(`H3`, na.rm = T),
                         H4 = mean(`H4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-3)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-3)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-3)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-3)] <- calculate_R2(synth_control)
    
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-3)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-3)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-3)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #FALSE PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H.1`, na.rm = T),
                         H2 = mean(`H.2`, na.rm = T),
                         H3 = mean(`H.3`, na.rm = T),
                         H4 = mean(`H.4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-2)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-2)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-2)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-2)] <- calculate_R2(synth_control)
    
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-2)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-2)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-2)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #ALL PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H1`, na.rm = T),
                         H2 = mean(`H2`, na.rm = T),
                         H3 = mean(`H3`, na.rm = T),
                         H4 = mean(`H4`, na.rm = T),
                         `H*1` = mean(`H.1`, na.rm = T),
                         `H*2` = mean(`H.2`, na.rm = T),
                         `H*3` = mean(`H.3`, na.rm = T),
                         `H*4` = mean(`H.4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-1)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-1)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-1)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-1)] <- calculate_R2(synth_control)
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-1)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-1)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-1)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #NO PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r)] = 1
      diagnostic$`p_value_<=0.1`[(4*r)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r)] <- calculate_R2(synth_control)
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4 * r)] <- NA
    diagnostic$`p_value_<=0.05`[(4 * r)] <- NA
    diagnostic$`p_value_<=0.1`[(4 * r)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  # Check if the current iteration is a multiple of 50
  if (r %% 50 == 0) {
    # Increment the iteration count
    iteration_count <- iteration_count + 1
    
    # Construct the filename based on the current range of r
    filename <- paste0("tcc_diagnosticc_3001_to_", (iteration_count * 50 + 3001), ".csv")
    
    # Save the diagnostic data frame to CSV
    write.csv(diagnostic, filename, row.names = FALSE)
  }
}

#FITH SERIES-----------

# Read the CSV file into a data frame
df <- read.csv("tcc_output_data_5.csv")

iteration_count <- 0

for(r in 1:(nrow(diagnostic)/4)) {
  print(r)
  
  series_used <- df[df$series_id == r,]
  
  #TRUE PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H1`, na.rm = T),
                         H2 = mean(`H2`, na.rm = T),
                         H3 = mean(`H3`, na.rm = T),
                         H4 = mean(`H4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-3)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-3)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-3)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-3)] <- calculate_R2(synth_control)
    
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-3)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-3)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-3)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #FALSE PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H.1`, na.rm = T),
                         H2 = mean(`H.2`, na.rm = T),
                         H3 = mean(`H.3`, na.rm = T),
                         H4 = mean(`H.4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-2)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-2)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-2)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-2)] <- calculate_R2(synth_control)
    
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-2)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-2)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-2)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #ALL PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Generate the aggregate predictors used to fit the weights
      generate_predictor(time_window = 38:50,
                         # V = mean(V, na.rm = T),
                         H1 = mean(`H1`, na.rm = T),
                         H2 = mean(`H2`, na.rm = T),
                         H3 = mean(`H3`, na.rm = T),
                         H4 = mean(`H4`, na.rm = T),
                         `H*1` = mean(`H.1`, na.rm = T),
                         `H*2` = mean(`H.2`, na.rm = T),
                         `H*3` = mean(`H.3`, na.rm = T),
                         `H*4` = mean(`H.4`, na.rm = T)) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r-1)] = 1
      diagnostic$`p_value_<=0.1`[(4*r-1)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r-1)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r-1)] <- calculate_R2(synth_control)
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4*r-1)] <- NA
    diagnostic$`p_value_<=0.05`[(4*r-1)] <- NA
    diagnostic$`p_value_<=0.1`[(4*r-1)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  #NO PREDICTORS
  tryCatch({
    V_out_predictors <-
      
      series_used %>%
      
      # initial the synthetic control object
      synthetic_control(outcome = V, # outcome
                        unit = unit, # unit index in the panel data
                        time = t, # time index in the panel data
                        i_unit = 1, # unit where the intervention occurred
                        i_time = 50, # time period when the intervention occurred
                        generate_placebos=T # generate placebo synthetic controls (for inference)
      ) %>%
      
      # Lagged V 
      generate_predictor(time_window = 38,
                         V_38 = V) %>%
      generate_predictor(time_window = 43,
                         V_43 = V) %>%
      generate_predictor(time_window = 49,
                         V_49 = V) %>%
      
      
      # Generate the fitted weights for the synthetic control
      generate_weights(optimization_window = 38:50, # time to use in the optimization task
      ) %>%
      
      # Generate the synthetic control
      generate_control()
    
    
    p_value_rank <- V_out_predictors %>% grab_signficance()
    
    if (p_value_rank$unit_name[1]== 1){
      diagnostic$`p_value_<=0.05`[(4*r)] = 1
      diagnostic$`p_value_<=0.1`[(4*r)] = 1
    } else if(p_value_rank$unit_name[2]== 1){
      diagnostic$`p_value_<=0.1`[(4*r)] = 1
    }
    
    synth_control <- V_out_predictors %>% grab_synthetic_control()
    
    synth_control <- synth_control[synth_control$time_unit >= 38 & synth_control$time_unit <= 50,]
    
    diagnostic$`R^2`[(4*r)] <- calculate_R2(synth_control)
  }, error = function(e) {
    # Handle the error by setting values to NA
    diagnostic$`R^2`[(4 * r)] <- NA
    diagnostic$`p_value_<=0.05`[(4 * r)] <- NA
    diagnostic$`p_value_<=0.1`[(4 * r)] <- NA
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  # Check if the current iteration is a multiple of 50
  if (r %% 50 == 0) {
    # Increment the iteration count
    iteration_count <- iteration_count + 1
    
    # Construct the filename based on the current range of r
    filename <- paste0("tcc_diagnosticc_4001_to_", (iteration_count * 50 + 4001), ".csv")
    
    # Save the diagnostic data frame to CSV
    write.csv(diagnostic, filename, row.names = FALSE)
  }
}
