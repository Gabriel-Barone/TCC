library(dplyr)
library(readxl)

setwd("C:/Users/gpime/OneDrive/Área de Trabalho/TCC")


# Read the CSV file into a data frame
df <- read.csv("DIAGNOSTIC_FILE_NAME")
df$R.2[is.na(df$R.2)] <- 1
df <- df[complete.cases(df), ]

# Print the updated data frame
print(df)

filtered_df <- df %>%
  filter(`p_value_<=0.1` == 1)

filtered_df_R2_85 <- filtered_df %>%
  filter(`R^2` >= 0.85)

df_R2_85 <- df %>%
  filter(`R^2` >= 0.85)

filtered_df_R2_95 <- filtered_df %>%
  filter(`R^2` >= 0.95)

df_R2_95 <- df %>%
  filter(`R^2`>= 0.95)

p_value_0.05_raw <- length(unique(filtered_df$series_id))/length(unique(df$series_id))
p_value_0.05_85 <- length(unique(filtered_df_R2_85$series_id))/length(unique(df_R2_85$series_id))
p_value_0.05_95 <- length(unique(filtered_df_R2_95$series_id))/length(unique(df_R2_95$series_id))

#------------noF-----------------------

df_noF <- df %>%
  filter(synth_control_type != "F_predictors")

filtered_df_noF <- filtered_df %>%
  filter(synth_control_type != "F_predictors")

filtered_df_R2_85_noF <- filtered_df_noF %>%
  filter(`R^2` >= 0.85)

df_R2_85_noF <- df_noF %>%
  filter(`R^2` >= 0.85)

filtered_df_R2_95_noF <- filtered_df_noF %>%
  filter(`R^2` >= 0.95)

df_R2_95_noF <- df_noF %>%
  filter(`R^2` >= 0.95)

p_value_0.05_raw_noF <- length(unique(filtered_df_noF$series_id))/length(unique(df_noF$series_id))
p_value_0.05_85_noF <- length(unique(filtered_df_R2_85_noF$series_id))/length(unique(df_R2_85_noF$series_id))
p_value_0.05_95_noF <- length(unique(filtered_df_R2_95_noF$series_id))/length(unique(df_R2_95_noF$series_id))

#----------noA------------------------------------------------------

df_noA <- df %>%
  filter(synth_control_type != "ALL_predictors")

filtered_df_noA <- filtered_df %>%
  filter(synth_control_type != "ALL_predictors")

filtered_df_R2_85_noA <- filtered_df_noA %>%
  filter(`R^2` >= 0.85)

df_R2_85_noA <- df_noA %>%
  filter(`R^2` >= 0.85)

filtered_df_R2_95_noA <- filtered_df_noA %>%
  filter(`R^2` >= 0.95)

df_R2_95_noA <- df_noA %>%
  filter(`R^2` >= 0.95)

p_value_0.05_raw_noA <- length(unique(filtered_df_noA$series_id))/length(unique(df_noA$series_id))
p_value_0.05_85_noA <- length(unique(filtered_df_R2_85_noA$series_id))/length(unique(df_R2_85_noA$series_id))
p_value_0.05_95_noA <- length(unique(filtered_df_R2_95_noA$series_id))/length(unique(df_R2_95_noA$series_id))

#----------noN------------------------------------------------------

df_noN <- df %>%
  filter(synth_control_type != "NO_predictors")

filtered_df_noN <- filtered_df %>%
  filter(synth_control_type != "NO_predictors")

filtered_df_R2_85_noN <- filtered_df_noN %>%
  filter(`R^2` >= 0.85)

df_R2_85_noN <- df_noN %>%
  filter(`R^2` >= 0.85)

filtered_df_R2_95_noN <- filtered_df_noN %>%
  filter(`R^2` >= 0.95)

df_R2_95_noN <- df_noN %>%
  filter(`R^2` >= 0.95)

p_value_0.05_raw_noN <- length(unique(filtered_df_noN$series_id))/length(unique(df_noN$series_id))
p_value_0.05_85_noN <- length(unique(filtered_df_R2_85_noN$series_id))/length(unique(df_R2_85_noN$series_id))
p_value_0.05_95_noN <- length(unique(filtered_df_R2_95_noN$series_id))/length(unique(df_R2_95_noN$series_id))


#----------noFA------------------------------------------------------

df_T <- df %>%
  filter(synth_control_type == "T_predictors" | synth_control_type == "NO_predictors" )

filtered_df_T <- filtered_df %>%
  filter(synth_control_type == "T_predictors" | synth_control_type == "NO_predictors")

filtered_df_R2_85_T <- filtered_df_T %>%
  filter(`R^2` >= 0.85)

df_R2_85_T <- df_T %>%
  filter(`R^2` >= 0.85)

filtered_df_R2_95_T <- filtered_df_T %>%
  filter(`R^2` >= 0.95)

df_R2_95_T <- df_T %>%
  filter(`R^2` >= 0.95)

p_value_0.05_raw_T <- length(unique(filtered_df_T$series_id))/length(unique(df_T$series_id))
p_value_0.05_85_T <- length(unique(filtered_df_R2_85_T$series_id))/length(unique(df_R2_85_T$series_id))
p_value_0.05_95_T <- length(unique(filtered_df_R2_95_T$series_id))/length(unique(df_R2_95_T$series_id))



























filtered_df <- df %>%
  filter(`p_value_..0.1` == 1)

filtered_df_R2_85 <- filtered_df %>%
  filter(`R.2` >= 0.85)

df_R2_85 <- df %>%
  filter(`R.2` >= 0.85)

filtered_df_R2_95 <- filtered_df %>%
  filter(`R.2` >= 0.95)

df_R2_95 <- df %>%
  filter(`R.2`>= 0.95)

p_value_0.05_raw <- length(unique(filtered_df$series_id))/length(unique(df$series_id))
p_value_0.05_85 <- length(unique(filtered_df_R2_85$series_id))/length(unique(df_R2_85$series_id))
p_value_0.05_95 <- length(unique(filtered_df_R2_95$series_id))/length(unique(df_R2_95$series_id))

#------------noF-----------------------

df_noF <- df %>%
  filter(synth_control_type != "F_predictors")

filtered_df_noF <- filtered_df %>%
  filter(synth_control_type != "F_predictors")

filtered_df_R2_85_noF <- filtered_df_noF %>%
  filter(`R.2` >= 0.85)

df_R2_85_noF <- df_noF %>%
  filter(`R.2` >= 0.85)

filtered_df_R2_95_noF <- filtered_df_noF %>%
  filter(`R.2` >= 0.95)

df_R2_95_noF <- df_noF %>%
  filter(`R.2` >= 0.95)

p_value_0.05_raw_noF <- length(unique(filtered_df_noF$series_id))/length(unique(df_noF$series_id))
p_value_0.05_85_noF <- length(unique(filtered_df_R2_85_noF$series_id))/length(unique(df_R2_85_noF$series_id))
p_value_0.05_95_noF <- length(unique(filtered_df_R2_95_noF$series_id))/length(unique(df_R2_95_noF$series_id))

#----------noA------------------------------------------------------

df_noA <- df %>%
  filter(synth_control_type != "ALL_predictors")

filtered_df_noA <- filtered_df %>%
  filter(synth_control_type != "ALL_predictors")

filtered_df_R2_85_noA <- filtered_df_noA %>%
  filter(`R.2` >= 0.85)

df_R2_85_noA <- df_noA %>%
  filter(`R.2` >= 0.85)

filtered_df_R2_95_noA <- filtered_df_noA %>%
  filter(`R.2` >= 0.95)

df_R2_95_noA <- df_noA %>%
  filter(`R.2` >= 0.95)

p_value_0.05_raw_noA <- length(unique(filtered_df_noA$series_id))/length(unique(df_noA$series_id))
p_value_0.05_85_noA <- length(unique(filtered_df_R2_85_noA$series_id))/length(unique(df_R2_85_noA$series_id))
p_value_0.05_95_noA <- length(unique(filtered_df_R2_95_noA$series_id))/length(unique(df_R2_95_noA$series_id))

#----------noN------------------------------------------------------

df_noN <- df %>%
  filter(synth_control_type != "NO_predictors")

filtered_df_noN <- filtered_df %>%
  filter(synth_control_type != "NO_predictors")

filtered_df_R2_85_noN <- filtered_df_noN %>%
  filter(`R.2` >= 0.85)

df_R2_85_noN <- df_noN %>%
  filter(`R.2` >= 0.85)

filtered_df_R2_95_noN <- filtered_df_noN %>%
  filter(`R.2` >= 0.95)

df_R2_95_noN <- df_noN %>%
  filter(`R.2` >= 0.95)

p_value_0.05_raw_noN <- length(unique(filtered_df_noN$series_id))/length(unique(df_noN$series_id))
p_value_0.05_85_noN <- length(unique(filtered_df_R2_85_noN$series_id))/length(unique(df_R2_85_noN$series_id))
p_value_0.05_95_noN <- length(unique(filtered_df_R2_95_noN$series_id))/length(unique(df_R2_95_noN$series_id))


#----------noFA------------------------------------------------------

df_T <- df %>%
  filter(synth_control_type == "T_predictors" | synth_control_type == "NO_predictors" )

filtered_df_T <- filtered_df %>%
  filter(synth_control_type == "T_predictors" | synth_control_type == "NO_predictors")

filtered_df_R2_85_T <- filtered_df_T %>%
  filter(`R.2` >= 0.85)

df_R2_85_T <- df_T %>%
  filter(`R.2` >= 0.85)

filtered_df_R2_95_T <- filtered_df_T %>%
  filter(`R.2` >= 0.95)

df_R2_95_T <- df_T %>%
  filter(`R.2` >= 0.95)

p_value_0.05_raw_T <- length(unique(filtered_df_T$series_id))/length(unique(df_T$series_id))
p_value_0.05_85_T <- length(unique(filtered_df_R2_85_T$series_id))/length(unique(df_R2_85_T$series_id))
p_value_0.05_95_T <- length(unique(filtered_df_R2_95_T$series_id))/length(unique(df_R2_95_T$series_id))

