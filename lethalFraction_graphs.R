library(tidyverse)
library(dplyr)
library(ggplot2)

#edit path to where master csv file is
master_csv_path <- "~/Desktop/Simulation/N=K=1000/master_csvs/master_N1000.csv"
sim_length <- 500 #number of generations in the sim
reporting_frequency <- 1 #report every generation

#calculate the expected number of rows for a full simulation
num_expected_rows_full_sim <- sim_length
print(paste("Expected number of rows for a full simulation (Gen 0 to Gen", sim_length, "):", num_expected_rows_full_sim))

#read master csv file into dataframe
df_master <- read_csv(master_csv_path)

###data processing###

df_processed <- df_master %>%
  group_by(L, R, seed) %>%
  summarise(
    num_rows_reported = n(),
#get the population size from the last row (last generation)
    final_population_size = last(psiz),
    .groups = 'drop'
  ) %>%
  mutate(
#create a factor variable for 'Survived' or 'Died'
    survival_status = factor(
      ifelse(
        final_population_size > 0 & num_rows_reported == num_expected_rows_full_sim,
        "Survived",
        "Died"
      )
    ),
#determine the extinction generation
    extinction_generation = ifelse(
      survival_status == "Survived",
      sim_length,
      num_rows_reported - 1
    )
  )

#debugging info to make sure the csv has been processed correctly
print("Head of processed data (df_processed) with revised survival logic:")
print(head(df_processed))
print(df_processed %>% count(survival_status))
print(df_processed %>% select(L, R, seed, num_rows_reported, survival_status) %>% head())
print("Summary of survival_status counts:")
print(df_processed %>% count(survival_status))

###first plot - generation of extinction###

plot_extinction_generation <- ggplot(df_processed, aes(x = as.factor(L), y = extinction_generation, color = survival_status, shape = survival_status)) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Died" = "red", "Survived" = "blue")) +
  scale_shape_manual(values = c("Died" = 4, "Survived" = 16)) +
  labs(
    title = "Impact of Lethal Fraction (L) on Virus Population Survival",
    x = "Lethal Fraction (L)",
    y = "Generation of Extinction (or Survival at Max Sim Length)",
    color = "Run Outcome",
    shape = "Run Outcome"
  ) +
  geom_hline(yintercept = sim_length, linetype = "dashed", color = "darkgreen", alpha = 0.6) +
  annotate("text", x = max(df_processed$L) + 0.1, y = sim_length + 10, label = "Max Simulation Length", hjust = 1, vjust = 0, color = "darkgreen", size = 3.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  facet_wrap(~ R, ncol = 2)

print(plot_extinction_generation)

##################################################################################

###impact of lethal rate on survival plot

num_expected_rows_full_sim <- sim_length
print(paste("Expected number of rows for a full simulation (Gen 0 to Gen", sim_length, "):", num_expected_rows_full_sim))

df_master <- read_csv(master_csv_path)

print("Column names after loading CSV:")
print(colnames(df_master))
print("Head of df_master after loading CSV:")
print(head(df_master))


#can probably get rid of this processing chunk. compare to prior block
df_processed <- df_master %>%
  group_by(L, R, seed) %>%
  summarise(
    num_rows_reported = n(),
    final_population_size = last(psiz),
    .groups = 'drop'
  ) %>%
  mutate(
    survival_status = factor(
      #changed && to & after issues. seems to work
      ifelse(
        final_population_size > 0 & num_rows_reported == num_expected_rows_full_sim,
        "Survived",
        "Died"
      )
    ),
    extinction_generation = ifelse(
      survival_status == "Survived",
      sim_length,
      num_rows_reported - 1
    )
  )

#more debugging
print("--- Debugging df_processed Data ---")
print("Counts of runs where final_population_size > 0:")
print(df_processed %>% filter(final_population_size > 0) %>% count())

print("Counts of runs where num_rows_reported == num_expected_rows_full_sim:")
print(df_processed %>% filter(num_rows_reported == num_expected_rows_full_sim) %>% count())

print("Counts of runs where *both* conditions are met (i.e., 'Survived' before factor conversion):")
print(df_processed %>% filter(final_population_size > 0 & num_rows_reported == num_expected_rows_full_sim) %>% count())

print("Actual counts of 'Survived' vs 'Died' from df_processed:")
print(df_processed %>% count(survival_status))


survival_rates <- df_processed %>%
  group_by(L, R) %>%
  summarise(
    num_survived = sum(survival_status == "Survived"),
    total_runs = n(),
    survival_percentage = (num_survived / total_runs) * 100,
    .groups = 'drop'
  )

#make sure that survival is occurring
print("Calculated survival rates per L and R:")
print(survival_rates)

#plot
plot_survival_rate <- ggplot(survival_rates, aes(x = as.factor(L), y = survival_percentage)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", survival_percentage)), vjust = -0.5, size = 3.5) +
  labs(
    title = "Percentage of Runs Surviving by Lethal Rate (L)",
    x = "Lethal Rate (L)",
    y = "Survival Percentage (%)"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "none" 
  ) +
  facet_wrap(~ R, ncol = 2)
print(plot_survival_rate)


