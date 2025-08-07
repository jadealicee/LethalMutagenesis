#load
library(tidyverse)
library(dplyr)
library(ggplot2)

#set up file path and simulation parameters
master_csv_path <- "~/Desktop/Simulation/N=K=1000/master_csvs/master_N1000.csv"
sim_length <- 500 #number of gens in the sim
reporting_frequency <- 1 #report every generation

num_expected_rows_full_sim <- sim_length
print(paste("Expected number of rows for a full simulation (Gen 1 to Gen", sim_length, "):", num_expected_rows_full_sim))

#read master CSV file into a dataframe
df_master <- read_csv(master_csv_path)

#process raw data to calculate survival status and extinction generation for each run
df_processed <- df_master %>%
  #group the data by the unique identifiers of each simulation run (L, R, U, and seed)
  group_by(L, R, U, seed) %>%
  #summarise each group to get key metrics
  summarise(
    #count the number of rows reported for each run
    num_rows_reported = n(),
    #get the population size from the last row (last generation)
    final_population_size = last(psiz),
    .groups = 'drop'
  ) %>%
  #add new columns based on the summarised data
  mutate(
    #create factor variable for 'Survived' or 'Died'
    survival_status = factor(
      ifelse(
        #run 'Survived' if final population > 0 AND it ran for the full simulation length
        final_population_size > 0 & num_rows_reported == num_expected_rows_full_sim,
        "Survived",
        "Died"
      )
    ),
    #determine extinction generation
    extinction_generation = ifelse(
      #if the run survived, extinction generation is full sim length of 500
      survival_status == "Survived",
      sim_length,
      #if it died, it's the generation before the last one reported
      num_rows_reported - 1
    )
  )

###Extinction Generation Plot###
plot_extinction_generation <- ggplot(df_processed, aes(x = as.factor(L), y = extinction_generation, color = survival_status, shape = survival_status)) +
  #jitter to visualise overlapping data points
  geom_jitter(width = 0.15, size = 3, alpha = 0.7) +
  #set colors and shapes for 'Died' and 'Survived' statuses
  scale_color_manual(values = c("Died" = "red", "Survived" = "blue")) +
  scale_shape_manual(values = c("Died" = 4, "Survived" = 16)) +
  #plot labels and title
  labs(
    title = "Impact of Lethality Rate (L) on Virus Population Survival",
    x = "Lethality Rate (L)",
    y = "Generation of Extinction (or Survival at Max Sim Length)",
    color = "Run Outcome",
    shape = "Run Outcome"
  ) +
  #add horizontal dashed line for the maximum sim length
  geom_hline(yintercept = sim_length, linetype = "dashed", color = "darkgreen", alpha = 0.6) +
  #annotation for max sim length
  annotate("text", x = max(df_processed$L) + 0.1, y = sim_length + 10, label = "Max Simulation Length", hjust = 1, vjust = 0, color = "darkgreen", size = 3.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  #facet the plot by both R and U
  facet_grid(R ~ U, labeller = "label_both")

print(plot_extinction_generation)

###Survival Rate Plot###

#calculate survival rates based on processed data
survival_rates <- df_processed %>%
  #group the data by L, R, and U to calculate a percentage for each combination
  group_by(L, R, U) %>%
  summarise(
    #count the number of runs that survived
    num_survived = sum(survival_status == "Survived"),
    #count total number of runs in the group
    total_runs = n(),
    #calculate survival percentage
    survival_percentage = (num_survived / total_runs) * 100,
    .groups = 'drop'
  )

print("Calculated survival rates per L, R, and U:")
print(survival_rates)

plot_survival_rate <- ggplot(survival_rates, aes(x = as.factor(L), y = survival_percentage)) +
  geom_col(fill = "steelblue", width = 0.7) +
  #add text labels on top of columns showing the exact percentage
  geom_text(aes(label = sprintf("%.1f%%", survival_percentage)), vjust = -0.5, size = 3.5) +
  #plot labels and title
  labs(
    title = "Percentage of Runs Surviving by Lethal Rate (L)",
    x = "Lethal Rate (L)",
    y = "Survival Percentage (%)"
  ) +
  #format the y-axis labels as percentages
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  #facet the plot by both R and U
  facet_grid(R ~ U, labeller = "label_both")

print(plot_survival_rate)