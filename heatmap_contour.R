#load
library(tidyverse)
library(dplyr)
library(ggplot2)

#edit path to where master CSV file is located
master_csv_path <- "~/Desktop/Simulation/N=K=1000/master_csvs/master_N1000.csv"
sim_length <- 500 # number of generations 

#500 total
num_expected_rows_full_sim <- sim_length

#read master CSV file into a dataframe
df_master <- read_csv(master_csv_path)

#process raw data to calculate survival status and extinction generation for each run
df_processed <- df_master %>%
  #group data by unique identifiers of each sim run
  group_by(L, R, U, seed) %>%
  #summarise each group
  summarise(
    num_rows_reported = n(),
    final_population_size = last(psiz),
    .groups = 'drop'
  ) %>%
  #add new columns based on summarised data
  mutate(
    survival_status = factor(
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

#calculate survival rates based on the processed data for plotting
survival_rates <- df_processed %>%
  #group data by L, R, and U to calculate a percentage for each combo
  group_by(L, R, U) %>%
  summarise(
    num_survived = sum(survival_status == "Survived"),
    total_runs = n(),
    survival_percentage = (num_survived / total_runs) * 100,
    .groups = 'drop'
  )

###Heatmap of Survival Percentage###
plot_heatmap <- ggplot(survival_rates, aes(x = as.factor(L), y = as.factor(R), fill = survival_percentage)) +
  geom_tile() +
  #colour scale for survival percentage
  scale_fill_viridis_c(name = "Survival (%)", option = "D", direction = -1) +
  #plot labels and title
  labs(
    title = "Heatmap of Survival Percentage by Lethal (L) and Reproductive (R) Rates",
    x = "Lethal Rate (L)",
    y = "Reproductive Rate (R)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  #facet plot by mutation rate (U)
  facet_wrap(~ U)

print(plot_heatmap)

###############################################################

###Contour Plot of Survival Percentage###
plot_contour <- ggplot(survival_rates, aes(x = L, y = R, z = survival_percentage)) +
  geom_contour_filled(aes(fill = after_stat(level)), alpha = 0.7) +
  geom_point(color = "black", size = 2) +
  #colour scale
  scale_fill_viridis_d(name = "Survival (%)") +
  #plot labels and title
  labs(
    title = "Contour Plot of Survival Percentage by Lethal (L) and Reproductive (R) Rates",
    x = "Lethal Rate (L)",
    y = "Reproductive Rate (R)"
  ) +
  #y-axis breaks to be the unique R values (increments of 1)
  scale_y_continuous(breaks = unique(survival_rates$R)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  #facet the plot by mutation rate (U)
  facet_wrap(~ U)

print(plot_contour)
