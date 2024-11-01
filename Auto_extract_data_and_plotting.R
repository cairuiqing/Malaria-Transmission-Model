library(ggplot2)
library(dplyr)
library(magrittr)
library(stringr)
library(rcartocolor)
library(see)
library(cowplot)
library(RColorBrewer)
extract_simulation_data <- function(file_name_human, file_name_mosquito, file_name_symptom, 
                                    file_name_eir, file_name_location,
                                    biting_scenario, parameter, value) {
  hap_age <- readRDS(file_name_human)
  eir <- readRDS(file_name_eir)
  symptom <- readRDS(file_name_symptom)
  location_data <- readRDS(file_name_location)
  
  num_persons <- dim(hap_age)[1]
  num_days <- 730
  num_simulations <- 1
  num_haplotypes <- dim(hap_age)[2] / num_days
  
  moi_array <- array(0, c(num_persons, num_days, num_simulations))
  symptom_array <- array(0, c(num_persons, num_days, num_simulations))
  
  # for (i in 1:num_simulations) {
  #   for (j in 1:num_days) {
  #     for (k in 1:num_persons) {
  #       moi_array[k, j, i] <- length(which(hap_age[k, (((j - 1) * num_haplotypes) + 1):(j * num_haplotypes), i] > 7))
  #       symptom_array[k, j, i] <- ifelse(is.na(symptom[k, j, i]), 0, symptom[k, j, i])
  #     }
  #   }
  # }
  # 
  # Verify dimensions before the loop
  print(dim(hap_age))
  print(num_persons)
  print(num_haplotypes)
  print(num_days)
  print(num_simulations)
  
  for (i in 1:num_simulations) {
    for (j in 1:num_days) {
      for (k in 1:num_persons) {
        
        
        moi_array[k, j, i] <- length(which(hap_age[k, (((j - 1) * num_haplotypes) + 1):(j * num_haplotypes), i] > 7))
        symptom_array[k, j, i] <- ifelse(is.na(symptom[k, j, i]), 0, symptom[k, j, i])
      }
    }
  }
  
  long_moi_symptom_status <- data.frame(
    MOI = as.vector(moi_array),
    Symptoms = as.vector(symptom_array),
    Location = as.vector(location_data),
    Simulation = rep(1:num_simulations, each = num_persons * num_days),
    Day = rep(1:num_days, each = num_persons, times = num_simulations),
    Person = rep(1:num_persons, times = num_days * num_simulations), 
    months = rep(c(rep(1:24, each = num_persons * 30), rep(25, each = num_persons * 10)), times = num_simulations)
  )
  
  monthly_moi <- long_moi_symptom_status %>%
    filter(Day %in% seq(1, 730, by = 30) | Symptoms == 1)
  
  monthly_moi_and_pos <- monthly_moi %>%
    mutate(pcr_pos = ifelse(MOI > 0, 1, 0))
 
  # Prevalence Calculation for the General Population
  monthly_prev <- monthly_moi_and_pos %>%
    group_by(months, Simulation) %>%
    summarize(
      num_pos_in_month = length(which(pcr_pos == 1)), 
      proportion_pos_in_month = length(which(pcr_pos == 1)) / n(), 
      number_sampled = n(),
      .groups = 'drop'
    )
  
  monthly_prev_summary <- monthly_prev %>%
    group_by(months) %>%
    summarize(
      median_num_pos = median(num_pos_in_month), 
      median_proportion_pos = median(proportion_pos_in_month),
      Biting = biting_scenario, 
      Parameter = parameter, 
      Value = value,
      .groups = 'drop'
    )
  
  # Prevalence Calculation by Location
  monthly_prev_location <- monthly_moi_and_pos %>%
    group_by(months, Location, Simulation) %>%
    summarize(
      num_pos_in_month = length(which(pcr_pos == 1)), 
      proportion_pos_in_month = length(which(pcr_pos == 1)) / n(), 
      number_sampled = n(),
      .groups = 'drop'
    )
  monthly_prev_location_summary <- monthly_prev_location %>%
    group_by(months, Location) %>%
    summarize(
      median_num_pos = median(num_pos_in_month), 
      median_proportion_pos = median(proportion_pos_in_month),
      Biting = biting_scenario, 
      Parameter = parameter, 
      Value = value,
      .groups = 'drop'
    )
  
  overall_MOI_sampled <- monthly_moi_and_pos %>%
    filter(months > 12) %>%
    group_by(Simulation) %>%
    reframe(Y = MOI, Biting = biting_scenario, Parameter = parameter, Value = value, Metric = "Human MOI") %>%
    filter(Y != 0)
  
  eir_length <- length(as.vector(eir))
  overall_eir <- data.frame(
    Y = as.vector(eir),
    Simulation = rep(1:num_simulations, each = eir_length / num_simulations),
    Biting = biting_scenario,
    Parameter = parameter,
    Value = value,
    Metric = "EIR"
  )
  
  mosquito_moi <- readRDS(file_name_mosquito)
  
  mosquito_sims_moi <- data.frame(
    Y = as.vector(mosquito_moi),
    Simulation = rep(1:num_simulations, each = length(as.vector(mosquito_moi)) / num_simulations),
    Biting = biting_scenario,
    Parameter = parameter,
    Value = value,
    Metric = "Mosquito MOI"
  )
  
  overall_mos_moi <- mosquito_sims_moi %>%
    filter(Y != 0 & !is.na(Y))
  
  sim_metric_distribution_df <- rbind(overall_MOI_sampled, overall_eir, overall_mos_moi)
  
  hum_MOI_in_range <- overall_MOI_sampled %>%
    mutate(in_range = ifelse(Y <= 6 & Y >= 1, 1, 0))
  
  hum_MOI_in_range_df <- data.frame(
    proportion_in_range = mean(hum_MOI_in_range$in_range), 
    Metric = "Human MOI",
    Biting = biting_scenario, 
    Parameter = parameter,
    Value = value
  )
  
  eir_in_range <- overall_eir %>%
    mutate(in_range = ifelse(Y < 3, 1, 0))
  
  eir_in_range_df <- data.frame(
    proportion_in_range = mean(eir_in_range$in_range), 
    Metric = "EIR",
    Biting = biting_scenario, 
    Parameter = parameter,
    Value = value
  )
  
  mos_MOI_in_range <- overall_mos_moi %>%
    mutate(in_range = ifelse(Y >= 4 & Y <= 9, 1, 0))
  
  mos_MOI_in_range_df <- data.frame(
    proportion_in_range = mean(mos_MOI_in_range$in_range), 
    Metric = "Mosquito MOI",
    Biting = biting_scenario, 
    Parameter = parameter,
    Value = value
  )
  
  in_range_df <- rbind(hum_MOI_in_range_df, eir_in_range_df, mos_MOI_in_range_df)
  
  pr_eir_df <- data.frame(
    EIR = median(as.vector(eir)),
    Prev = median(monthly_prev_summary$median_proportion_pos),
    Biting = biting_scenario,
    Parameter = parameter,
    Value = value
  )
  
  return(list(monthly_prev_summary, monthly_prev_location_summary, sim_metric_distribution_df, in_range_df, pr_eir_df))
}

# Define the scenario name and folder path
scenario_name <- "Prop01Prob03EvenTrips"
folder_path <- scenario_name

# Extract data using the scenario name
scenario_data <- extract_simulation_data(
  file.path(folder_path, paste0("haplotype_age_", scenario_name)),
  file.path(folder_path, paste0("mosquito_MOI_", scenario_name)),
  file.path(folder_path, paste0("symptom_status_", scenario_name)),
  file.path(folder_path, paste0("eir_", scenario_name)),
  file.path(folder_path, paste0("location_", scenario_name)),
  scenario_name, "UnevenTrips", 0.1
)

saveRDS(scenario_data, file = file.path(folder_path, paste0(scenario_name, "_plotfile.rds")))

monthly_prev_summary <- scenario_data[[1]]
monthly_prev_by_loc_summary <- scenario_data[[2]]
sim_metric_distribution_df <- scenario_data[[3]]
in_range_df <- scenario_data[[4]]
pr_eir_df <- scenario_data[[5]]

# Plotting Monthly Prevalence for General Population
ggplot(monthly_prev_summary, aes(x = months, y = median_proportion_pos)) +
  geom_line() +
  geom_point() +
  labs(title = "Monthly Prevalence for General Population",
       x = "Months",
       y = "Median Proportion Positive") + 
  theme_classic()

# Plotting Monthly Prevalence by Location
ggplot(monthly_prev_by_loc_summary, aes(x = months, y = median_proportion_pos, color = as.factor(Location))) +
  geom_line() +
  geom_point() +
  labs(title = "Monthly Prevalence by Location",
       x = "Months",
       y = "Median Proportion Positive",
       color = "Location") + 
  theme_classic() +
  scale_fill_brewer(palette = "Set3")

# Plotting Violin Plot for Simulation Metrics
sim_metric_distribution_df$Metric <- factor(sim_metric_distribution_df$Metric, levels = c("Human MOI", "EIR", "Mosquito MOI"))

ggplot(sim_metric_distribution_df, aes(x = Metric, y = Y, fill = Metric)) +
  geom_violin(trim = FALSE, scale='width') +
  geom_boxplot(color='black', outlier.shape = NA, width=0.1) +
  labs(title = "Distribution of Simulation Metrics for Even Trips",
       x = "Metric",
       y = "Y",
       fill = "Metric") +
  theme_classic() +
  scale_fill_brewer(palette = "Set3") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Extract and Plot Prevalence vs MOI Relationship
human_MOI_data <- sim_metric_distribution_df %>%
  filter(Metric == "Human MOI") %>%
  mutate(MOI_Type = "Human MOI")

mosquito_MOI_data <- sim_metric_distribution_df %>%
  filter(Metric == "Mosquito MOI") %>%
  mutate(MOI_Type = "Mosquito MOI")

MOI_combined <- rbind(human_MOI_data, mosquito_MOI_data)

prevalence_MOI_df <- monthly_prev_summary %>%
  mutate(MOI_Human = human_MOI_data$Y[1:nrow(monthly_prev_summary)],
         MOI_Mosquito = mosquito_MOI_data$Y[1:nrow(monthly_prev_summary)])

ggplot(prevalence_MOI_df) +
  geom_point(aes(x = MOI_Human, y = median_proportion_pos, color = "Human MOI"), size = 2) +
  geom_point(aes(x = MOI_Mosquito, y = median_proportion_pos, color = "Mosquito MOI"), size = 2) +
  labs(title = "Prevalence vs. MOI (Human and Mosquito)",
       x = "MOI",
       y = "Median Proportion Positive") +
  theme_classic() +
  scale_color_manual(name = "MOI Type", values = c("Human MOI" = "#9966FF", "Mosquito MOI" = "#FFCC00")) +
  theme(plot.title = element_text(hjust = 0.5))