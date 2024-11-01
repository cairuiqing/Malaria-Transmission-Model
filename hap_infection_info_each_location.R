library(ggplot2)
library(dplyr)
library(magrittr)
library(stringr)
library(rcartocolor)
library(see)
library(cowplot)
library(RColorBrewer)

# extract_imported_haplotypes_from_model_output <- function(file_name_human, file_name_location, file_name_mosquito, file_name_eir, num_simulations, num_days) {
#   
#   hap_age <- readRDS(file_name_human)
#   location_data <- readRDS(file_name_location)
#   mosquito_MOI <- readRDS(file_name_mosquito)
#   eir_data <- readRDS(file_name_eir)
#   
#   num_persons <- dim(location_data)[1]
#   num_haplotypes <- dim(hap_age)[2] / num_days
#   
#   imported_haplotypes <- list()  
#   
#   for (sim in 1:num_simulations) {
#     for (day in 2:num_days) {
#       for (person in 1:num_persons) {
#         
#         current_location <- location_data[person, day, sim]
#         previous_location <- location_data[person, day - 1, sim]
#         
#         # Check if the person moved to a new location
#         if (current_location != previous_location) {
#           # Identify the haplotypes this person is carrying
#           hap_carried <- which(hap_age[person, (((day - 1) * num_haplotypes) + 1):(day * num_haplotypes), sim] > 7)
#           
#           # For each haplotype the person is carrying, check for potential transmission to mosquitoes
#           for (hap in hap_carried) {
#             # Check if the person was bitten by mosquitoes in the new location using EIR (if there's transmission)
#             if (eir_data[sim, person] > 0) {  # Adjust to use sim index
#               # Track the mosquito that was infected by this haplotype and the new location
#               mos_infected <- mosquito_MOI[sim, person]  # Adjust to use sim index
#               if (!is.na(mos_infected) && mos_infected > 0) {
#                 # Record the haplotype import event
#                 imported_haplotypes[[current_location]] <- append(
#                   imported_haplotypes[[current_location]],
#                   list(hap = hap, time = day, person = person, simulation = sim)
#                 )
#               }
#             }
#           }
#         }
#       }
#     }
#   }
#   
#   # Convert the list of imported haplotypes into a dataframe
#   imported_haplotypes_df <- do.call(rbind, lapply(imported_haplotypes, function(x) {
#     data.frame(
#       Location = names(imported_haplotypes),
#       Haplotype = unlist(lapply(x, [[, "hap")),
#       Time = unlist(lapply(x, [[, "time")),
#       Person = unlist(lapply(x, [[, "person")),
#       Simulation = unlist(lapply(x, [[, "simulation"))
#       )
#     }))
#   
#   return(imported_haplotypes_df)
#   }
# 
# # Example of calling the function with file paths and number of simulations and days
# scenario_name <- "Prop01Prob03EvenTrips"
# folder_path <- scenario_name
# num_simulations <- 1  # Replace with the actual number of simulations
# num_days <- 730  # Replace with the actual number of days
# 
# scenario_imported_hap <- extract_imported_haplotypes_from_model_output(
#   file.path(folder_path, paste0("haplotype_age_", scenario_name)),
#   file.path(folder_path, paste0("location_", scenario_name)),
#   file.path(folder_path, paste0("mosquito_MOI_", scenario_name)),
#   file.path(folder_path, paste0("eir_", scenario_name)),
#   num_simulations,
#   num_days
# )
# 
# 
# print(scenario_imported_hap)

#############################################
##### Import hap info for whom traveled #####
#############################################
extract_imported_haplotypes_from_model_output <- function(file_name_human, file_name_location, num_simulations, num_days) {
  
  hap_age <- readRDS(file_name_human)
  location_data <- readRDS(file_name_location)
  
  num_persons <- dim(location_data)[1]
  num_haplotypes <- dim(hap_age)[2] / num_days
  
  imported_haplotypes <- list()  # To store the imported haplotypes for each location and time point
  
  for (sim in 1:num_simulations) {
    for (day in 2:num_days) {
      for (person in 1:num_persons) {
        # Get the current and previous day's locations for each simulation
        current_location <- location_data[person, day, sim]
        previous_location <- location_data[person, day - 1, sim]
        
        # Check if the person moved to a new location
        if (current_location != previous_location) { 
          # Identify the haplotypes this person is carrying
          hap_carried <- which(hap_age[person, (((day - 1) * num_haplotypes) + 1):(day * num_haplotypes), sim] > 7)
          
          if (length(hap_carried) > 0) {
            # Initialize the list for the current location if it doesn't exist yet
            if (is.null(imported_haplotypes[[as.character(current_location)]])) {
              imported_haplotypes[[as.character(current_location)]] <- data.frame(
                Haplotype = integer(),
                Time = integer(),
                Person = integer(),
                Simulation = integer()
              )
            }
            # Append the imported haplotypes directly to a dataframe
            imported_haplotypes[[as.character(current_location)]] <- rbind(
              imported_haplotypes[[as.character(current_location)]],
              data.frame(
                Haplotype = hap_carried,
                Time = day,
                Person = person,
                Simulation = sim
              )
            )
          }
        }
      }
    }
  }
  
  # Combine all location dataframes into one
  imported_haplotypes_df <- do.call(rbind, lapply(names(imported_haplotypes), function(location) {
    data.frame(
      Location_to = location,
      imported_haplotypes[[location]]
    )
  }))
  imported_haplotypes_df$Location_to <- as.numeric(imported_haplotypes_df$Location_to)
  return(imported_haplotypes_df)
}

# Example of calling the function with file paths and number of simulations and days
scenario_name <- "Prop01Prob03EvenTrips"
folder_path <- scenario_name
num_simulations <- 1  # Replace with the actual number of simulations
num_days <- 730  # Replace with the actual number of days

scenario_imported_hap <- extract_imported_haplotypes_from_model_output(
  file.path(folder_path, paste0("haplotype_age_", scenario_name)),
  file.path(folder_path, paste0("location_", scenario_name)),
  num_simulations,
  num_days
)

print(scenario_imported_hap)


haplotype_summary <- scenario_imported_hap %>%
  group_by(Location_to) %>%
  summarise(Haplotype_Count = n())
haplotype_summary$Location <- factor(haplotype_summary$Location_to, levels = as.character(1:10))

ggplot(haplotype_summary, aes(x = Location_to, y = Haplotype_Count)) +
  geom_bar(stat = "identity", fill = "#33CC99", width = 0.75) +
  theme_classic() +
  labs(title = "Haplotypes Flowing into Locations Over Time",
       x = "Importation Location",
       y = "Counts of Imported Haplotypes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black")) +
  scale_x_discrete(limits = as.character(1:10)) 




#############################################
#############################################
################ debug ######################
#############################################
#############################################
extract_imported_haplotypes_from_model_output <- function(file_name_human, file_name_location, num_simulations, num_days) {
  
  hap_age <- readRDS(file_name_human)
  location_data <- readRDS(file_name_location)
  
  num_persons <- dim(location_data)[1]
  num_haplotypes <- dim(hap_age)[2] / num_days
  
  imported_haplotypes <- list()  # To store the imported haplotypes for each location and time point
  
  for (sim in 1:num_simulations) {
    for (day in 2:num_days) {
      for (person in 1:num_persons) {
        # Get the current and previous day's locations for each simulation
        current_location <- location_data[person, day, sim]
        previous_location <- location_data[person, day - 1, sim]
        
        # Check if the person moved to a new location
        if (current_location != previous_location) { 
          # Identify the haplotypes this person is carrying
          hap_carried <- which(hap_age[person, (((day - 1) * num_haplotypes) + 1):(day * num_haplotypes), sim] > 7)
          
          if (length(hap_carried) > 0) {
            # Initialize the list for the current location if it doesn't exist yet
            if (is.null(imported_haplotypes[[as.character(current_location)]])) {
              imported_haplotypes[[as.character(current_location)]] <- list()
            }
            # Append the imported haplotypes
            imported_haplotypes[[as.character(current_location)]] <- append(
              imported_haplotypes[[as.character(current_location)]],
              list(hap = hap_carried, time = day, person = person, simulation = sim)
            )
          }
        }
      }
    }
  }
  
  # Convert the list of imported haplotypes into a dataframe
  imported_haplotypes_df <- do.call(rbind, lapply(names(imported_haplotypes), function(location) {
    data.frame(
      Location = location,
      Haplotype = unlist(lapply(imported_haplotypes[[location]], `[[`, "hap")),
      Time = unlist(lapply(imported_haplotypes[[location]], `[[`, "time")),
      Person = unlist(lapply(imported_haplotypes[[location]], `[[`, "person")),
      Simulation = unlist(lapply(imported_haplotypes[[location]], `[[`, "simulation"))
    )
  }))
  
  return(imported_haplotypes_df)
}

# Example of calling the function with file paths and number of simulations and days
scenario_name <- "Prop01Prob03EvenTrips"
folder_path <- scenario_name
num_simulations <- 1  # Replace with the actual number of simulations
num_days <- 730  # Replace with the actual number of days

scenario_imported_hap <- extract_imported_haplotypes_from_model_output(
  file.path(folder_path, paste0("haplotype_age_", scenario_name)),
  file.path(folder_path, paste0("location_", scenario_name)),
  num_simulations,
  num_days
)

# Display or save the dataframe of imported haplotypes
print(scenario_imported_hap)

