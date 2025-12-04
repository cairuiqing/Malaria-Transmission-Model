library(dplyr)

compute_importation_metrics_simple <- function(
    scenario_name,
    base_dir = ".",
    burnin_day = 365,
    write_outputs = TRUE
) {
  # Path to the transmission-chain log
  chain_path <- file.path(
    base_dir,
    scenario_name,
    paste0("transmission_chain_log_", scenario_name, ".csv")
  )
  
  if (!file.exists(chain_path)) {
    stop("File not found: ", chain_path)
  }
  
  chain <- read.csv(chain_path, stringsAsFactors = FALSE)
  
  if (!nrow(chain)) {
    stop("Transmission chain log is empty: ", chain_path)
  }
  
  # Keep only post-burn-in days (year 2 by default)
  chain <- chain %>%
    dplyr::filter(Day > burnin_day)
  
  if (!nrow(chain)) {
    stop("No transmission events after burn-in day = ", burnin_day)
  }
  
  # Define which SourceType categories count as "traveler involvement"
  traveler_types <- c(
    "Traveler_Recent",
    "Traveler_LongTerm",
    "Returnee_Imported",
    "Returnee_Recent"
  )
  
  chain <- chain %>%
    mutate(
      TransmissionType   = as.character(TransmissionType),
      SourceType         = as.character(SourceType),
      SourceWasTraveler  = SourceType %in% traveler_types
    )
  
  ## 1. Overall metrics -------------------------------------------
  overall <- chain %>%
    summarise(
      n_events          = n(),
      n_imported        = sum(TransmissionType == "Imported", na.rm = TRUE),
      importation_rate  = ifelse(n_events > 0, n_imported / n_events, NA_real_),
      n_travel_source   = sum(SourceWasTraveler, na.rm = TRUE),
      p_source_traveler = ifelse(n_events > 0, n_travel_source / n_events, NA_real_)
    )
  
  ## 2. By target location ----------------------------------------
  by_target <- chain %>%
    group_by(TargetAddress) %>%
    summarise(
      n_events          = n(),
      n_imported        = sum(TransmissionType == "Imported", na.rm = TRUE),
      importation_rate  = n_imported / n_events,
      n_travel_source   = sum(SourceWasTraveler, na.rm = TRUE),
      p_source_traveler = n_travel_source / n_events,
      .groups           = "drop"
    )
  
  ## 3. By SourceType (Local vs traveler categories) --------------
  by_source_type <- chain %>%
    group_by(SourceType) %>%
    summarise(
      n_events          = n(),
      n_imported        = sum(TransmissionType == "Imported", na.rm = TRUE),
      importation_rate  = n_imported / n_events,
      .groups           = "drop"
    ) %>%
    mutate(
      frac_of_all_events = n_events / sum(n_events)
    )
  
  ## 4. By simulation ---------------------------------------------
  by_sim <- chain %>%
    group_by(Simulation) %>%
    summarise(
      n_events          = n(),
      n_imported        = sum(TransmissionType == "Imported", na.rm = TRUE),
      importation_rate  = n_imported / n_events,
      n_travel_source   = sum(SourceWasTraveler, na.rm = TRUE),
      p_source_traveler = n_travel_source / n_events,
      .groups           = "drop"
    )
  
  ## 5. Time series: overall by day -------------------------------
  by_day <- chain %>%
    group_by(Day) %>%
    summarise(
      n_events          = n(),
      n_imported        = sum(TransmissionType == "Imported", na.rm = TRUE),
      importation_rate  = n_imported / n_events,
      n_travel_source   = sum(SourceWasTraveler, na.rm = TRUE),
      p_source_traveler = n_travel_source / n_events,
      .groups           = "drop"
    )
  
  ## 6. Time series: by day and target location -------------------
  by_day_target <- chain %>%
    group_by(TargetAddress, Day) %>%
    summarise(
      n_events          = n(),
      n_imported        = sum(TransmissionType == "Imported", na.rm = TRUE),
      importation_rate  = n_imported / n_events,
      n_travel_source   = sum(SourceWasTraveler, na.rm = TRUE),
      p_source_traveler = n_travel_source / n_events,
      .groups           = "drop"
    )
  
  # Optionally write summaries out to CSV in the scenario folder
  if (write_outputs) {
    out_dir <- file.path(base_dir, scenario_name)
    
    write.csv(overall,
              file = file.path(out_dir, "importation_overall_summary.csv"),
              row.names = FALSE)
    write.csv(by_target,
              file = file.path(out_dir, "importation_by_target.csv"),
              row.names = FALSE)
    write.csv(by_source_type,
              file = file.path(out_dir, "importation_by_source_type.csv"),
              row.names = FALSE)
    write.csv(by_sim,
              file = file.path(out_dir, "importation_by_simulation.csv"),
              row.names = FALSE)
    write.csv(by_day,
              file = file.path(out_dir, "importation_by_day.csv"),
              row.names = FALSE)
    write.csv(by_day_target,
              file = file.path(out_dir, "importation_by_day_target.csv"),
              row.names = FALSE)
  }
  
  # Return everything as a list for direct use
  list(
    overall        = overall,
    by_target      = by_target,
    by_source_type = by_source_type,
    by_sim         = by_sim,
    by_day         = by_day,
    by_day_target  = by_day_target
  )
}



res <- compute_importation_metrics("tracking test",
                                   base_dir = "H:/Johns Hopkins Bloomberg School of Public Health/Research Amy/Malaria Sequencing Project/multi-job_idd_server_submission/11-19-2025_after_debug")

res$overall
res$by_target
res$by_od