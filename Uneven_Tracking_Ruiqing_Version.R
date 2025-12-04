library(dplyr)
library(extraDistr)
library(microbenchmark)
library(Rcpp)

sourceCpp("get_biting_status.cpp")

######################################################
######## Fixed Parameters and Re-used Functions ####
######################################################
n_p <- c(20, 10, 10)
n_m <- c(3000, 1500, 1500)

num_loc <- length(n_p)

# Day numbering, sample days, etc.
mos_sample_days <- seq(from = 365, to = 730, by = 7)
hum_sample_days <- seq(from = 365, to = 730, by = 30)

dry_days <- c(1:59, 305:424, 670:730)
rainy_days <- c(60:151, 425:516)
moderate_days <- c(152:304, 517:669)

# Read haplotype frequencies from CSV
hap_table_csv <- read.csv("haplotype_frequencies.csv")

human_freq_data <- hap_table_csv %>%
  filter(gene == "CSP" & species == "Human")

human_haps <- human_freq_data$haplotype_number
human_freq <- human_freq_data$frequency

mos_freq_data <- hap_table_csv %>%
  filter(gene == "CSP" & species == "Mosquito")

mos_haps <- mos_freq_data$haplotype_number
mos_freq <- mos_freq_data$frequency

mos_gt_df <- data.frame(haps = mos_haps, freq = mos_freq)

csp_total <- hap_table_csv %>%
  filter(gene == "CSP") %>%
  group_by(haplotype_number) %>%
  summarize(overall_frequency = sum(frequency)) %>%
  mutate(freq_category = ifelse(overall_frequency <= median(overall_frequency), 1,
                                ifelse(overall_frequency > quantile(overall_frequency, 0.75), 3, 2)))

gt_df <- data.frame(hap = csp_total$haplotype_number, 
                    freq = csp_total$overall_frequency, 
                    freq_cat = csp_total$freq_category)

haps <- unique(csp_total$haplotype_number)

# Function to get infection status based on mosquito age
get_infection <- function(x) {
  inf <- ifelse(x < 2, 0,
                ifelse(x < 3, rbinom(1, size = 1, p = 0.05),
                       ifelse(x < 4, rbinom(1, size = 1, p = 0.1),
                              ifelse(x < 5, rbinom(1, size = 1, p = 0.15), rbinom(1, size = 1, p = 0.25)))))
  return(inf)
}

# Function to get multiplicity of infection (MOI) for humans based on age
get_moi <- function(x) {
  moi <- ifelse(x < 3, rtpois(1, 4, a = 0, b = 5),
                ifelse(x < 4, rtpois(1, 6, a = 0, b = 10),
                       ifelse(x < 5, rtpois(1, 6, a = 0, b = 15), rtpois(1, 6, a = 0, b = 17))))
  return(moi)
}

# Function to reconstruct sample function to deal with only one haplotype case
resample <- function(x, size, replace = FALSE, prob = NULL) {
  if(length(x) == 1) {
    if (size > 1 && replace == FALSE) {
      stop("Cannot sample more elements than available without replacement.")
    }
    if (replace == TRUE) {
      return(rep(x, size))
    } else {
      return(x)
    }
  } else {
    return(sample(x, size = size, replace = replace, prob = prob))
  }
}

# Function to assign haplotypes to a human infection based on MOI
get_pers_infec <- function(x, haps, freq) {
  haps_index <- resample(haps, size = x, prob = freq)
  return(haps_index)
}

# Function to get old infections for humans (if any haplotype age >= 7)
get_old_p_infec2 <- function(x) {
  return((rowSums(x >= 7) > 0) * 1)
}

# Function to simulate mosquito death based on age
get_mos_death3 <- function(x) {
  thresholds <- c(3, 6, 9, 12, 15, Inf)
  probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.7)
  return(rbinom(length(x), size = 1, p = probs[sapply(x, function(i) sum(i > thresholds)) + 1]))
}

# Function to remove zero values and take minimum of non-zero entries
remove_0_values_take_min <- function(x) {
  if(length(x[x > 0]) > 0) {
    min_val <- min(x[x > 0])
  } else {
    min_val <- 0
  }
  return(min_val)
}

.format_dur <- function(secs) {
  h <- floor(secs/3600); m <- floor((secs %% 3600)/60); s <- round(secs %% 60, 1)
  sprintf("%02dh:%02dm:%04.1fs", h, m, s)
}

#############################
#### Simulation Function ####
#############################
run_biting_sim <- function(pr_symp_infec, pr_symp_non_infec, pr_clear, pr_off_feed, pr_on_feed_rainy, pr_on_feed_dry, 
                           pr_on_feed_moderate, pr_hum_to_mos, pr_mos_to_hum, num_loc, 
                           pr_num_biting, n_m, n_p, scenario_name, n_sim, 
                           proportion_suceptible, pr_suceptibility, pr_nonSuceptibility, n_days, 
                           proportion_mobile, pr_move, prob_matrix) {
  
  total_start_time <- Sys.time()
  time_per_sim <- numeric(n_sim)
  
  # storage for per-day elapsed time since total_start_time
  day_elapsed_sec <- matrix(NA_real_, nrow = n_sim, ncol = n_days)
  day_elapsed_hr  <- matrix(NA_real_, nrow = n_sim, ncol = n_days)
  
  mosquito_MOI_df <- matrix(NA_integer_, nrow = n_sim, ncol = 30 * n_days)
  eir_df <- matrix(NA, nrow = n_sim, ncol = sum(n_p))
  age_mos_df <- matrix(NA, nrow = n_sim, ncol = sum(n_m))
  age_human_haps_array <- array(NA, c(sum(n_p), length(haps) * n_days, n_sim))
  symptoms <- array(NA, c(sum(n_p), n_days, n_sim))
  location <- array(NA, c(sum(n_p), n_days, n_sim))
  initial_locs_matrix <- matrix(NA, nrow = sum(n_p), n_sim)
  
  ## Initialize global tracking logs (accumulated across simulations)
  # New: movement_log structure now includes TripID and LegID
  global_movement_log <- data.frame(
    PersonID    = integer(),
    Origin      = character(),
    Destination = character(),
    Day         = integer(),
    Simulation  = integer(),
    TripID      = integer(),   # trip episode number (per person)
    LegID       = integer(),   # 1,2,... within trip; 0 = return home
    stringsAsFactors = FALSE
  )
  
  global_human_to_mos_log <- data.frame(
    PersonID    = integer(),
    MosquitoID  = integer(),
    Location    = character(),
    HaplotypeID = character(),
    Day         = integer(),
    Simulation  = integer(),
    SourceStatus = character(),
    stringsAsFactors = FALSE
  )
  
  global_trans_chain_log <- data.frame(
    OriginAddress    = character(),
    TargetAddress    = character(),
    TargetPersonID   = integer(),
    MosquitoID       = integer(),
    HaplotypeID      = character(),
    Day              = integer(),
    Simulation       = integer(),
    SourceType       = character(),
    TransmissionType = character(),
    stringsAsFactors = FALSE
  )
  
  # Pre-read haplotype frequency data for mosquito infection
  hap_freq <- gt_df$freq
  hap_ids <- as.character(gt_df$hap)
  
  for(q in 1:n_sim) {
    sim_start_time <- Sys.time()  # Start time for current simulation
    
    ## Initialize simulation-specific tracking logs and mosquito infection origin vector
    movement_log <- data.frame(
      PersonID    = integer(),
      Origin      = character(),
      Destination = character(),
      Day         = integer(),
      Simulation  = integer(),
      TripID      = integer(),
      LegID       = integer(),
      stringsAsFactors = FALSE
    )
    
    human_to_mos_log <- data.frame(
      PersonID    = integer(),
      MosquitoID  = integer(),
      Location    = character(),
      HaplotypeID = character(),
      Day         = integer(),
      Simulation  = integer(),
      SourceStatus = character(),
      stringsAsFactors = FALSE
    )
    
    trans_chain_log <- data.frame(
      OriginAddress    = character(),
      TargetAddress    = character(),
      TargetPersonID   = integer(),
      MosquitoID       = integer(),
      HaplotypeID      = character(),
      Day              = integer(),
      Simulation       = integer(),
      SourceType       = character(),
      TransmissionType = character(),
      stringsAsFactors = FALSE
    )
    
    mosquito_origin <- rep(NA_character_, times = sum(n_m))
    mosquito_source_type <- rep(NA_character_, times = sum(n_m))
    
    # New: imported_hap_flag[i, h] = 1 if person i currently carries hap h that was acquired
    # from a different home location; 0 otherwise.
    imported_hap_flag <- matrix(0L, nrow = sum(n_p), ncol = length(haps))
    # has_imported_hap[i] = 1 if person i currently has at least one imported hap.
    has_imported_hap  <- rep(0L, sum(n_p))
    
    # Initialize people's locations
    init_locs_p <- rep(1:num_loc, n_p)
    initial_locs_matrix[, q] <- init_locs_p
    
    # Initialize mosquitoes' locations
    init_locs_m <- rep(1:num_loc, n_m)
    
    human_locs <- init_locs_p
    
    age_m <- rtpois(sum(n_m), 4, a = 0, b = 14)
    
    # Set starting infection status for people and mosquitoes
    inf_p <- rbinom(sum(n_p), size = 1, p = 0.3)
    idx <- which(inf_p == 1)
    moi_p <- integer(sum(n_p))
    if (length(idx)) moi_p[idx] <- rtpois(length(idx), 2, a = 0, b = 16)
    
    inf_m <- sapply(age_m, get_infection)
    idx <- which(inf_m == 1)
    moi_m <- integer(sum(n_m))
    if (length(idx)) moi_m[idx] <- sapply(age_m[idx], get_moi)
    
    # Distribute haplotypes to locations
    haps_per_loc_low <- c(rep(floor(length(which(gt_df$freq_cat == 1)) / length(n_p)), length(n_p) - 1),
                          length(which(gt_df$freq_cat == 1)) - sum(rep(floor(length(which(gt_df$freq_cat == 1)) / length(n_p)), length(n_p) - 1)))
    haps_per_loc_med <- c(rep(floor(length(which(gt_df$freq_cat == 2)) / length(n_p)), length(n_p) - 1),
                          length(which(gt_df$freq_cat == 2)) - sum(rep(floor(length(which(gt_df$freq_cat == 2)) / length(n_p)), length(n_p) - 1)))
    haps_per_loc_high <- c(rep(floor(length(which(gt_df$freq_cat == 3)) / length(n_p)), length(n_p) - 1),
                           length(which(gt_df$freq_cat == 3)) - sum(rep(floor(length(which(gt_df$freq_cat == 3)) / length(n_p)), length(n_p) - 1)))
    
    hap_loc_list <- vector(mode = 'list', length = length(n_p))
    gt_df_new <- gt_df
    
    for(i in 1:length(n_p)) {
      hap_loc_list[[i]] <- c(resample(gt_df_new$hap[which(gt_df_new$freq_cat == 1)],
                                      size = haps_per_loc_low[i], replace = F),
                             resample(gt_df_new$hap[which(gt_df_new$freq_cat == 2)],
                                      size = haps_per_loc_med[i], replace = F),
                             resample(gt_df_new$hap[which(gt_df_new$freq_cat == 3)],
                                      size = haps_per_loc_high[i], replace = F))
      gt_df_new <- gt_df_new[-which(gt_df_new$hap %in% unlist(hap_loc_list[[i]])), ]
    }
    
    # Initialize human infection haplotypes
    infec_p <- list()
    for(i in 1:num_loc) {
      mois <- moi_p[which(init_locs_p == i)]
      inf_p_loc <- list()
      for(j in 1:n_p[i]) {
        i_p <- list(get_pers_infec(mois[j], hap_loc_list[[i]], gt_df$freq[which(gt_df$hap %in% hap_loc_list[[i]])]))
        inf_p_loc <- c(inf_p_loc, i_p)
      }
      infec_p <- c(infec_p, inf_p_loc)
    }
    
    age_haps_p <- matrix(0, nrow = sum(n_p), ncol = length(haps))
    for(i in 1:sum(n_p)) {
      if(length(infec_p[[i]]) > 0) {
        age_haps_p[i, infec_p[[i]]] <- rpois(length(infec_p[[i]]), 20)
      }
    }
    
    # Initialize mosquito infection haplotypes
    infec_m <- list()
    for(i in 1:num_loc) {
      mois <- moi_m[which(init_locs_m == i)]
      inf_m_loc <- list()
      inf_m_loc <- assign_mos_haps(
        mois, 
        gt_df$freq[which(gt_df$hap %in% hap_loc_list[[i]])], 
        hap_loc_list[[i]]
      )
      infec_m <- c(infec_m, inf_m_loc)
    }
    
    age_haps_m <- matrix(0, nrow = sum(n_m), ncol = length(haps))
    for(i in 1:sum(n_m)) {
      if(length(infec_m[[i]]) > 0) {
        age_haps_m[i, infec_m[[i]]] <- resample(1:age_m[i], size = length(infec_m[[i]]), replace = T)
      }
    }
    
    old_pers_infec <- get_old_p_infec2(age_haps_p)
    
    symp_index <- rep(0, sum(n_p))
    symp_index[old_pers_infec == 1] <- rbinom(length(which(old_pers_infec == 1)), 1, pr_symp_infec)
    symp_index[old_pers_infec == 0] <- rbinom(length(which(old_pers_infec == 0)), 1, pr_symp_non_infec)
    symp_age <- ifelse(symp_index == 1 & old_pers_infec == 1, 1, 0)
    
    num_infec_bites <- rep(0, sum(n_p))
    min_age_haps <- apply(age_haps_m, 1, remove_0_values_take_min)
    
    prob_bit_last_3_days <- ifelse(age_m < 2, 0,
                                   ifelse(min_age_haps == 0 | (age_m - min_age_haps) >= 3, pr_on_feed_dry, pr_off_feed))
    bit_last_3_days <- rbinom(length(prob_bit_last_3_days), 1, prob_bit_last_3_days)
    
    suceptible_people <- resample(1:sum(n_p), size = proportion_suceptible * sum(n_p), replace = F)
    suceptible_prob <- ifelse(1:sum(n_p) %in% suceptible_people, pr_suceptibility, pr_nonSuceptibility)
    
    mobile_humans <- ifelse(1:sum(n_p) %in% resample(1:sum(n_p), size = proportion_mobile * sum(n_p), replace = FALSE), 1, 0)
    
    days_away <- rep(0, sum(n_p))
    last_day <- rep(0, sum(n_p))
    length_trip <- rep(0, sum(n_p))
    
    days_since_return <- rep(1000, sum(n_p))
    
    last_trip_loc <- rep(NA_character_, sum(n_p))
    
    # New: add trip bookkeeping per person
    trip_id  <- rep(0L, sum(n_p))  # counts how many trips this person has started
    trip_leg <- rep(0L, sum(n_p))  # counts legs within the current trip (1,2,3,..., 0 when at home)
    
    ##############################
    #### Start of Daily Loop #####
    ##############################
    for(r in 1:n_days) {
      # Increment mosquito age by one day
      age_m <- age_m + 1
      age_haps_m <- age_haps_m + 1
      age_haps_m[age_haps_m == 1] <- 0
      
      # Increment human haplotype ages
      age_haps_p <- age_haps_p + 1
      age_haps_p[age_haps_p == 1] <- 0
      
      days_since_return <- days_since_return + 1
      
      # Process mosquito deaths and update corresponding states
      death_m <- get_mos_death3(age_m)
      mosquito_origin[which(death_m == 1)] <- NA_character_
      mosquito_source_type[which(death_m == 1)] <- NA_character_
      
      age_m[which(death_m == 1)] <- 0
      inf_m[which(death_m == 1)] <- 0
      moi_m[which(death_m == 1)] <- 0
      infec_m[which(death_m == 1)] <- rep(list(integer(0)), sum(death_m == 1))
      age_haps_m[death_m == 1, ] <- 0
      bit_last_3_days[death_m == 1] <- 0
      
      # Process symptomatic infections (increment symptom age, clear if over threshold)
      symp_age[which(symp_age != 0)] <- symp_age[which(symp_age != 0)] + 1
      symp_age[symp_age > 14] <- 0
      
      
      # Clear human parasites on 2nd symptomatic day when treatment onsite
      treated_idx <- which(symp_age == 2)
      if (length(treated_idx) > 0) {
        infec_p[treated_idx] <- rep(list(integer(0)), length(treated_idx))     # wipe the haplotype list
        age_haps_p[treated_idx,] <- 0 # zero out all hap ages
        
        # New: When treatment clears infections, also clear imported flags
        imported_hap_flag[treated_idx, ] <- 0L
        has_imported_hap[treated_idx] <- 0L
        
        symp_age[treated_idx] <- 0    # clear symptoms
      }
      
      old_pers_infec <- get_old_p_infec2(age_haps_p)
      
      symp_index[old_pers_infec == 1] <- rbinom(length(old_pers_infec[old_pers_infec == 1]), 1, pr_symp_infec)
      symp_index[old_pers_infec == 0] <- rbinom(length(old_pers_infec[old_pers_infec == 0]), 1, pr_symp_non_infec)
      
      new_onset <- which(symp_age == 0 & symp_index == 1 & old_pers_infec == 1)
      if (length(new_onset) > 0) symp_age[new_onset] <- 1
      
      symptoms[which(symp_age == 1), r, q] <- 1
      
      # New: update clearing of human parasites 
      # Clear human parasites by age of parasites
      for(i in 1:sum(n_p)) {
        if(length(na.omit(infec_p[[i]])) > 0) {           
          for(j in na.omit(infec_p[[i]])) {
            clear <- ifelse(age_haps_p[i, j] >= 90, rbinom(1, 1, 0.95),
                            ifelse(age_haps_p[i, j] >= 65, rbinom(1, 1, 0.85),
                                   ifelse(age_haps_p[i, j] >= 30, rbinom(1, 1, pr_clear), 0)))
            if(clear == 1) {
              infec_p[[i]] <- infec_p[[i]][!is.na(infec_p[[i]]) & infec_p[[i]] != j]
              age_haps_p[i, j] <- 0
              imported_hap_flag[i, j] <- 0L   # clear imported flag when hap clears
            }
          }
        }
        # New: update whether person i still has any imported hap
        has_imported_hap[i] <- as.integer(any(imported_hap_flag[i, ] == 1L))
      }
      
      ################################
      ######## Personnel Movement ####
      ################################
      # New: To revise the movement logic to better track trips for all the travelers instead of only infected person ----------
      # 1) Process returns from trips whose last_day was flagged yesterday
      for(p in 1:sum(n_p)) {
        if(last_day[p] == 1) {
          # Only log if they are actually somewhere else (just to be safe)
          if(human_locs[p] != init_locs_p[p]) {
            movement_log <- rbind(
              movement_log,
              data.frame(
                PersonID = p,
                Origin = as.character(human_locs[p]),
                Destination = as.character(init_locs_p[p]),
                Day = r,
                Simulation = q,
                TripID = trip_id[p],
                LegID = 0L,
                stringsAsFactors = FALSE
              )
            )
            
            last_trip_loc[p] <- as.character(human_locs[p])
          }
          
          # New: Update state to "back home"
          human_locs[p]      <- init_locs_p[p]
          days_away[p]       <- 0
          length_trip[p]     <- 0
          days_since_return[p] <- 0
          trip_leg[p]        <- 0L    # reset leg counter when home
        }
      }
      
      # 2) People currently away from home (days_away > 0) can move between non-home locations
      #    but cannot go back to their home except via the explicit return movement above.
      away_idx <- which(days_away > 0 & length_trip > 0)
      if (length(away_idx) > 0) {
        for(i in away_idx) {
          current_loc <- human_locs[i]
          home_loc    <- init_locs_p[i]
          
          # Skip if somehow already at home (shouldn't happen due to return logic)
          if (current_loc == home_loc) next
          
          # Decide whether this away person moves today
          move_today <- rbinom(1, 1, pr_move[current_loc])
          if (move_today == 1) {
            # Candidate destinations: cannot stay in place, cannot go home
            cand_locs <- setdiff(1:num_loc, c(current_loc, home_loc))
            if (length(cand_locs) > 0) {
              dest_probs <- prob_matrix[current_loc, cand_locs]
              # Handle any NA from prob_matrix (e.g., if not fully specified)
              if (all(is.na(dest_probs))) {
                # Fall back to uniform if row is all NA
                dest_probs <- rep(1, length(cand_locs))
              }
              new_loc <- resample(cand_locs, size = 1, replace = FALSE, prob = dest_probs)
              
              # Update leg counter within the same trip
              trip_leg[i] <- trip_leg[i] + 1L
              
              # Log this movement
              movement_log <- rbind(
                movement_log,
                data.frame(
                  PersonID    = i,
                  Origin      = as.character(current_loc),
                  Destination = as.character(new_loc),
                  Day         = r,
                  Simulation  = q,
                  TripID      = trip_id[i],
                  LegID       = trip_leg[i],
                  stringsAsFactors = FALSE
                )
              )
              
              # Update location
              human_locs[i] <- new_loc
            }
          }
        }
      }
      
      # 3) People at home and not currently on a trip can start a new trip
      new_trip_today <- rep(0L, sum(n_p))
      for(loc in 1:num_loc) {
        home_candidates <- which(
          human_locs == loc &                # physically in this location
            init_locs_p == loc &             # this location is their true home
            mobile_humans == 1 &             # mobile individuals only
            days_away == 0 &                 # not currently away
            length_trip == 0                 # no active trip planned
        )
        
        if (length(home_candidates) > 0) {
          # Decide who starts a trip today from this location
          move_flags <- rbinom(length(home_candidates), 1, pr_move[loc])
          starters   <- home_candidates[move_flags == 1]
          
          if (length(starters) > 0) {
            for(i in starters) {
              home_loc <- init_locs_p[i]
              
              # Destination cannot be home; use prob_matrix row for routing
              cand_locs <- setdiff(1:num_loc, home_loc)
              dest_probs <- prob_matrix[home_loc, cand_locs]
              if (all(is.na(dest_probs))) {
                dest_probs <- rep(1, length(cand_locs))
              }
              new_loc <- resample(cand_locs, size = 1, replace = FALSE, prob = dest_probs)
              
              # New trip: increment trip_id and set first leg
              trip_id[i]  <- trip_id[i] + 1L
              trip_leg[i] <- 1L
              
              # Set travel duration for this trip
              length_trip[i] <- ceiling(rexp(1, 0.125))
              days_away[i]   <- 1
              
              # Log departure from home with leg 1
              movement_log <- rbind(
                movement_log,
                data.frame(
                  PersonID    = i,
                  Origin      = as.character(home_loc),
                  Destination = as.character(new_loc),
                  Day         = r,
                  Simulation  = q,
                  TripID      = trip_id[i],
                  LegID       = trip_leg[i],
                  stringsAsFactors = FALSE
                )
              )
              
              # Update location
              human_locs[i] <- new_loc
              new_trip_today[i] <- 1L
            }
          }
        }
      }
      
      # 4) Update days_away and recompute last_day based on updated trips
      if (any(days_away >= 1)) {
        days_away[days_away >= 1] <- days_away[days_away >= 1] + 1
      }
      # Ensure new trips for today are set to exactly 1 day away
      if (any(new_trip_today == 1L)) {
        days_away[new_trip_today == 1L] <- 1
      }
      
      # Recompute last_day for use tomorrow (when last_day==1, they will return home)
      last_day <- rep(0, sum(n_p))
      last_day[which(days_away > 0 & days_away == length_trip)] <- 1
      
      # Store current locations for this day
      location[, r, q] <- human_locs
      
      # End New edition --------------------------------------------------------- 
      
      # humans_moving <- rep(0, sum(n_p))
      # for(i in 1:length(n_p)) {
      #   idx <- which(human_locs == i & mobile_humans == 1 & (length_trip - days_away == 0))
      #   if(length(idx) > 0)
      #     humans_moving[idx] <- rbinom(length(idx), 1, pr_move[i])
      # }
      # 
      # length_trip[which(humans_moving == 1)] <- ceiling(rexp(length(which(humans_moving == 1)), 0.125))
      # days_away[which(days_away >= 1)] <- days_away[which(days_away >= 1)] + 1
      # days_away[which(humans_moving == 1)] <- 1
      # last_day <- rep(0, sum(n_p))
      # last_day[which(days_away > 0 & days_away == length_trip)] <- 1
      # 
      # # Process personnel movement: if an individual moves, update their location and record the event if they are infected.
      # for(i in 1:sum(n_p)) {
      #   if(humans_moving[i] == 1) {
      #     current_loc <- human_locs[i]
      #     # Randomly select a new location excluding current one
      #     new_loc <- resample((1:num_loc)[-current_loc], size = 1, prob = prob_matrix[current_loc, -current_loc])
      #     # If the individual is infected (here we check via symptom index or old infection flag), record the movement event
      #     if((r > 1) && (symp_index[i] == 1 || old_pers_infec[i] == 1)) {
      #       movement_log <- rbind(movement_log, data.frame(
      #         PersonID = i, Origin = as.character(current_loc),
      #         Destination = as.character(new_loc), Day = r, Simulation = q,
      #         stringsAsFactors = FALSE
      #       ))
      #     }
      #     human_locs[i] <- new_loc
      #   }
      # }
      # location[, r, q] <- human_locs
      
      ################################
      #### Mosquito Biting & Transmission
      ################################
      mos_bite <- matrix(0, sum(n_m), sum(n_p))
      which_mos_bite <- rep(0, sum(n_m))
      which_hum_bite <- rep(0, sum(n_p))
      person_bitten <- rep(0, sum(n_p))
      
      bit_last_3_days[bit_last_3_days != 0] <- bit_last_3_days[bit_last_3_days != 0] + 1
      bit_last_3_days[bit_last_3_days > 3] <- 0
      
      mos_biting_probs <- rep(pr_off_feed, length(age_m))
      mos_biting_probs[age_m < 2] <- 0 # No bite from Mosquito due to age < 2 days, see `get_infection` function
      ready <- (bit_last_3_days < 1) & (age_m >= 2)
      
      if (r %in% rainy_days)     mos_biting_probs[ready] <- pr_on_feed_rainy
      if (r %in% moderate_days)  mos_biting_probs[ready] <- pr_on_feed_moderate
      if (r %in% dry_days)       mos_biting_probs[ready] <- pr_on_feed_dry
      
      bites <- rbinom(length(mos_biting_probs), 1, mos_biting_probs)
      which_mos_bite <- (bites == 1)
      
      for(i in which(which_mos_bite)) {
        num_biting <- sample(c(1,2,3,4,5,6,7), size = 1, prob = pr_num_biting)
        people_at_loc <- which(human_locs == init_locs_m[i])
        
        if(length(people_at_loc) > 0) {
          actual_bites <- min(num_biting, length(people_at_loc))
          which_people_bite <- resample(people_at_loc, 
                                        size = actual_bites, 
                                        replace = FALSE, 
                                        prob = suceptible_prob[people_at_loc])
          
          person_bitten[which_people_bite] <- 1
          mos_bite[i, which_people_bite] <- 1
          which_hum_bite[which_people_bite] <- 1
        }
      }
      bit_last_3_days[which_mos_bite == 1] <- 1
      
      # Process interactions: For each bitten person, simulate transmission events
      for(i in which(person_bitten == 1)) {
        mos_index1 <- mos_bite[, i] == 1
        mos_index <- which(mos_index1)
        inf_bites <- rep(0, length(mos_index))
        
        for(j in 1:length(mos_index)) {
          old_haps_p_i <- which(age_haps_p[i, ] >= 14)
          
          if(length(old_haps_p_i) > 0 && symp_age[i] < 1) {
            transfer_haps <- rep(NA, length(old_haps_p_i))
            
            # New: Pre-compute travel / return / imported status for this person on this day
            is_traveler <- (human_locs[i] != init_locs_p[i])
            is_recent_returnee <- (human_locs[i] == init_locs_p[i] && days_since_return[i] <= 14)
            has_imported <- (has_imported_hap[i] == 1L)
            
            for(k in 1:length(old_haps_p_i)) {
              prob <- pr_hum_to_mos
              transfer <- rbinom(1, 1, prob)
              if(transfer == 1) {
                transfer_haps[k] <- old_haps_p_i[k]
                
                source_status_label <- "Local"
                parasite_origin_addr <- as.character(init_locs_p[i])
                
                # New: update the detail categories for travelers / returnees
                if (is_traveler) {
                  # away from home
                  if (days_away[i] <= 14) {
                    # Importation days while away (up to 14 days)
                    source_status_label <- "Traveler_Recent" # Count as Importation (Source)
                  } else {
                    # After 14 days away, treat as long-term traveler (no longer counted as importation)
                    source_status_label <- "Traveler_LongTerm" # Need to discuss
                  }
                } else if (is_recent_returnee) {
                  # At home and recently returned (<= 14 days since return)
                  if (has_imported) {
                    # Has at least one imported hap
                    source_status_label <- "Returnee_Imported"
                    if (!is.na(last_trip_loc[i])) {
                      parasite_origin_addr <- last_trip_loc[i]
                    } else {
                      parasite_origin_addr <- as.character(init_locs_p[i])
                    }
                  } else {
                    # No imported hap, but still within 14-day post-travel window
                    source_status_label <- "Returnee_Recent"
                    parasite_origin_addr <- as.character(init_locs_p[i])
                  }
                } else if (has_imported && human_locs[i] == init_locs_p[i] && days_since_return[i] < 1000) {
                  # Returned long ago but still carries imported hap – traveler with past trip
                  source_status_label <- "Returnee_Imported"
                  if (!is.na(last_trip_loc[i])) {
                    parasite_origin_addr <- last_trip_loc[i]
                  } else {
                    parasite_origin_addr <- as.character(init_locs_p[i])
                  }
                }
                
                mosquito_origin[mos_index[j]] <- parasite_origin_addr
                mosquito_source_type[mos_index[j]] <- source_status_label
                
                human_to_mos_log <- rbind(human_to_mos_log, data.frame(
                  PersonID = i, 
                  MosquitoID = mos_index[j], 
                  Location = as.character(human_locs[i]),
                  HaplotypeID = as.character(transfer_haps[k]), 
                  Day = r, 
                  Simulation = q,
                  SourceStatus = source_status_label, 
                  stringsAsFactors = FALSE
                ))
              }
            }
            new_haps <- setdiff(na.omit(transfer_haps), infec_m[[mos_index[j]]])
            infec_m[[mos_index[j]]] <- c(infec_m[[mos_index[j]]], new_haps)
            age_haps_m[mos_index[j], new_haps] <- 1
          }
          
          
          # Mosquito-to-human transmission: if the mosquito has old haplotypes (>=9 days) and the human is susceptible,
          # then transfer haplotype(s) from mosquito to human.
          old_haps_m_j <- which(age_haps_m[mos_index[j], ] >= 9)
          
          if(length(old_haps_m_j) > 0 && symp_age[i] == 0) {
            transfer_haps_m <- rep(NA, length(old_haps_m_j))
            for(l in 1:length(old_haps_m_j)) {
              prob_m <- pr_mos_to_hum
              transfer_m <- rbinom(1, 1, prob_m)
              if(transfer_m == 1) {
                transfer_haps_m[l] <- old_haps_m_j[l]
              }
            }
            new_haps_m <- setdiff(na.omit(transfer_haps_m), infec_p[[i]])
            infec_p[[i]] <- c(infec_p[[i]], new_haps_m)
            age_haps_p[i, new_haps_m] <- 1
            
            # New: If this infections originates from a different home location than the person's,
            # mark hap(s) as imported for this person.
            if(length(new_haps_m) > 0 && !is.na(mosquito_origin[mos_index[j]])) {
              origin_loc_for_person <- as.character(mosquito_origin[mos_index[j]])
              home_loc_for_person   <- as.character(init_locs_p[i])
              if(origin_loc_for_person != home_loc_for_person) {
                imported_hap_flag[i, new_haps_m] <- 1L
              }
              has_imported_hap[i] <- as.integer(any(imported_hap_flag[i, ] == 1L))
            }
            
            if(any(!is.na(transfer_haps_m))) {
              inf_bites[j] <- 1
              # If the mosquito was infected from a mobile individual (mosquito_origin not NA), record the full transmission chain
              if(!is.na(mosquito_origin[mos_index[j]])) {
                if(length(new_haps_m) > 0) {
                  origin_loc <- as.character(mosquito_origin[mos_index[j]])
                  target_loc <- as.character(human_locs[i])
                  source_type <- as.character(mosquito_source_type[mos_index[j]])
                  
                  is_imported_case <- FALSE
                  if(source_type %in% c("Traveler_Recent", "Returnee_Imported", "Returnee_Recent")) {
                    is_imported_case <- TRUE
                  }
                  
                  
                  trans_chain_log <- rbind(trans_chain_log, data.frame(
                    OriginAddress = origin_loc,
                    TargetAddress = target_loc,
                    TargetPersonID = i, 
                    MosquitoID = mos_index[j],
                    HaplotypeID = as.character(new_haps_m[1]), 
                    Day = r, 
                    Simulation = q,
                    SourceType = source_type,
                    TransmissionType = ifelse(is_imported_case, "Imported", "Local"),
                    stringsAsFactors = FALSE
                  ))
                }
              }
            }
          }
        }
        if(r >= 365) {
          num_infec_bites[i] <- num_infec_bites[i] + sum(inf_bites)
        }
      }  # End for each bitten person
      
      # If today is a mosquito sampling day, process sampling (keep original logic)
      if(r %in% mos_sample_days) {
        sample_index <- resample(1:sum(n_m), size = 30)
        mos_moi <- rep(NA, 30)
        for(t in 1:30) {
          mos_moi[t] <- length(na.omit(unlist(infec_m[[sample_index[t]]])))
        }
        mosquito_MOI_df[q, (1 + (30 * (r - 1))):(30 * r)] <- mos_moi
        mosquito_origin[sample_index] <- NA_character_
        mosquito_source_type[sample_index] <- NA_character_
        age_m[sample_index] <- 0
        inf_m[sample_index] <- 0
        moi_m[sample_index] <- 0
        infec_m[sample_index] <- rep(list(integer(0)), length(sample_index))
        bit_last_3_days[sample_index] <- 0
        age_haps_m[sample_index, ] <- 0
      }
      
      # Store daily human haplotype ages
      age_human_haps_array[, ((1 + (length(haps) * (r - 1))):(length(haps) * r)), q] <- age_haps_p
      
      # progress print + log (place at the END of the day loop)
      now <- Sys.time()
      elapsed_sec <- as.numeric(difftime(now, total_start_time, units = "secs"))
      day_elapsed_sec[q, r] <- elapsed_sec
      day_elapsed_hr[q, r]  <- elapsed_sec / 3600
      
      cat(sprintf(
        "[Progress] Sim %d/%d — Day %d/%d done | %s since start\n",
        q, n_sim, r, n_days, .format_dur(elapsed_sec)
      ))
      if (interactive()) flush.console()
      
    }  # End of daily loop
    
    eir_df[q, ] <- num_infec_bites
    age_mos_df[q, ] <- age_m
    sim_end_time <- Sys.time()
    time_per_sim[q] <- as.numeric(difftime(sim_end_time, sim_start_time, units = "hours"))
    print(paste("Time taken for simulation", q, ":", time_per_sim[q], "hours"))
    
    ## Add simulation-specific tracking logs with simulation number, and accumulate to global logs
    global_movement_log <- rbind(global_movement_log, movement_log)
    global_human_to_mos_log <- rbind(global_human_to_mos_log, human_to_mos_log)
    global_trans_chain_log <- rbind(global_trans_chain_log, trans_chain_log)
  }  # End for each simulation replicate
  
  total_end_time <- Sys.time()
  total_duration <- as.numeric(difftime(total_end_time, total_start_time, units = "hours"))
  print(paste("Total time taken for all simulations:", total_duration, "hours"))
  
  folder_path <- scenario_name
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
    cat("Folder created:", folder_path, "\n")
  } else {
    cat("Folder already exists:", folder_path, "\n")
  }
  
  saveRDS(mosquito_MOI_df, file = file.path(folder_path, paste0("mosquito_MOI_", scenario_name)))
  saveRDS(eir_df, file = file.path(folder_path, paste0("eir_", scenario_name)))
  saveRDS(age_mos_df, file = file.path(folder_path, paste0("mos_age_", scenario_name)))
  saveRDS(symptoms, file = file.path(folder_path, paste0("symptom_status_", scenario_name)))
  saveRDS(age_human_haps_array, file = file.path(folder_path, paste0("haplotype_age_", scenario_name)))
  saveRDS(location, file = file.path(folder_path, paste0("location_", scenario_name)))
  saveRDS(initial_locs_matrix, file = file.path(folder_path, paste0("initial_locs_", scenario_name)))
  
  # Write timing log to CSV file
  timing_log_days <- data.frame(
    Level          = "day",
    Simulation     = rep(1:n_sim, each = n_days),
    Day            = rep(1:n_days, times = n_sim),
    ElapsedSeconds = as.vector(t(day_elapsed_sec)),
    ElapsedHours   = as.vector(t(day_elapsed_hr))
  )
  
  timing_log_sims <- data.frame(
    Level          = "simulation",
    Simulation     = 1:n_sim,
    Day            = NA_integer_,
    ElapsedSeconds = round(time_per_sim * 3600, 3),
    ElapsedHours   = round(time_per_sim, 4)
  )
  
  timing_log_total <- data.frame(
    Level          = "total",
    Simulation     = "All",
    Day            = NA_integer_,
    ElapsedSeconds = round(as.numeric(difftime(total_end_time, total_start_time, units = "secs")), 3),
    ElapsedHours   = round(total_duration, 4)
  )
  
  timing_log <- rbind(timing_log_days, timing_log_sims, timing_log_total)
  write.csv(timing_log, file = file.path(folder_path, "timing_log.csv"), row.names = FALSE)
  
  # Write global tracking logs to CSV files
  write.csv(global_movement_log, file = file.path(folder_path, paste0("movement_log_", scenario_name, ".csv")), row.names = FALSE)
  write.csv(global_human_to_mos_log, file = file.path(folder_path, paste0("human_to_mos_log_", scenario_name, ".csv")), row.names = FALSE)
  write.csv(global_trans_chain_log, file = file.path(folder_path, paste0("transmission_chain_log_", scenario_name, ".csv")), row.names = FALSE)
  
  ## Create an origin-destination (OD) matrix based on the global_trans_chain_log
  # Rows: Origin; Columns: Destination; Diagonals set to 0.
  od_matrix <- matrix(0, nrow = num_loc, ncol = num_loc)
  rownames(od_matrix) <- as.character(1:num_loc)
  colnames(od_matrix) <- as.character(1:num_loc)
  
  if(nrow(global_trans_chain_log) > 0) {
    for(i in 1:nrow(global_trans_chain_log)) {
      origin <- as.character(global_trans_chain_log$OriginAddress[i])
      destination <- as.character(global_trans_chain_log$TargetAddress[i])
      if(origin != destination) {  # Ignore self-transmission
        if(origin %in% rownames(od_matrix) && destination %in% colnames(od_matrix)){
          od_matrix[origin, destination] <- od_matrix[origin, destination] + 1
        }
      }
    }
  }
  # Write OD matrix to CSV file
  write.csv(od_matrix, file = file.path(folder_path, paste0("OD_matrix_", scenario_name, ".csv")), row.names = TRUE)
  
  print(scenario_name)
  print(q)
}

generate_mobility_matrix <- function(n_cities = 1, cc = 0.2, 
                                     tc = 0.45, tt = 0.05, ct = 0.1) {
  n_towns = num_loc - n_cities
  prob_matrix <- matrix(NA, nrow = num_loc, ncol = num_loc)
  
  for (i in 1:num_loc) {
    for (j in 1:num_loc) {
      if (i != j) {
        # Determine probabilities based on city/town relationships
        if (i <= n_cities && j <= n_cities) {
          prob_matrix[i, j] <- cc # City-to-city
        } else if (i > n_cities && j <= n_cities) {
          prob_matrix[i, j] <- tc  # Town-to-city
        } else if (i > n_cities && j > n_cities) {
          prob_matrix[i, j] <- tt # Town-to-town
        } else {
          prob_matrix[i, j] <- ct # City-to-town
        }
      }
    }
  }
  
  return(prob_matrix)
}

prob_matrix <- generate_mobility_matrix()

run_biting_sim(
  pr_symp_infec=0.05, pr_symp_non_infec=0.05, pr_clear=0.85,
  pr_off_feed=0.007,
  pr_on_feed_rainy=0.174,
  pr_on_feed_dry=0.05*0.174/0.15,
  pr_on_feed_moderate=0.1*0.174/0.15*1.08,
  pr_hum_to_mos=0.60, pr_mos_to_hum=0.30,
  pr_num_biting=c(0.60, 0.33, 0.030, 0.007, 0.002, 0, 0),
  n_m=n_m, n_p=n_p, num_loc=num_loc,
  proportion_suceptible=0.20, pr_suceptibility=0.01, pr_nonSuceptibility=0.005,
  proportion_mobile=0.20,
  pr_move=c(rep(0.03,1), rep(0.06,2)),
  n_days=730, n_sim=1, scenario_name="tracking test",
  prob_matrix=prob_matrix
)