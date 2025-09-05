library(dplyr)
library(extraDistr)
library(microbenchmark)
library(Rcpp)

sourceCpp("get_biting_status.cpp")

######################################################
######## Fixed Parameters and Re-used Functions ####
######################################################
n_p <- c(400, 400, 200, 200, 200, 200)
n_m <- c(60000, 60000, 30000, 30000, 30000, 30000)

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

# Function to get infection status for humans based on age
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

# Function to assign haplotypes to a human infection based on MOI
get_pers_infec <- function(x, haps, freq) {
  haps_index <- sample(haps, size = x, prob = freq)
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
  
  mosquito_MOI_df <- matrix(NA, nrow = n_sim, ncol = 30 * n_days)
  eir_df <- matrix(NA, nrow = n_sim, ncol = sum(n_p))
  age_mos_df <- matrix(NA, nrow = n_sim, ncol = sum(n_m))
  age_human_haps_array <- array(NA, c(sum(n_p), length(haps) * n_days, n_sim))
  symptoms <- array(NA, c(sum(n_p), n_days, n_sim))
  location <- array(NA, c(sum(n_p), n_days, n_sim))
  initial_locs_matrix <- matrix(NA, nrow = sum(n_p), n_sim)
  
  ## Initialize global tracking logs (accumulated across simulations)
  global_movement_log    <- data.frame(PersonID = integer(), Origin = character(), 
                                       Destination = character(), Day = integer(), Simulation = integer())
  global_human_to_mos_log <- data.frame(PersonID = integer(), MosquitoID = integer(),
                                        Location = character(), HaplotypeID = character(), Day = integer(), Simulation = integer())
  global_trans_chain_log <- data.frame(OriginAddress = character(), TargetAddress = character(),
                                       TargetPersonID = integer(), MosquitoID = integer(),
                                       HaplotypeID = character(), Day = integer(), Simulation = integer())
  
  # Pre-read haplotype frequency data for mosquito infection
  hap_freq <- gt_df$freq
  hap_ids <- as.character(gt_df$hap)
  
  for(q in 1:n_sim) {
    sim_start_time <- Sys.time()  # Start time for current simulation
    
    ## Initialize simulation-specific tracking logs and mosquito infection origin vector
    movement_log    <- data.frame(PersonID = integer(), Origin = character(), 
                                  Destination = character(), Day = integer())
    human_to_mos_log <- data.frame(PersonID = integer(), MosquitoID = integer(),
                                   Location = character(), HaplotypeID = character(), Day = integer())
    trans_chain_log <- data.frame(OriginAddress = character(), TargetAddress = character(),
                                  TargetPersonID = integer(), MosquitoID = integer(),
                                  HaplotypeID = character(), Day = integer())
    mosquito_origin <- rep(NA_character_, times = sum(n_m))
    
    # Initialize people's locations
    init_locs_p <- rep(1:num_loc, n_p)
    initial_locs_matrix[, q] <- init_locs_p
    
    # Initialize mosquitoes' locations
    init_locs_m <- rep(1:num_loc, n_m)
    
    human_locs <- init_locs_p
    
    age_m <- rtpois(sum(n_m), 4, a = 0, b = 14)
    
    # Set starting infection status for people and mosquitoes
    inf_p <- rbinom(sum(n_p), size = 1, p = 0.3)
    moi_p <- ifelse(inf_p == 1, rtpois(1, 2, a = 0, b = 16), 0)
    
    inf_m <- sapply(age_m, get_infection)
    moi_m <- ifelse(inf_m == 0, 0, sapply(age_m[inf_m == 1], get_moi))
    
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
      hap_loc_list[[i]] <- c(sample(gt_df_new$hap[which(gt_df_new$freq_cat == 1)],
                                    size = haps_per_loc_low[i], replace = F),
                             sample(gt_df_new$hap[which(gt_df_new$freq_cat == 2)],
                                    size = haps_per_loc_med[i], replace = F),
                             sample(gt_df_new$hap[which(gt_df_new$freq_cat == 3)],
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
      for(j in 1:n_m[i]) {
        i_m <- assign_mos_haps(mois[j], gt_df$freq[which(gt_df$hap %in% hap_loc_list[[i]])], hap_loc_list[[i]])
        inf_m_loc <- c(inf_m_loc, i_m)
      }
      infec_m <- c(infec_m, inf_m_loc)
    }
    
    age_haps_m <- matrix(0, nrow = sum(n_m), ncol = length(haps))
    for(i in 1:sum(n_m)) {
      if(length(infec_m[[i]]) > 0) {
        age_haps_m[i, infec_m[[i]]] <- sample(1:age_m[i], size = length(infec_m[[i]]), replace = T)
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
    bit_last_3_days <- ifelse(prob_bit_last_3_days == pr_on_feed_dry, 
                              rbinom(length(prob_bit_last_3_days), 1, pr_on_feed_dry),
                              rbinom(length(prob_bit_last_3_days), 1, pr_off_feed))
    
    suceptible_people <- sample(1:sum(n_p), size = proportion_suceptible * sum(n_p), replace = F)
    suceptible_prob <- ifelse(1:sum(n_p) %in% suceptible_people, pr_suceptibility, pr_nonSuceptibility)
    
    mobile_humans <- ifelse(1:sum(n_p) %in% sample(1:sum(n_p), size = proportion_mobile * sum(n_p), replace = FALSE), 1, 0)
    
    days_away <- rep(0, sum(n_p))
    last_day <- rep(0, sum(n_p))
    length_trip <- rep(0, sum(n_p))
    
    # For simplicity, assume all mosquitoes are infectious at biting time
    mosquito_infectious <- rep(TRUE, sum(n_m))
    # Pre-allocate a vector to store the haplotype carried by each mosquito (initially NA)
    mosquito_hap_id <- rep(NA_character_, sum(n_m))
    
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
      
      # Process mosquito deaths and update corresponding states
      death_m <- get_mos_death3(age_m)
      age_m[which(death_m == 1)] <- 0
      inf_m[which(death_m == 1)] <- 0
      moi_m[which(death_m == 1)] <- 0
      infec_m[which(death_m == 1)] <- NA
      age_haps_m[death_m == 1, ] <- 0
      bit_last_3_days[death_m == 1] <- 0
      
      # Process symptomatic infections (increment symptom age, clear if over threshold)
      symp_age[which(symp_age != 0)] <- symp_age[which(symp_age != 0)] + 1
      symp_age[symp_age > 14] <- 0
      
      
      # Clear human parasites on 2nd symptomatic day when treatment onsite
      treated_idx <- which(symp_age == 2)
      if (length(treated_idx) > 0) {
        infec_p[treated_idx] <-NA     # wipe the haplotype list
        age_haps_p[treated_idx,] <- 0 # zero out all hap ages
      }
      
      old_pers_infec <- get_old_p_infec2(age_haps_p)
      
      symp_index[old_pers_infec == 1] <- rbinom(length(old_pers_infec[old_pers_infec == 1]), 1, pr_symp_infec)
      symp_index[old_pers_infec == 0] <- rbinom(length(old_pers_infec[old_pers_infec == 0]), 1, pr_symp_non_infec)
      symp_age[which(symp_index == 1 & old_pers_infec == 1)] <- 1
      
      symptoms[which(symp_age == 1), r, q] <- 1
      
      # Clear human parasites by age of parasites
      for(i in 1:sum(n_p)) {
        if(sum(!is.na(infec_p[[i]])) > 0) {           
          for(j in na.omit(infec_p[[i]])) {
            clear <- ifelse(age_haps_p[i, j] >= 90, rbinom(1, 1, 0.95),
                            ifelse(age_haps_p[i, j] >= 65, rbinom(1, 1, 0.85),
                                   ifelse(age_haps_p[i, j] >= 30, rbinom(1, 1, pr_clear), 0)))
            if(clear == 1) {
              infec_p[[i]][which(infec_p[[i]] == j)] <- NA
              age_haps_p[i, j] <- 0
            }
          }
        }
      }
      
      ################################
      ######## Personnel Movement ####
      ################################
      # Reset: for individuals whose movement ended in the previous day, return to initial location
      for(p in 1:sum(n_p)) {
        if(last_day[p] == 1) {
          human_locs[p] <- init_locs_p[p]
          days_away[p] <- 0
          length_trip[p] <- 0
        }
      }
      
      humans_moving <- rep(0, sum(n_p))
      for(i in 1:length(n_p)) {
        idx <- which(human_locs == i & mobile_humans == 1 & (length_trip - days_away == 0))
        if(length(idx) > 0)
          humans_moving[idx] <- rbinom(length(idx), 1, pr_move[i])
      }
      
      length_trip[which(humans_moving == 1)] <- ceiling(rexp(length(which(humans_moving == 1)), 0.125))
      days_away[which(days_away >= 1)] <- days_away[which(days_away >= 1)] + 1
      days_away[which(humans_moving == 1)] <- 1
      last_day <- rep(0, sum(n_p))
      last_day[which(length_trip - days_away == 1 | length_trip == 1)] <- 1
      
      # Process personnel movement: if an individual moves, update their location and record the event if they are infected.
      for(i in 1:sum(n_p)) {
        if(humans_moving[i] == 1) {
          current_loc <- human_locs[i]
          # Randomly select a new location excluding current one
          new_loc <- sample((1:num_loc)[-current_loc], size = 1, prob = prob_matrix[current_loc, -current_loc])
          # If the individual is infected (here we check via symptom index or old infection flag), record the movement event
          if((r > 1) && (symp_index[i] == 1 || old_pers_infec[i] == 1)) {
            movement_log <- rbind(movement_log, data.frame(
              PersonID = i,
              Origin = as.character(current_loc),
              Destination = as.character(new_loc),
              Day = r
            ))
          }
          human_locs[i] <- new_loc
        }
      }
      location[, r, q] <- human_locs
      
      ################################
      #### Mosquito Biting & Transmission
      ################################
      mos_bite <- matrix(0, sum(n_m), sum(n_p))
      which_mos_bite <- rep(0, sum(n_m))
      which_hum_bite <- rep(0, sum(n_p))
      person_bitten <- rep(0, sum(n_p))
      
      bit_last_3_days[bit_last_3_days != 0] <- bit_last_3_days[bit_last_3_days != 0] + 1
      bit_last_3_days[bit_last_3_days > 3] <- 0
      
      mos_biting_probs <- rep(0, sum(n_m))
      if(r %in% c(rainy_days, moderate_days)) {
        mos_biting_probs[bit_last_3_days < 1 & age_m >= 2] <- rbinom(length(mos_biting_probs[bit_last_3_days < 1 & age_m >= 2]), 1, pr_on_feed_rainy)
      } else {
        mos_biting_probs[bit_last_3_days < 1 & age_m >= 2] <- rbinom(length(mos_biting_probs[bit_last_3_days < 1 & age_m >= 2]), 1, pr_on_feed_dry)
      }
      mos_biting_probs[bit_last_3_days >= 1 | age_m < 2] <- rbinom(length(mos_biting_probs[bit_last_3_days >= 1 | age_m < 2]), 1, pr_off_feed)
      
      bites <- rbinom(sum(n_m), 1, mos_biting_probs)
      which_mos_bite <- (bites == 1)
      
      for(i in which(which_mos_bite)) {
        num_biting <- sample(c(1,2,3,4,5,6,7), size = 1, prob = pr_num_biting)
        which_people_bite <- sample(which(human_locs == init_locs_m[i]), size = num_biting, replace = F, 
                                    prob = suceptible_prob[which(human_locs == init_locs_m[i])])
        person_bitten[which_people_bite] <- 1
        mos_bite[i, which_people_bite] <- 1
        which_hum_bite[which_people_bite] <- 1
      }
      bit_last_3_days[which_mos_bite == 1] <- 1
      
      # Process interactions: For each bitten person, simulate transmission events
      for(i in which(person_bitten == 1)) {
        mos_index1 <- mos_bite[, i] == 1
        mos_index <- which(mos_index1)
        inf_bites <- rep(0, length(mos_index))
        for(j in 1:length(mos_index)) {
          # Human-to-mosquito transmission: if the human has old haplotypes (>=14 days) and meets conditions,
          # then transfer haplotype(s) from human to mosquito.
          old_haps_p_i <- which(age_haps_p[i, ] >= 14)
          if(length(old_haps_p_i) > 0 && symp_age[i] < 1) {
            transfer_haps <- rep(NA, length(old_haps_p_i))
            for(k in 1:length(old_haps_p_i)) {
              prob <- pr_hum_to_mos
              transfer <- rbinom(1, 1, prob)
              if(transfer == 1) {
                transfer_haps[k] <- old_haps_p_i[k]
                # If the human moved today, record the human-to-mosquito event
                if(any(movement_log$PersonID == i & movement_log$Day == r)) {
                  origin_addr <- movement_log$Origin[movement_log$PersonID == i & movement_log$Day == r][1]
                  human_to_mos_log <- rbind(human_to_mos_log, data.frame(
                    PersonID = i,
                    MosquitoID = mos_index[j],
                    Location = as.character(human_locs[i]),  # New location
                    HaplotypeID = as.character(transfer_haps[k]),
                    Day = r
                  ))
                  # If the mosquito has not been assigned an origin yet, assign the human's original address as its origin.
                  if(is.na(mosquito_origin[mos_index[j]])) {
                    mosquito_origin[mos_index[j]] <- origin_addr
                  }
                }
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
            if(length(new_haps_m) > 0) {
              inf_bites[j] <- 1
              # If the mosquito was infected from a mobile individual (mosquito_origin not NA), record the full transmission chain
              if(!is.na(mosquito_origin[mos_index[j]])) {
                trans_chain_log <- rbind(trans_chain_log, data.frame(
                  OriginAddress = as.character(mosquito_origin[mos_index[j]]),
                  TargetAddress = as.character(human_locs[i]),
                  TargetPersonID = i,
                  MosquitoID = mos_index[j],
                  HaplotypeID = as.character(new_haps_m[1]),
                  Day = r
                ))
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
        sample_index <- sample(1:sum(n_m), size = 30)
        mos_moi <- rep(NA, 30)
        for(t in 1:30) {
          mos_moi[t] <- length(na.omit(unlist(infec_m[[sample_index[t]]])))
        }
        mosquito_MOI_df[q, (1 + (30 * (r - 1))):(30 * r)] <- mos_moi
        
        age_m[sample_index] <- 0
        inf_m[sample_index] <- 0
        moi_m[sample_index] <- 0
        infec_m[sample_index] <- NA
        bit_last_3_days[sample_index] <- 0
        for(i in 1:sum(n_m)) {
          if(i %in% sample_index) {
            age_haps_m[i, ] <- 0
          }
        }
      }
      
      # Store daily human haplotype ages
      age_human_haps_array[, ((1 + (length(haps) * (r - 1))):(length(haps) * r)), q] <- age_haps_p
    }  # End of daily loop
    
    eir_df[q, ] <- num_infec_bites
    age_mos_df[q, ] <- age_m
    sim_end_time <- Sys.time()
    time_per_sim[q] <- as.numeric(difftime(sim_end_time, sim_start_time, units = "hours"))
    print(paste("Time taken for simulation", q, ":", time_per_sim[q], "hours"))
    
    ## Add simulation-specific tracking logs with simulation number, and accumulate to global logs
    movement_log$Simulation <- q
    human_to_mos_log$Simulation <- q
    trans_chain_log$Simulation <- q
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
  timing_log <- data.frame(
    Simulation = 1:n_sim,
    TimeTaken = time_per_sim
  )
  timing_log <- rbind(timing_log, data.frame(Simulation = "Total", TimeTaken = total_duration))
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
        od_matrix[origin, destination] <- od_matrix[origin, destination] + 1
      }
    }
  }
  # Write OD matrix to CSV file
  write.csv(od_matrix, file = file.path(folder_path, paste0("OD_matrix_", scenario_name, ".csv")), row.names = TRUE)
  
  print(scenario_name)
  print(q)
}

generate_mobility_matrix <- function(n_cities = 2) {
  n_towns = num_loc - n_cities
  prob_matrix <- matrix(NA, nrow = num_loc, ncol = num_loc)
  
  for (i in 1:num_loc) {
    for (j in 1:num_loc) {
      if (i != j) {
        # Determine probabilities based on city/town relationships
        if (i <= n_cities && j <= n_cities) {
          prob_matrix[i, j] <- 0.2 # City-to-city
        } else if (i > n_cities && j <= n_cities) {
          prob_matrix[i, j] <- 0.45  # Town-to-city
        } else if (i > n_cities && j > n_cities) {
          prob_matrix[i, j] <- 0.05 # Town-to-town
        } else {
          prob_matrix[i, j] <- 0.1 # City-to-town
        }
      }
    }
  }
  
  return(prob_matrix)
}

prob_matrix <- generate_mobility_matrix()

run_biting_sim(pr_symp_infec = 0.05, pr_symp_non_infec = 0.05, pr_clear = 0.85, 
               pr_off_feed = 0.01, pr_on_feed_rainy = 0.135, 
               pr_on_feed_dry = 0.05*0.135/0.15, pr_on_feed_moderate = 0.1*0.135/0.15, 
               pr_hum_to_mos = 0.6, pr_mos_to_hum = 0.3, 
               pr_num_biting = c(0.6, 0.34, 0.03, 0.003, 0, 0, 0), 
               n_m=n_m, proportion_suceptible = 0.2, pr_suceptibility=0.01,
               pr_nonSuceptibility=0.005, n_p= n_p, proportion_mobile = 0.1, 
               pr_move = c(rep(0.03, 2), rep(0.06, 4)), num_loc=num_loc, 
               n_days = 730, scenario_name = "Tracking_uneven_final", 
               n_sim=5, prob_matrix=prob_matrix)