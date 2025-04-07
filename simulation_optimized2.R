library(dplyr)
library(extraDistr)
library(microbenchmark)
library(Rcpp)
sourceCpp("get_biting_status_optimized.cpp")

total_mem_log <- data.frame(
  Day = as.integer(c()),
  Memory_Used = as.numeric(c()),
  Max_Used = as.numeric(c()),
  stringsAsFactors = FALSE
)
######################################################
######## Fixed Parameters and re-used functions#######
#####################################################
# day numbering, could feed every 3 days sample mosquitoes every 7 days, sample humans every month (except sick visits)
mos_sample_days <- seq(from = 365, to = 730, by = 7)
hum_sample_days <- seq(from = 365, to = 730, by = 30)

dry_days <- c(1:59, 305:424, 670:730)
rainy_days <- c(60:151, 425:516)
moderate_days <- c(152:304, 517:669)



# all haplotypesß with frequencies
hap_table_csv <- read.csv("haplotype_frequencies.csv")

human_freq_data <- hap_table_csv %>%
  filter(gene == "CSP" & species == "Human")

human_haps <- human_freq_data$haplotype_number
human_freq <- human_freq_data$frequency


mos_freq_data <- hap_table_csv %>%
  filter(gene == "CSP" & species == "Mosquito")

mos_haps <- mos_freq_data$haplotype_number
mos_freq <- mos_freq_data$frequency

mos_gt_df <- data.frame(
  haps = mos_haps,
  freq = mos_freq
)


csp_total <- hap_table_csv %>%
  filter(gene == "CSP") %>%
  group_by(haplotype_number) %>%
  summarize(overall_frequency = sum(frequency)) %>%
  mutate(freq_category = ifelse(overall_frequency <= median(overall_frequency), 1,
    ifelse(overall_frequency > quantile(overall_frequency, 0.75), 3, 2)
  ))


gt_df <- data.frame(
  hap = csp_total$haplotype_number,
  freq = csp_total$overall_frequency,
  freq_cat = csp_total$freq_category
)

haps <- unique(csp_total$haplotype_number)
# functions to get infection status and haplotype compostion for humans and mosquitoes

get_infection <- function(x) {
  inf <- ifelse(x < 2, 0,
    ifelse(x < 3, rbinom(1, size = 1, p = 0.05),
      ifelse(x < 4, rbinom(1, size = 1, p = 0.1),
        ifelse(x < 5, rbinom(1, size = 1, 0.15),
          rbinom(1, size = 1, p = 0.25)
        )
      )
    )
  )
  return(inf)
}

get_infection_optimized <- function(x) {
  # Define probability breaks and values
  prob_table <- data.frame(
    upper = c(2, 3, 4, 5, Inf),
    prob = c(0, 0.05, 0.1, 0.15, 0.25)
  )

  # Vectorized probability assignment
  probs <- prob_table$prob[findInterval(x, prob_table$upper)]

  # Vectorized random generation
  rbinom(length(x), size = 1, prob = probs)
}

get_moi <- function(x) {
  moi <- ifelse(x < 3, rtpois(1, 4, a = 0, b = 5),
    ifelse(x < 4, rtpois(1, 6, a = 0, b = 10),
      ifelse(x < 5, rtpois(1, 6, a = 0, b = 15),
        rtpois(1, 6, a = 0, b = 17)
      )
    )
  )
  return(moi)
}

get_moi_optimized <- function(x) {
  # Vectorized parameter assignment
  lambdas <- rep(6, length(x))
  lambdas[x < 3] <- 4

  upper_bounds <- rep(17, length(x))
  upper_bounds[x < 4] <- 10
  upper_bounds[x < 3] <- 5
  upper_bounds[x < 5] <- 15 # Note: This must come after x < 4

  # Vectorized truncated Poisson
  rtpois(length(x), lambda = lambdas, a = 0, b = upper_bounds)
}

get_pers_infec <- function(x, haps, freq) {
  haps_index <- sample(haps, size = x, prob = freq)
  return(haps_index)
}


get_old_p_haps <- function(x) {
  old_index <- which(x >= 14)
  return(old_index)
}


get_old_p_infec <- function(x) {
  old_p_infec <- rep(0, n_p)
  for (i in 1:n_p) {
    if (any(x[i, ] >= 7)) {
      old_p_infec[i] <- 1
    }
  }
  return(old_p_infec)
}

get_old_p_infec2 <- function(x) {
  return((rowSums(x >= 7) > 0) * 1)
}

# old_haps_p<-apply(age_haps_p,1,get_old_p_haps)

get_old_m_haps <- function(x) {
  old_index <- which(x >= 9)
  return(old_index)
}



# function to get the mosquito death as a function of age
get_mos_death <- function(x) {
  death_m <- ifelse(x < 3, rbinom(1, size = 1, p = 0.01),
    ifelse(x < 6, rbinom(1, size = 1, p = 0.05),
      ifelse(x < 9, rbinom(1, size = 1, p = 0.1),
        ifelse(x < 12, rbinom(1, size = 1, p = 0.25),
          ifelse(x < 15, rbinom(1, size = 1, p = 0.5), rbinom(1, size = 1, p = 0.7))
        )
      )
    )
  )
  return(death_m)
}


get_mos_death2 <- function(x) {
  thresholds <- c(3, 6, 9, 12, 15, Inf)
  probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.7)

  return(rbinom(1, size = 1, p = probs[sum(x > thresholds) + 1]))
}


get_mos_death3 <- function(x) {
  thresholds <- c(3, 6, 9, 12, 15, Inf)
  probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.7)
  return(rbinom(length(x), size = 1, p = probs[sapply(x, function(i) sum(i > thresholds)) + 1]))
}


remove_0_values_take_min <- function(x) {
  if (length(x[x > 0]) > 0) {
    min <- min(x[x > 0])
  } else {
    min <- 0
  }
}

create_result_folder <- function(scenario_name) {
  folder_path <- scenario_name
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
    cat("Folder created:", folder_path, "\n")
  } else {
    cat("Folder already exists:", folder_path, "\n")
  }
}

create_simulation_data_structures <- function(sim_round, n_days, people_total, mosquito_total, haps) {
  simulation_data <- list()
  simulation_data$mosquito_MOI_df <- matrix(NA, nrow = 1, ncol = 30 * n_days)
  simulation_data$eir_df <- matrix(NA, nrow = 1, ncol = people_total)
  simulation_data$age_mos_df <- matrix(NA, nrow = 1, ncol = mosquito_total)
  simulation_data$age_human_haps_arr <- array(NA, dim = c(people_total, length(haps) * n_days, 1))
  simulation_data$symptoms <- array(NA, dim = c(people_total, n_days, 1))
  simulation_data$location <- array(NA, dim = c(people_total, n_days, 1))
  simulation_data$initial_locs_matrix <- matrix(NA, nrow = people_total, ncol = 1)
  return(simulation_data)
}

#' Distribute haplotypes evenly across locations
#'
#' @param gt_df Dataframe containing haplotype data
#' @param freq_category Frequency category to distribute (1=low, 2=medium, 3=high)
#' @param n_locations Number of locations to distribute across
#' @return Vector with counts per location
distribute_haplotypes <- function(gt_df, freq_category, n_locations) {
  total <- length(which(gt_df$freq_cat == freq_category))
  res <- c(
    rep(floor(total / n_locations), n_locations - 1),
    total - sum(rep(floor(total / n_locations), n_locations - 1))
  )
  return(res)
}

combine_human_and_mos_haplotypes <- function(gt_df_new, n_p, haps_per_loc_low, haps_per_loc_med, haps_per_loc_high) {
  hap_loc_list <- vector("list", length(n_p))

  hap_by_freq <- split(gt_df_new$hap, gt_df_new$freq_cat) 

  used <- logical(nrow(gt_df_new))

  for (i in seq_along(n_p)) {
    # Sample from each category
    sampled_low <- sample(
      setdiff(hap_by_freq[["1"]], gt_df$hap[used]),
      haps_per_loc_low[i],
      replace = FALSE
    )
    sampled_med <- sample(
      setdiff(hap_by_freq[["2"]], gt_df$hap[used]),
      haps_per_loc_med[i],
      replace = FALSE
    )
    sampled_high <- sample(
      setdiff(hap_by_freq[["3"]], gt_df$hap[used]),
      haps_per_loc_high[i],
      replace = FALSE
    )
    
    # Store results
    hap_loc_list[[i]] <- c(sampled_low, sampled_med, sampled_high)
    
    # Update used haplotypes (faster than modifying gt_df repeatedly)
    used[gt_df$hap %in% hap_loc_list[[i]]] <- TRUE
  }
  # Return results (optionally, subset gt_df at the end if needed)
  list(
    hap_loc_list = hap_loc_list,
    updated_gt_df = gt_df[!used, ]
  )
}



run_biting_sim <- function(pr_symp_infec, pr_symp_non_infec, pr_clear, pr_off_feed, pr_on_feed_rainy, pr_on_feed_dry,
                           pr_on_feed_moderate, pr_hum_to_mos, pr_mos_to_hum, num_loc,
                           pr_num_biting, n_m, n_p, scenario_name, n_sim,
                           proportion_suceptible, pr_suceptibility, pr_nonSuceptibility, n_days,
                           proportion_mobile, pr_move, prob_matrix) {
  create_result_folder(scenario_name)
  #   mosquito_MOI_df<- matrix(NA, nrow=n_sim, ncol=30*n_days)
  #   eir_df<- matrix(NA, nrow=n_sim,ncol=sum(n_p))
  #   age_mos_df<- matrix(NA, nrow=n_sim, ncol=sum(n_m))
  #   age_human_haps_array<- array(NA, c(sum(n_p),length(haps)*n_days, n_sim))
  #   symptoms<-array(NA, c(sum(n_p), n_days,n_sim) )
  #   location<- array(NA, c(sum(n_p),n_days,n_sim))
  #   initial_locs_matrix<- matrix(NA, nrow=sum(n_p), n_sim)
  people_total <- sum(n_p)
  mosquito_total <- sum(n_m)
  for (sim_round in 1:n_sim) {
    cat(sprintf("Simulation Round %d start.\n", sim_round))
    simulation_data <- create_simulation_data_structures(sim_round, n_days, people_total, mosquito_total, haps)
    init_locs_p <- rep(1:num_loc, n_p)
    init_locs_m <- rep(1:num_loc, n_m)
    simulation_data$initial_locs_matrix[, 1] <- init_locs_p
    human_locs <- init_locs_p
    age_m <- rtpois(mosquito_total, 4, a = 0, b = 14)

    # starting infection status
    inf_p <- rbinom(people_total, size = 1, p = 0.3)
    moi_p <- ifelse(inf_p == 1, rtpois(1, 2, a = 0, b = 16), 0)
    # starting infection status for mosquito these are mosquito-age dependent older =more likely to be infected
    inf_m <- get_infection_optimized(age_m)
    moi_m <- numeric(length(age_m))
    moi_m[inf_m == 1] <- get_moi_optimized(age_m[inf_m == 1])
    moi_m[inf_m == 0] <- 0
    # starting infection status for mosquito these are mosquito-age dependent older =more likely to be infected
    haps_per_loc_low <- distribute_haplotypes(gt_df, 1, length(n_p))
    haps_per_loc_med <- distribute_haplotypes(gt_df, 2, length(n_p))
    haps_per_loc_high <- distribute_haplotypes(gt_df, 3, length(n_p))
    gt_df_new<- gt_df
    combined_human_and_mos_results <- combine_human_and_mos_haplotypes(gt_df_new, n_p, haps_per_loc_low, haps_per_loc_med, haps_per_loc_high)
    hap_loc_list <- combined_human_and_mos_results$hap_loc_list
    gt_df_new <- combined_human_and_mos_results$updated_gt_df
  }
}


n_p <- c(400, 400, 400, 200, 200, 200, 200, 200, 200, 200)
n_m <- c(60000, 60000, 60000, 30000, 30000, 30000, 30000, 30000, 30000, 30000)
prob_matrix <- matrix(c(
  NA, 0.18927454, 0.2172186, 0.08067833, 0.05954918, 0.08407578, 0.08209966, 0.08445692, 0.10393687, 0.09871013,
  0.2043349, NA, 0.2310639, 0.07063723, 0.04596831, 0.07313788, 0.07271722, 0.0762319, 0.11342066, 0.112488,
  0.21535214, 0.20712752, NA, 0.082191, 0.04834044, 0.08459333, 0.07636434, 0.06562914, 0.09783988, 0.1225622,
  0.49747389, 0.37058961, 0.07852352, NA, 0.00520151, 0.00746747, 0.00926119, 0.00463082, 0.01273598, 0.01411601,
  0.53759056, 0.35854353, 0.08910926, 0.00606895, NA, 0.00257361, 0.00253124, 0.00090631, 0.00075796, 0.00291857,
  0.53708057, 0.30871531, 0.12214735, 0.00678717, 0.00542284, NA, 0.00730853, 0.00124424, 0.00107299, 0.01022,
  0.49646372, 0.30731489, 0.14811581, 0.00341724, 0.00582395, 0.00789625, NA, 0.01285844, 0.00722665, 0.01088206,
  0.51845289, 0.30501799, 0.1175132, 0.00612425, 0.00563752, 0.00885674, 0.01020168, NA, 0.01002848, 0.01816726,
  0.47241647, 0.35579904, 0.12343639, 0.0086103, 0.00518435, 0.00520621, 0.00714952, 0.00980665, NA, 0.01239007,
  0.5076943, 0.3310045, 0.1304243, 0.00781055, 0.0061867, 0.00572986, 0.00814565, 0.0088086, 0.01539554, NA
), nrow = 10)

num_loc <- length(n_p)
run_biting_sim(
  pr_symp_infec = 0.05, pr_symp_non_infec = 0.05, pr_clear = 0.85,
  pr_off_feed = 0.01, pr_on_feed_rainy = 0.15,
  pr_on_feed_dry = 0.05, pr_on_feed_moderate = 0.1,
  pr_hum_to_mos = 0.6, pr_mos_to_hum = 0.3,
  pr_num_biting = c(0.6, 0.35, 0.04, 0.01, 0, 0, 0),
  n_m = n_m, proportion_suceptible = 0.2, pr_suceptibility = 0.01,
  pr_nonSuceptibility = 0.005, n_p = n_p, proportion_mobile = 0.1,
  pr_move = c(rep(0.03, 3), rep(0.05, 7)), num_loc = num_loc,
  n_days = 730, scenario_name = "Prop01Prob0305City3Town7",
  n_sim = 5, prob_matrix = prob_matrix
)
