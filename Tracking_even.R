library(dplyr)
library(extraDistr)
library(microbenchmark)
library(Rcpp)

sourceCpp("get_biting_status.cpp")

######################################################
########Fixed Parameters and re-used functions#######
#####################################################
n_p <- c(20, 20, 20)
n_m <- c(3000, 3000, 3000)

num_loc <- length(n_p)

# day numbering, sample days etc.
mos_sample_days <- seq(from = 365, to = 730, by = 7)
hum_sample_days <- seq(from = 365, to = 730, by = 30)

dry_days <- c(1:59, 305:424, 670:730)
rainy_days <- c(60:151, 425:516)
moderate_days <- c(152:304, 517:669)

# haplotype frequencies
hap_table_csv <- read.csv("haplotype_frequencies.csv")

human_freq_data <- hap_table_csv %>%
  filter(gene == "CSP" & species == "Human")

human_haps <- human_freq_data$haplotype_number
human_freq <- human_freq_data$frequency

mos_freq_data <- hap_table_csv %>%
  filter(gene == "CSP" & species == "Mosquito")

mos_haps <- mos_freq_data$haplotype_number
mos_freq <- mos_freq_data$frequency

mos_gt_df <- data.frame(haps = mos_haps, 
                        freq = mos_freq)

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
# functions to get infection status and haplotype composition for humans and mosquitoes

get_infection <- function(x) {
  inf <- ifelse(x < 2, 0,
                ifelse(x < 3, rbinom(1, size = 1, p = 0.05),
                       ifelse(x < 4, rbinom(1, size = 1, p = 0.1),
                              ifelse(x < 5, rbinom(1, size = 1, p = 0.15), rbinom(1, size = 1, p = 0.25)))))
  return(inf)
}

get_moi <- function(x) {
  moi <- ifelse(x < 3, rtpois(1, 4, a = 0, b = 5),
                ifelse(x < 4, rtpois(1, 6, a = 0, b = 10),
                       ifelse(x < 5, rtpois(1, 6, a = 0, b = 15), rtpois(1, 6, a = 0, b = 17))))
  return(moi)
}

get_pers_infec <- function(x, haps, freq) {
  haps_index <- sample(haps, size = x, prob = freq)
  return(haps_index)
}

get_old_p_infec2 <- function(x) {
  return((rowSums(x >= 7) > 0) * 1)
}

get_mos_death3 <- function(x) {
  thresholds <- c(3, 6, 9, 12, 15, Inf)
  probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.7)
  return(rbinom(length(x), size = 1, p = probs[sapply(x, function(i) sum(i > thresholds)) + 1]))
}

remove_0_values_take_min <- function(x) {
  if(length(x[x > 0]) > 0) {
    min <- min(x[x > 0])
  } else {
    min <- 0
  }
  return(min)
}

#############################
#### Simulation function ####
############################

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
  
  ## 初始化全局追踪日志（跨simulation累积）
  global_movement_log    <- data.frame(PersonID = integer(), Origin = character(), 
                                       Destination = character(), Day = integer(), Simulation = integer())
  global_human_to_mos_log <- data.frame(PersonID = integer(), MosquitoID = integer(),
                                        Location = character(), HaplotypeID = character(), Day = integer(), Simulation = integer())
  global_trans_chain_log <- data.frame(OriginAddress = character(), TargetAddress = character(),
                                       TargetPersonID = integer(), MosquitoID = integer(),
                                       HaplotypeID = character(), Day = integer(), Simulation = integer())
  
  # 预先读取单倍型数据，用于蚊子感染
  hap_freq <- gt_df$freq
  hap_ids <- as.character(gt_df$hap)
  
  for(q in 1:n_sim) {
    sim_start_time <- Sys.time() # Starting time for current simulation
    
    ## 每个simulation独立的追踪日志和蚊子感染来源向量（总蚊子数量 = sum(n_m)）
    movement_log    <- data.frame(PersonID = integer(), Origin = character(), 
                                  Destination = character(), Day = integer())
    human_to_mos_log <- data.frame(PersonID = integer(), MosquitoID = integer(),
                                   Location = character(), HaplotypeID = character(), Day = integer())
    trans_chain_log <- data.frame(OriginAddress = character(), TargetAddress = character(),
                                  TargetPersonID = integer(), MosquitoID = integer(),
                                  HaplotypeID = character(), Day = integer())
    mosquito_origin <- rep(NA_character_, times = sum(n_m))
    
    # initial locations for people
    init_locs_p <- rep(1:num_loc, n_p)
    initial_locs_matrix[, q] <- init_locs_p
    
    # initial locations for mosquitoes
    init_locs_m <- rep(1:num_loc, n_m)
    
    human_locs <- init_locs_p
    
    age_m <- rtpois(sum(n_m), 4, a = 0, b = 14)
    
    # starting infection status for people and mosquitoes
    inf_p <- rbinom(sum(n_p), size = 1, p = 0.3)
    moi_p <- ifelse(inf_p == 1, rtpois(1, 2, a = 0, b = 16), 0)
    
    inf_m <- sapply(age_m, get_infection)
    moi_m <- ifelse(inf_m == 0, 0, sapply(age_m[inf_m == 1], get_moi))
    
    # 分配haplotypes到各个地点
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
    
    # 初始化人群感染haplotypes
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
    
    # 初始化蚊子haplotypes
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
    
    # 为简化示例，假设mosquito_infectious为所有蚊子在叮咬时均具传染性（实际可根据年龄等条件调整）
    mosquito_infectious <- rep(TRUE, sum(n_m))
    # 同时预设一个用于存储蚊子携带的单倍型（初始设为NA）
    mosquito_hap_id <- rep(NA_character_, sum(n_m))
    
    ##############################
    #### Start of daily loop #####
    #############################
    for(r in 1:n_days) {
      # 更新蚊子和人haplotype年龄
      age_m <- age_m + 1
      age_haps_m <- age_haps_m + 1
      age_haps_m[age_haps_m == 1] <- 0
      age_haps_p <- age_haps_p + 1
      age_haps_p[age_haps_p == 1] <- 0
      
      # 蚊子死亡及更新
      death_m <- get_mos_death3(age_m)
      age_m[which(death_m == 1)] <- 0
      inf_m[which(death_m == 1)] <- 0
      moi_m[which(death_m == 1)] <- 0
      infec_m[which(death_m == 1)] <- NA
      age_haps_m[death_m == 1, ] <- 0
      bit_last_3_days[death_m == 1] <- 0
      
      # 症状处理
      symp_age[which(symp_age != 0)] <- symp_age[which(symp_age != 0)] + 1
      symp_age[symp_age > 14] <- 0
      symp_index[old_pers_infec == 1] <- rbinom(length(old_pers_infec[old_pers_infec == 1]), 1, pr_symp_infec)
      symp_index[old_pers_infec == 0] <- rbinom(length(old_pers_infec[old_pers_infec == 0]), 1, pr_symp_non_infec)
      symp_age[which(symp_index == 1 & old_pers_infec == 1)] <- 1
      
      # 清除人群体内的parasites（省略部分细节）
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
      ######## 人员移动阶段 #########
      ################################
      # 重置部分：处理上次移动结束返回原位
      for(p in 1:sum(n_p)) {
        if(last_day[p] == 1) {
          human_locs[p] <- init_locs_p[p]
          days_away[p] <- 0
          length_trip[p] <- 0
        }
      }
      
      humans_moving <- rep(0, sum(n_p))
      for(i in 1:length(n_p)) {
        # 对于同一地点、且具备移动资格的人
        idx <- which(human_locs == i & mobile_humans == 1 & (length_trip - days_away == 0))
        if(length(idx) > 0)
          humans_moving[idx] <- rbinom(length(idx), 1, pr_move[i])
      }
      
      length_trip[which(humans_moving == 1)] <- ceiling(rexp(length(which(humans_moving == 1)), 0.125))
      days_away[which(days_away >= 1)] <- days_away[which(days_away >= 1)] + 1
      days_away[which(humans_moving == 1)] <- 1
      last_day <- rep(0, sum(n_p))
      last_day[which(length_trip - days_away == 1 | length_trip == 1)] <- 1
      
      # 人员移动：如果某人发生移动，则更新其位置，并记录移动事件（仅记录感染者的移动）
      for(i in 1:sum(n_p)) {
        # 假设根据模型规则，每人每天最多移动一次
        # 如果该人发生移动，则更新其位置
        if(humans_moving[i] == 1) {
          current_loc <- human_locs[i]
          # 随机选择一个新地点，排除原地点
          new_loc <- sample((1:num_loc)[-current_loc], size = 1, prob = prob_matrix[current_loc, -current_loc])
          # 如果该人已感染，则记录移动事件
          if( (r > 1) && (symp_index[i] == 1 || old_pers_infec[i] == 1) ) {
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
      #### 蚊虫叮咬与传播阶段 #######
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
      
      # 蚊虫与人体的交互：遍历每个被叮咬的人
      for(i in which(person_bitten == 1)) {
        mos_index1 <- mos_bite[, i] == 1
        mos_index <- which(mos_index1)
        inf_bites <- rep(0, length(mos_index))
        for(j in 1:length(mos_index)) {
          # 人 -> 蚊 传播：如果该人体内有足够老的haplotypes（>=14天），并且其症状状态符合条件，则进行单倍型转移
          old_haps_p_i <- which(age_haps_p[i, ] >= 14)
          if(length(old_haps_p_i) > 0 && symp_age[i] < 1) {
            transfer_haps <- rep(NA, length(old_haps_p_i))
            for(k in 1:length(old_haps_p_i)) {
              prob <- pr_hum_to_mos
              transfer <- rbinom(1, 1, prob)
              if(transfer == 1) {
                transfer_haps[k] <- old_haps_p_i[k]
                # 如果该人在今天发生了移动，则记录人传蚊事件
                if(any(movement_log$PersonID == i & movement_log$Day == r)) {
                  origin_addr <- movement_log$Origin[movement_log$PersonID == i & movement_log$Day == r][1]
                  human_to_mos_log <- rbind(human_to_mos_log, data.frame(
                    PersonID = i,
                    MosquitoID = mos_index[j],
                    Location = as.character(human_locs[i]),  # 新地址
                    HaplotypeID = as.character(transfer_haps[k]),
                    Day = r
                  ))
                  # 若该蚊子尚未记录感染来源，则记录之
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
          # 蚊 -> 人 传播：如果该蚊子体内有足够老的haplotypes（>=9天），且目标人当前未表现出症状，则进行单倍型转移
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
              # 如果该蚊子原先已因移动感染者而获得单倍型（mosquito_origin非NA），则记录完整传播链事件
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
      } # end for each bitten person
      
      # 若为蚊子抽样日，处理样本采集（保持原逻辑）
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
      
      # 更新每天的人类haplotype年龄存储
      age_human_haps_array[, ((1 + (length(haps) * (r - 1))):(length(haps) * r)), q] <- age_haps_p
    }  # end of daily loop
    
    eir_df[q, ] <- num_infec_bites
    age_mos_df[q, ] <- age_m
    sim_end_time <- Sys.time()
    time_per_sim[q] <- as.numeric(difftime(sim_end_time, sim_start_time, units = "hours"))
    print(paste("Time taken for simulation", q, ":", time_per_sim[q], "hours"))
    
    ## 将当前simulation的追踪日志添加上simulation编号后，累积到全局日志中
    movement_log$Simulation <- q
    human_to_mos_log$Simulation <- q
    trans_chain_log$Simulation <- q
    global_movement_log <- rbind(global_movement_log, movement_log)
    global_human_to_mos_log <- rbind(global_human_to_mos_log, human_to_mos_log)
    global_trans_chain_log <- rbind(global_trans_chain_log, trans_chain_log)
    
  }  # end for each simulation replicate
  
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
  
  # 写入时间日志
  timing_log <- data.frame(
    Simulation = 1:n_sim,
    TimeTaken = time_per_sim
  )
  timing_log <- rbind(timing_log, data.frame(Simulation = "Total", TimeTaken = total_duration))
  write.csv(timing_log, file = file.path(folder_path, "timing_log.csv"), row.names = FALSE)
  
  # 将追踪日志写入CSV文件
  write.csv(global_movement_log, file = file.path(folder_path, paste0("movement_log_", scenario_name, ".csv")), row.names = FALSE)
  write.csv(global_human_to_mos_log, file = file.path(folder_path, paste0("human_to_mos_log_", scenario_name, ".csv")), row.names = FALSE)
  write.csv(global_trans_chain_log, file = file.path(folder_path, paste0("transmission_chain_log_", scenario_name, ".csv")), row.names = FALSE)
  
  print(scenario_name)
  print(q)
}

# 当前参数调用示例
prob_matrix <- matrix(0.45, nrow = num_loc, ncol = num_loc)
diag(prob_matrix) <- NA
run_biting_sim(pr_symp_infec = 0.05, pr_symp_non_infec = 0.05, pr_clear = 0.85, pr_off_feed = 0.01, 
               pr_on_feed_rainy = 0.135, pr_on_feed_dry = 0.05 * 0.135 / 0.15, pr_on_feed_moderate = 0.1 * 0.135 / 0.15, 
               pr_hum_to_mos = 0.6, pr_mos_to_hum = 0.3, num_loc = num_loc, 
               pr_num_biting = c(0.6, 0.34, 0.03, 0.003, 0, 0, 0), n_m = n_m, proportion_suceptible = 0.2, 
               pr_suceptibility = 0.01, pr_nonSuceptibility = 0.005, n_p = n_p, proportion_mobile = 0.1, 
               pr_move = rep(0.03, num_loc), n_days = 730, scenario_name = "even_tracking_test_2", 
               n_sim = 1, prob_matrix = prob_matrix)
