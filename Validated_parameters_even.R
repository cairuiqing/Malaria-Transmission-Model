library(dplyr)
library(extraDistr)
library(microbenchmark)
library(Rcpp)

sourceCpp("get_biting_status.cpp")
######################################################
########Fixed Parameters and re-used functions#######
#####################################################
n_p<- c(200, 200, 200, 200, 400, 400)
n_m<- c(30000, 30000, 30000, 30000, 60000, 60000)

num_loc <- length(n_p)

#day numbering, could feed every 3 days sample mosquitoes every 7 days, sample humans every month (except sick visits)
mos_sample_days<- seq(from=365,to=730,by=7)
hum_sample_days<- seq(from=365, to=730, by = 30)

dry_days<- c(1:59, 305:424, 670:730)
rainy_days<- c(60:151, 425:516)
moderate_days<- c(152:304, 517:669)



#all haplotypesß with frequencies
hap_table_csv<- read.csv("haplotype_frequencies.csv")

human_freq_data<- hap_table_csv%>%
  filter(gene=="CSP"& species=="Human")

human_haps<- human_freq_data$haplotype_number
human_freq<- human_freq_data$frequency


mos_freq_data<- hap_table_csv%>%
  filter(gene=="CSP"& species=="Mosquito")

mos_haps<- mos_freq_data$haplotype_number
mos_freq<- mos_freq_data$frequency

mos_gt_df<- data.frame(haps=mos_haps, 
                       freq=mos_freq)


csp_total<- hap_table_csv%>%
  filter(gene=="CSP")%>%
  group_by(haplotype_number)%>%
  summarize(overall_frequency= sum(frequency))%>%
  mutate(freq_category= ifelse(overall_frequency<=median(overall_frequency),1,
                               ifelse(overall_frequency>quantile(overall_frequency,0.75),3,2)))


gt_df<- data.frame(hap=csp_total$haplotype_number, 
                   freq=csp_total$overall_frequency, 
                   freq_cat= csp_total$freq_category)

haps<- unique(csp_total$haplotype_number)
#functions to get infection status and haplotype compostion for humans and mosquitoes

get_infection<- function(x){
  inf<-ifelse(x<2,0,ifelse(x<3,rbinom(1,size=1,p=0.05),ifelse(x<4,rbinom(1,size=1,p=0.1),ifelse(x<5,rbinom(1,size=1,0.15),rbinom(1,size=1,p=0.25)))))
  
  return(inf)
}
get_moi<- function(x){
  moi<- ifelse(x<3,rtpois(1,4,a=0,b=5),ifelse(x<4, rtpois(1,6,a=0,b=10), ifelse(x<5,rtpois(1,6,a=0,b=15),rtpois(1,6,a=0,b=17))))
  return(moi)
}



get_pers_infec<- function(x,haps,freq){
  haps_index<-sample(haps, size=x, prob=freq)
  return(haps_index)
}

#get_mos_infec<- function(x,haps,freq){
#haps_index<-sample(haps, size=x, prob=freq)
#return(haps_index)
#}

#haplotype exchange between human and mosquitoes if biting is happening

get_old_p_haps<- function(x){
  old_index<- which(x>=14)
  return(old_index)
}


get_old_p_infec<- function(x){
  old_p_infec<- rep(0,n_p)
  for(i in 1:n_p){
    if(any(x[i,]>=7)){
      old_p_infec[i]<-1
    }
  }
  return(old_p_infec)
}

get_old_p_infec2 <- function(x){
  return((rowSums(x >= 7) > 0) * 1)
}

#old_haps_p<-apply(age_haps_p,1,get_old_p_haps)

get_old_m_haps<- function(x){
  old_index<- which(x>=9)
  return(old_index)
}



#function to get the mosquito death as a function of age
get_mos_death<- function(x){
  death_m<- ifelse(x<3,rbinom(1,size=1,p=0.01),
                   ifelse(x<6,rbinom(1,size=1,p=0.05),
                          ifelse(x<9,rbinom(1,size=1,p=0.1),
                                 ifelse(x<12,rbinom(1,size = 1,p=0.25),
                                        ifelse(x<15,rbinom(1,size=1,p=0.5),rbinom(1,size=1,p=0.7))))))
  return(death_m)
}


get_mos_death2 <- function(x){
  thresholds <- c(3, 6, 9, 12, 15, Inf)
  probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.7)
  
  return(rbinom(1, size=1, p=probs[sum(x > thresholds) + 1]))
}


get_mos_death3 <- function(x){
  thresholds <- c(3, 6, 9, 12, 15, Inf)
  probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.7)
  return(rbinom(length(x), size=1, p=probs[sapply(x, function(i) sum(i > thresholds)) + 1]))
}


remove_0_values_take_min<- function(x){
  if(length(x[x>0])>0){
    min<- min(x[x>0])
  }
  else{
    min<-0
  }
}


#############################
####Simulation function#####
############################

run_biting_sim<- function(pr_symp_infec, pr_symp_non_infec, pr_clear, pr_off_feed, pr_on_feed_rainy, pr_on_feed_dry, 
                          pr_on_feed_moderate, pr_hum_to_mos, pr_mos_to_hum, num_loc, 
                          pr_num_biting, n_m,n_p, scenario_name, n_sim, 
                          proportion_suceptible, pr_suceptibility, pr_nonSuceptibility, n_days, 
                          proportion_mobile, pr_move, prob_matrix){
  
  # Setting timer for the simulation
  total_start_time <- Sys.time()
  time_per_sim <- numeric(n_sim)
  
  mosquito_MOI_df<- matrix(NA, nrow=n_sim, ncol=30*n_days)
  eir_df<- matrix(NA, nrow=n_sim,ncol=sum(n_p))
  age_mos_df<- matrix(NA, nrow=n_sim, ncol=sum(n_m))
  age_human_haps_array<- array(NA, c(sum(n_p),length(haps)*n_days, n_sim))
  symptoms<-array(NA, c(sum(n_p), n_days,n_sim) )
  location<- array(NA, c(sum(n_p),n_days,n_sim))
  initial_locs_matrix<- matrix(NA, nrow=sum(n_p), n_sim)
  
  for(q in 1:n_sim){
    sim_start_time <- Sys.time() # Starting time for current simulation
    
    #initial location:
    init_locs_p<- c()
    # for(i in 1:length(n_p)){
    #   locs<-rep(i, n_p[i])
    #   init_locs_p<- c(init_locs_p, locs)
    # }
    # modification: =============
    init_locs_p <- rep(1:num_loc, n_p)
    initial_locs_matrix[,q]<- init_locs_p
    
    # init_locs_m<- c()
    # for(i in 1:length(n_m)){
    #   locs<-rep(i, n_m[i])
    #   init_locs_m<- c(init_locs_m, locs)
    # }
    # modification: =============
    init_locs_m <- rep(1:num_loc, n_m)
    
    #mosquito ages
    
    human_locs<- init_locs_p
    
    age_m<- rtpois(sum(n_m), 4,a=0,b=14)
    
    
    
    #starting infection status and MOIs for people can change depending on location
    
    inf_p<-rbinom(sum(n_p), size=1, p=0.3)
    
    moi_p<- rep(NA, sum(n_p))
    # for(i in 1:length(inf_p)) moi_p[i]<-ifelse(inf_p[i]==1,rtpois(1,2,a=0,b=16),0)
    moi_p <- ifelse(inf_p == 1,rtpois(1, 2, a = 0 , b = 16), 0)
    
    
    #starting infection status for mosquito these are mosquito-age dependent older =more likely to be infected
    
    inf_m<- sapply(age_m, get_infection)
    
    
    #starting MOIs for mosquitoes, these are also age dependent,older=higher MOI
    
    # moi_m<- ifelse(inf_m==0,0,NA)
    # 
    # 
    # 
    # moi_m[which(is.na(moi_m))]<-sapply(age_m[which(is.na(moi_m))],get_moi)
    moi_m <- ifelse(inf_m == 0, 0, sapply(age_m[inf_m == 1], get_moi))
    
    #split up haps into locations
    
    haps_per_loc_low<- c(rep(floor(length(which(gt_df$freq_cat==1))/length(n_p)),length(n_p)-1),
                         length(which(gt_df$freq_cat==1))-sum(rep(floor(length(which(gt_df$freq_cat==1))/length(n_p)),length(n_p)-1)))
    
    haps_per_loc_med<- c(rep(floor(length(which(gt_df$freq_cat==2))/length(n_p)),length(n_p)-1),
                         length(which(gt_df$freq_cat==2))-sum(rep(floor(length(which(gt_df$freq_cat==2))/length(n_p)),length(n_p)-1)))
    
    haps_per_loc_high<- c(rep(floor(length(which(gt_df$freq_cat==3))/length(n_p)),length(n_p)-1),
                          length(which(gt_df$freq_cat==3))-sum(rep(floor(length(which(gt_df$freq_cat==3))/length(n_p)),length(n_p)-1)))
    
    
    
    
    hap_loc_list<- vector(mode= 'list', length=length(n_p))
    
    gt_df_new<- gt_df
    
    
    ## combine humand and mos haplotypes into 1 pool
    #assign those
    
    for(i in 1:length(n_p)){
      hap_loc_list[[i]]<- c(sample(gt_df_new$hap[which(gt_df_new$freq_cat==1)],
                                   size=haps_per_loc_low[i], replace=F),
                            sample(gt_df_new$hap[which(gt_df_new$freq_cat==2)],
                                   size=haps_per_loc_med[i], replace=F),
                            sample(gt_df_new$hap[which(gt_df_new$freq_cat==3)],
                                   size=haps_per_loc_high[i], replace=F))
      gt_df_new<- gt_df_new[-which(gt_df_new$hap%in%unlist(hap_loc_list[[i]])),]
      
    }
    
    #starting haplotypes for people and age of those haplotypes (time from infectious bite)
    infec_p<- list()
    
    
    for(i in 1:num_loc){
      
      mois<- moi_p[which(init_locs_p==i)]
      inf_p<- list()
      for(j in 1:n_p[i]){
        i_p<- list(get_pers_infec(mois[j],hap_loc_list[[i]],gt_df$freq[which(gt_df$hap%in%hap_loc_list[[i]])]))
        inf_p<- c(inf_p, i_p)
      }
      infec_p<-c(infec_p, inf_p)
    }
    
    
    age_haps_p<- matrix(0,nrow=sum(n_p), ncol=length(haps))
    
    for(i in 1:sum(n_p)){
      if(length(infec_p[[i]])>0){
        age_haps_p[i,infec_p[[i]]]<- rpois(length(infec_p[[i]]),20)
      }
      
    }
    #starting haplotypes for mosquitoes  and age of those haplotypes
    
    
    infec_m<- list()
    
    
    for(i in 1:num_loc){
      
      mois<- moi_m[which(init_locs_m==i)]
      inf_m_loc<- list()
      for(j in 1:n_m[i]){
        i_m<- assign_mos_haps(mois[j],gt_df$freq[which(gt_df$hap%in%hap_loc_list[[i]])],hap_loc_list[[i]])
        inf_m_loc<- c(inf_m_loc, i_m)
      }
      infec_m<-c(infec_m, inf_m_loc)
    }
    
    #for(i in 1:(n_m-n_m_a1)){
    #infec_m[[which(mos_loc=="A2")[i]]]<- get_mos_infec(moi_m[which(mos_loc=="A2")[i]],
    #mos_haps_a2,mos_freq_a2)
    #}
    
    age_haps_m<- matrix(0,nrow=sum(n_m), ncol=length(haps))
    
    for(i in 1:sum(n_m)){
      if(length(infec_m[[i]])>0){
        age_haps_m[i,infec_m[[i]]]<- sample(1:age_m[i],size=length(infec_m[[i]]), replace=T)
      }
    }
    
    
    old_pers_infec<- get_old_p_infec2(age_haps_p)
    
    symp_index<- rep(0, sum(n_p))
    
    symp_index[old_pers_infec==1]<-rbinom(length(which(old_pers_infec==1)),1,pr_symp_infec)
    symp_index[old_pers_infec==0]<-rbinom(length(which(old_pers_infec==0)),1,pr_symp_non_infec)
    
    
    # symp_age<-rep(0,sum(n_p))
    # 
    # 
    # symp_age[which(symp_index==1&old_pers_infec==1)]<- 1
    # modification ==============
    symp_age <- ifelse(symp_index == 1 & old_pers_infec == 1, 1, 0)
    
    num_infec_bites<- rep(0,sum(n_p))
    
    
    
    
    min_age_haps<- apply(age_haps_m,1,remove_0_values_take_min)
    
    
    prob_bit_last_3_days<- ifelse(age_m<2,0,ifelse(min_age_haps==0 | (age_m-min_age_haps)>=3, pr_on_feed_dry, pr_off_feed))
    
    # bit_last_3_days<- rep(0,sum(n_m))
    # 
    # 
    # bit_last_3_days[prob_bit_last_3_days==pr_on_feed_dry]<-rbinom(length(prob_bit_last_3_days[prob_bit_last_3_days==pr_on_feed_dry]),1,pr_on_feed_dry)
    # bit_last_3_days[prob_bit_last_3_days==pr_off_feed]<-rbinom(length(prob_bit_last_3_days[prob_bit_last_3_days==pr_off_feed]),1,pr_off_feed)
    # Modification: =========
    bit_last_3_days <- ifelse(prob_bit_last_3_days == pr_on_feed_dry, 
                              rbinom(length(prob_bit_last_3_days), 1, pr_on_feed_dry),
                              rbinom(length(prob_bit_last_3_days), 1, pr_off_feed))
    
    suceptible_people<- sample(1:sum(n_p), size=proportion_suceptible*sum(n_p), replace=F)
    # suceptible_prob<- rep(0, sum(n_p))
    # suceptible_prob[suceptible_people]<- pr_suceptibility
    # suceptible_prob[-suceptible_people]<- pr_nonSuceptibility
    # Modification: =======
    suceptible_prob <- ifelse(1:sum(n_p) %in% suceptible_people, pr_suceptibility, pr_nonSuceptibility)
    
    # mobile_humans<- rep(0,sum(n_p))
    # mobile_humans[sample(1:sum(n_p), size=proportion_mobile*sum(n_p), replace=F)]<-1
    # Modification: ========
    mobile_humans <- ifelse(1:sum(n_p) %in% sample(1:sum(n_p), 
                                                   size = proportion_mobile * sum(n_p), 
                                                   replace = FALSE), 1, 0)
    
    days_away<- rep(0,sum(n_p))
    last_day<- rep(0,sum(n_p))
    length_trip<-rep(0,sum(n_p))
    ##############################
    ######Start of timed sim#####
    #############################
    
    for(r in 1:n_days){
      
      #age 1 day 
      age_m<- age_m+1
      
      # age_haps_m[age_haps_m!=0]<- age_haps_m[age_haps_m!=0]+1
      # 
      # age_haps_p[age_haps_p!=0]<- age_haps_p[age_haps_p!=0]+1
      
      age_haps_m <- age_haps_m + 1
      age_haps_m[age_haps_m == 1] <- 0
      
      age_haps_p <- age_haps_p + 1
      age_haps_p[age_haps_p == 1] <- 0
      
      #is the mosquito dead
      # death_m<- sapply(age_m,get_mos_death)
      death_m <- get_mos_death3(age_m)
      #if so replace it with a new mosquito(newborn)
      
      age_m[which(death_m==1)]<- 0
      
      
      inf_m[which(death_m==1)]<- 0
      
      
      moi_m[which(death_m==1)]<- 0
      
      infec_m[which(death_m==1)]<- NA
      
      
      age_haps_m[death_m == 1, ] <- 0
      
      bit_last_3_days[death_m==1]<-0
      
      # for(i in 1:n_m){
      #   if(death_m[i]==1 ){
      #     
      #     
      #     age_haps_m[i,]<- 0
      #     
      #     
      #     
      #   }
      # }
      
      
      
      # symptomatic infections get  sampled (sick visit) and cleared immediately (drugs)
      
      symp_age[which(symp_age!=0)]<- symp_age[which(symp_age!=0)]+1
      # symp_age<- ifelse(symp_age>14,0,symp_age)
      symp_age[symp_age > 14] <- 0
      
      # ifelse(symp_age>14,0,symp_age)
      
      infec_p[which(symp_age==2)]<- NA
      age_haps_p[which(symp_age==2),]<-0
      
      
      
      
      #if(symp_index[i] == 0){
      #if(max(age_haps_p[i,])!=0){
      #if(length(which(age_haps_p[i,]>=7))>0& length(which(age_haps_p[i,]<=30))>0){
      old_infec_p<-get_old_p_infec2(age_haps_p)
      
      symp_index[old_infec_p==1]<- rbinom(length(old_infec_p[old_infec_p==1]),1,pr_symp_infec)
      symp_index[old_infec_p==0]<-rbinom(length(old_infec_p[old_infec_p==0]),1,pr_symp_non_infec)
      
      
      
      symp_age[which(symp_index==1&old_infec_p==1)]<-1
      #}
      #}
      #}
      
      
      
      # if(length(which(age_haps_p[i,]>=7))>0& length(which(age_haps_p[i,]<=30))>0&max(age_haps_p[i,])!=0
      #    &symp_index[i]==0){
      #   new_index<-rbinom(1,1,pr_symp)
      #   if(new_index==1){
      #     symp_index[i]<-1
      #     symp_age[i]<- 1
      #   }
      # }
      
      
      
      
      #symptomatics<- which(symp_age==1)
      #if(length(symptomatics)>0){
      #symp_moi<- rep(NA, length(symptomatics))
      #for(i in 1:length(symptomatics)){
      #symp_moi[i]<- length(na.omit(unlist(infec_p[[symptomatics[i]]])))
      #}
      #if(sum(symp_moi,na.rm=T)>0){
      #symptomatic_MOI_df[q,(1+(n_p*(r-1))):(((n_p*(r-1)))+length(symp_moi))]<- symp_moi
      #}
      #}
      
      symptoms[which(symp_age==1),r,q]<- 1
      
      
      
      
      
      #then clear human parasites by age of parasites
      
      
      
      for(i in 1:sum(n_p)){
        # if(length(na.omit(infec_p[[i]]))>0){      
        if(sum(!is.na(infec_p[[i]]))>0){           
          for(j in c(na.omit(infec_p[[i]]))){
            clear<- ifelse(age_haps_p[i,j]>=90,rbinom(1,1,0.95),ifelse(age_haps_p[i,j]>=65,rbinom(1,1,0.85),ifelse(age_haps_p[i,j]>=30,rbinom(1,1,pr_clear),0)))
          }
          
          if(clear==1){
            infec_p[[i]][which(infec_p[[i]]==j)]<- NA
            age_haps_p[i,j]<- 0
          }
        }
      }
      
      ################################
      #############Movement##########
      ################################
      
      
      for(p in 1:sum(n_p)){
        if(last_day[p]==1){
          human_locs[p]<-init_locs_p[p]
          days_away[p]<-0
          length_trip[p]<-0
        }
      }
      
      
      humans_moving<- rep(0, sum(n_p))
      for(i in 1:length(n_p)){
        humans_moving[human_locs==i& mobile_humans==1 & length_trip-days_away==0] <- rbinom(length(which(human_locs==i&mobile_humans==1&length_trip-days_away==0)), 1, pr_move[i])
      }
      
      length_trip[which(humans_moving==1)]<- ceiling(rexp(length(which(humans_moving==1)),0.125))
      
      days_away[which(days_away>=1)]<-days_away[which(days_away>=1)]+1
      days_away[which(humans_moving==1)]<-1
      last_day<- rep(0, sum(n_p))
      last_day[which(length_trip-days_away==1|length_trip==1)]<-1
      
      for(i in 1:sum(n_p)){
        for(j in 1:length(n_p)){
          if(humans_moving[i]==1&human_locs[i]==j){
            human_locs[i]<- sample(c(1:length(n_p))[-j], size=1, prob= prob_matrix[j, -j])
          }
        }
      }
      # humans_moving<- rep(0,n_p)
      # for(i in 1:n_p){
      #   if(humans_loc[i]=="A1"){
      #     humans_moving[i]<-rbinom(1,1,pr_move_a1)
      #   }
      #   else{
      #     humans_moving[i]<-rbinom(1,1,pr_move_a2)
      #   }
      # }
      
      
      
      location[,r,q]<- human_locs
      
      # for(i in 1:n_p){
      #   if(humans_loc[i]=="A1"&humans_moving[i]==1){
      #     humans_loc[i]<-"A2"
      #   }
      #   else if(humans_loc[i]=="A1"&humans_moving[i]==0){
      #     humans_loc[i]<-"A1"
      #   }
      #   else if(humans_loc[i]=="A2"&humans_moving[i]==1){
      #     humans_loc[i]<-"A1"
      #   }
      #   else{
      #     humans_loc[i]<-"A2"
      #   }
      # }
      
      #if its a feeding day
      
      
      
      
      
      
      mos_bite<- matrix(0, sum(n_m),sum(n_p))
      which_mos_bite <- rep(0, sum(n_m))
      which_hum_bite <- rep(0, sum(n_p))
      person_bitten <- rep(0, sum(n_p))
      
      bit_last_3_days[bit_last_3_days!=0]<-bit_last_3_days[bit_last_3_days!=0]+1
      bit_last_3_days[bit_last_3_days>3] = 0
      
      mos_biting_probs<-rep(0,sum(n_m))
      
      if(r%in% c(rainy_days,moderate_days)){
        mos_biting_probs[bit_last_3_days<1&age_m>=2]<-rbinom(length(mos_biting_probs[bit_last_3_days<1&age_m>=2]),1,pr_on_feed_rainy)
      }
      else{
        mos_biting_probs[bit_last_3_days<1&age_m>=2]<-rbinom(length(mos_biting_probs[bit_last_3_days<1&age_m>=2]),1,pr_on_feed_dry)
      }
      
      mos_biting_probs[bit_last_3_days>=1|age_m<2]<-rbinom(length(mos_biting_probs[bit_last_3_days>=1|age_m<2]),1,pr_off_feed)
      
      
      
      bites <- rbinom(sum(n_m), 1, mos_biting_probs)
      which_mos_bite <- bites == 1
      
      for(i in which(which_mos_bite)){
        num_biting<- sample(c(1,2,3,4,5,6,7),size=1,prob=pr_num_biting)
        
        ###CHANGED
        which_people_bite <- sample(which(human_locs==init_locs_m[i]),size=num_biting,replace=F, prob=suceptible_prob[which(human_locs==init_locs_m[i])])
        person_bitten[which_people_bite] <- 1
        mos_bite[i,which_people_bite] <- 1
        which_hum_bite[which_people_bite] <- 1
        
      }
      
      bit_last_3_days[which_mos_bite==1]<-1
      
      # for(i in 1:n_m){
      #   if(age_m[i]>=2){
      #     if(min(age_haps_m[i,])<=(age_m[i]-3)){
      #       bite <- rbinom(1,1,pr_off_feed)
      #     }else{
      #       bite <- rbinom(1,1,pr_on_feed)
      #     }
      # 
      #     # bite<- ifelse(min(age_haps_m[i,])<=(age_m[i]-3),rbinom(1,1,pr_off_feed),rbinom(1,1,pr_on_feed))
      # 
      #     if(bite==1){
      #       num_biting<- sample(c(1,2,3,4,5,6,7),size=1,prob=pr_num_biting)
      # 
      #       ###CHANGED
      #       which_people_bite <- sample(which(humans_loc==mos_loc[i]),size=num_biting,replace=F)
      #       person_bitten[which_people_bite] <- 1
      #       mos_bite[i,which_people_bite]<-1
      #       which_hum_bite[which_people_bite] <- 1
      #       # mos_bite[i,sample(1:200,size=num_biting,replace=F)]<-1
      #       which_mos_bite[i] <- 1
      #     }
      #   }
      # 
      # }
      
      #haplotype exchange if biting is happening
      
      # old_haps_p<-apply(age_haps_p,1,get_old_p_haps)
      # 
      # 
      # old_haps_m<-apply(age_haps_m,1,get_old_m_haps)
      # 
      # ##########################
      # # which_mos_bite<- ifelse(apply(mos_bite,1,sum)>0,1,0)
      # 
      # 
      # ##CHANGED
      # # which_mos_bite<- (rowSums(mos_bite)>0) * 1
      # 
      # old_hap_mosquitoes<- rep(0,n_m)
      # 
      # for(i in 1:n_m){
      #   if(length(old_haps_m)>0){
      #     if(length(old_haps_m[[i]])>0){
      #       old_hap_mosquitoes[i]<-1
      #     }
      #   }
      # }
      # 
      # hap_transfer_mos<- rep(0,n_m)
      # 
      # for(i in 1:n_m){
      #   if(old_hap_mosquitoes[i]==1&which_mos_bite[i]==1){
      #     hap_transfer_mos[i]<-1
      #   }
      # }
      
      
      ####################
      
      # # which_hum_bite<- ifelse(apply(mos_bite,2,sum)>0,1,0)
      # # which_hum_bite<- (colSums(mos_bite)>0)*1
      # 
      # old_hap_hum<- rep(0,n_p)
      # 
      # 
      # for(i in 1:n_p){
      #   if(length(old_haps_p)>0){
      #     if(length(old_haps_p[[i]])>0){
      #       old_hap_hum[i]<-1
      #     }
      #   }
      # }
      # 
      # hap_transfer_hum<- rep(0,n_p)
      # 
      # for(i in 1:n_p){
      #   if(old_hap_hum[i]==1&which_hum_bite[i]==1){
      #     hap_transfer_hum[i]<-1
      #   }
      # }
      # 
      ##########################
      for(i in which(person_bitten == 1)){
        mos_index1<- mos_bite[,i]==1
        mos_index<- which(mos_index1)
        inf_bites<- rep(0,length(mos_index))
        for(j in 1:length(mos_index)){
          
          
          old_haps_p_i <- which(age_haps_p[i, ] >= 14)
          if(length(old_haps_p_i)>0&symp_age[i]<1){
            transfer_haps<- rep(NA,length(old_haps_p_i))
            for(k in 1:length(old_haps_p_i)){
              prob<- pr_hum_to_mos
              transfer<- rbinom(1,1,prob)
              if(transfer==1){
                transfer_haps[k]<- old_haps_p_i[k]
              }
            }
            new_haps<- setdiff(na.omit(transfer_haps),infec_m[[mos_index[j]]])
            infec_m[[mos_index[j]]]<-c(infec_m[[mos_index[j]]],new_haps)
            age_haps_m[mos_index[j],new_haps]<-1
          }
          
          
          
          old_haps_m_j <- which(age_haps_m[mos_index[j], ] >= 9)
          if(length(old_haps_m_j)>0&symp_age[i]==0){
            transfer_haps_m<- rep(NA,length(old_haps_m_j))
            for(l in 1:length(old_haps_m_j)){
              prob_m<- pr_mos_to_hum
              transfer_m<- rbinom(1,1,prob_m)
              if(transfer_m==1){
                transfer_haps_m[l]<- old_haps_m_j[l]
              }
            }
            new_haps_m<- setdiff(na.omit(transfer_haps_m), infec_p[[i]])
            infec_p[[i]]<- c(infec_p[[i]],new_haps_m)
            age_haps_p[i,new_haps_m]<- 1
            if(length(new_haps_m)>0){
              inf_bites[j]<-1
            }
          }
          # }
          
          
          
          
        }
        if(r>=365){
          num_infec_bites[i]<- num_infec_bites[i]+sum(inf_bites)
        }
        # }
        
      }
      
      
      # for(i in 1:n_p){
      #   mos_index<- which(mos_bite[,i]==1)
      #   if(length(mos_index)>0){
      #     inf_bites<- rep(0,length(mos_index))
      #     for(j in 1:length(mos_index)){
      #       if(length(old_haps_p)>0){
      #         if(length(old_haps_p[[i]])>0&symp_age[i]<2){
      #           transfer_haps<- rep(NA,length(old_haps_p[[i]]))
      #           for(k in 1:length(old_haps_p[[i]])){
      #             prob<- pr_hum_to_mos
      #             transfer<- rbinom(1,1,prob)
      #             if(transfer==1){
      #               transfer_haps[k]<- old_haps_p[[i]][k]
      #             }
      #           }
      #           new_haps<- setdiff(na.omit(transfer_haps),infec_m[[mos_index[j]]])
      #           infec_m[[mos_index[j]]]<-c(infec_m[[mos_index[j]]],new_haps)
      #           age_haps_m[mos_index[j],new_haps]<-1
      #         }
      #       }
      #       
      #       if(length(old_haps_m)>0){
      #         if(length(old_haps_m[[mos_index[j]]])>0&symp_index[i]<1){
      #           transfer_haps_m<- rep(NA,length(old_haps_m[[mos_index[j]]]))
      #           for(l in 1:length(old_haps_m[[mos_index[j]]])){
      #             prob_m<- pr_mos_to_hum
      #             transfer_m<- rbinom(1,1,prob_m)
      #             if(transfer_m==1){
      #               transfer_haps_m[l]<- old_haps_m[[mos_index[j]]][l]
      #             }
      #           }
      #           new_haps_m<- setdiff(na.omit(transfer_haps_m), infec_p[[i]])
      #           infec_p[[i]]<- c(infec_p[[i]],new_haps_m)
      #           age_haps_p[i,new_haps_m]<- 1
      #           inf_bites[j]<-1
      #         }
      #         
      #         
      #         
      #       }
      #     }
      #     if(r>=365){
      #       num_infec_bites[i]<- num_infec_bites[i]+sum(inf_bites)
      #     }
      #   }
      #   
      # }
      
      
      #if its only a mosquito sample day
      if(r%in% mos_sample_days){
        sample_index<- sample(1:sum(n_m),size=30)
        mos_moi<- rep(NA,30)
        for(t in 1:30){
          mos_moi[t]<- length(na.omit(unlist(infec_m[[sample_index[t]]])))
        }
        mosquito_MOI_df[q,(1+(30*(r-1))):(30*r)]<- mos_moi
        
        
        
        age_m[sample_index]<- 0
        #rtpois(length(sample_index), 4,a=1,b=14)
        
        inf_m[sample_index]<- 0
        #sapply(age_m[sample_index],get_infection)
        
        moi_m[c(sample_index)]<- 0
        #,which(inf_m==1)
        #sapply(age_m[c(sample_index,which(inf_m==1))],get_moi)
        infec_m[c(sample_index)]<- NA
        #,which(inf_m==1&moi_m>0)
        #sapply(moi_m[c(sample_index,which(inf_m==1&moi_m>0))],get_mos_infec)
        
        bit_last_3_days[sample_index]<-0
        for(i in 1:sum(n_m)){
          if(i%in%sample_index ){
            
            #& inf_m[i]==1 & moi_m[i]>0
            age_haps_m[i,]<- 0
            #c(infec_m[[i]])
            #sample(1:age_m[i],size=length(infec_m[[i]]),replace=T)
          }
        }
        
        
      }
      #if its only a human sample day
      
      age_human_haps_array[,((1+(length(haps)*(r-1))):(length(haps)*r)),q]<- age_haps_p
      
      
    }
    
    eir_df[q,]<- num_infec_bites
    age_mos_df[q,]<- age_m
    
    
    sim_end_time <- Sys.time() # Ending time for current simulation
    time_per_sim[q] <- as.numeric(difftime(sim_end_time, sim_start_time, units="hours"))
    print(paste("Time taken for simulation", q, ":", time_per_sim[q], "hours"))
  }
  
  total_end_time <- Sys.time()
  total_duration <- as.numeric(difftime(total_end_time, total_start_time, units="hours"))
  print(paste("Total time taken for all simulations:", total_duration, "hours"))
  
  
  folder_path <- scenario_name
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
    cat("Folder created:", folder_path, "\n")
  } else {
    cat("Folder already exists:", folder_path, "\n")
  }
  saveRDS(mosquito_MOI_df, file = file.path(folder_path, 
                                            paste0("mosquito_MOI_",scenario_name)))
  saveRDS(eir_df, file = file.path(folder_path,
                                   paste0("eir_", scenario_name)))
  saveRDS(age_mos_df, file = file.path(folder_path,
                                       paste0("mos_age_",scenario_name)))
  saveRDS(symptoms, file = file.path(folder_path,
                                     paste0("symptom_status_",scenario_name)))
  saveRDS(age_human_haps_array, file = file.path(folder_path,
                                                 paste0("haplotype_age_",scenario_name)))
  saveRDS(location, file = file.path(folder_path,
                                     paste0("location_",scenario_name)))
  saveRDS(initial_locs_matrix, file = file.path(folder_path,
                                                paste0("initial_locs_",scenario_name)))
  # Write timing logs to a file
  timing_log <- data.frame(
    Simulation = 1:n_sim,
    TimeTaken = time_per_sim
  )
  timing_log <- rbind(timing_log, data.frame(Simulation = "Total", TimeTaken = total_duration))
  write.csv(timing_log, file = file.path(folder_path, "timing_log.csv"), row.names = FALSE)
  print(scenario_name)
  print(q)
}



# Current parameters
prob_matrix <- matrix(0.45, nrow = num_loc, ncol = num_loc)
diag(prob_matrix) <- NA
run_biting_sim(pr_symp_infec = 0.05, pr_symp_non_infec = 0.05, pr_clear = 0.85, pr_off_feed = 0.0075, 
               pr_on_feed_rainy = 0.14, pr_on_feed_dry = 0.05*0.14/0.15, pr_on_feed_moderate = 0.1*0.14/0.15, 
               pr_hum_to_mos = 0.6, pr_mos_to_hum = 0.3, num_loc = num_loc, 
               pr_num_biting = c(0.6, 0.33, 0.033, 0.008, 0, 0, 0), n_m = n_m, proportion_suceptible = 0.185, 
               pr_suceptibility = 0.008, pr_nonSuceptibility = 0.005, n_p = n_p, proportion_mobile = 0.1, 
               pr_move = rep(0.03, num_loc), n_days = 730, scenario_name = "Validated_parameters_even", 
               n_sim = 5, prob_matrix = prob_matrix
)