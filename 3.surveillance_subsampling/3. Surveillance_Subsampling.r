#==Load libraries==
library(Biostrings)
library(dplyr)
library(stringr)

#==Taking simulation 100 as an example==
for (q in 100) {
  ############################
  #0. Read simulated data
  ############################
  seq <- readDNAStringSet(paste0("../2.seq_generation/sim",q,"_sequences.fasta"))
  epi_data <- read.delim(paste0("../2.seq_generation/sim",q,"_transmission_chain.txt")) %>% 
    mutate(samp.dates1 = as.Date("2009-01-01")+infectious.time) %>%
    filter(samp.dates1 <= as.Date("2009-06-11"))
  run <- readRDS(paste0("../1.epi_model/res_v",q,".rds")); run$size
  run <- as.data.frame(t(run$dynamics$inf_s)) %>% mutate(date = seq(as.Date("2009-01-01")+0, as.Date("2009-01-01")+160, 1))
  
  ############################
  #1. Homogeneous surveillance
  ############################
  tmp2 <- list()
  set.seed(9)
  for (j in 1:5) { #For each simulation, we generated five sub-sampling replicates but only one was used to ensure the replicates were performed at the simulation level
    tmp <- c()
    for (i in 1:max(epi_data$patch)) {
      tmp_data <- epi_data %>% filter(samp.dates1 >= as.Date("2009-04-01")) %>% filter(patch == i) #surveillance was initiated after 1 April
      tmp_samp1 <- sample_n(tmp_data, size = ceiling(nrow(tmp_data)*0.1)) #10% infections in each region would be randomly sampled to become cases
      tmp_samp2 <- sample_n(tmp_samp1, size = ceiling(nrow(tmp_samp1)*0.05)) #5% cases in each region would be randomly sequenced
      tmp_samp3 <- sample_n(tmp_samp2, size = ceiling(nrow(tmp_samp2)*0.5)) #50% of sequences would be randomly shared to the public repository
      tmp1 <- tmp_samp3$traits
      tmp <- c(tmp, tmp1)
    }
    
    epi_data_tmp <- epi_data %>% filter(traits %in% tmp)
    while (epi_data_tmp$regions[which.min(epi_data_tmp$samp.dates1)] != "Mexico") { #To ensure the first case was sampled from Mexico
      print(paste0("re-run", j))
      tmp <- c()
      for (i in 1:max(epi_data$patch)) {
        tmp_data <- epi_data %>% filter(samp.dates1 >= as.Date("2009-04-01")) %>% filter(patch == i)
        tmp_samp1 <- sample_n(tmp_data, size = ceiling(nrow(tmp_data)*0.1))
        tmp_samp2 <- sample_n(tmp_samp1, size = ceiling(nrow(tmp_samp1)*0.05))
        tmp_samp3 <- sample_n(tmp_samp2, size = ceiling(nrow(tmp_samp2)*0.5))
        tmp1 <- tmp_samp3$traits
        tmp <- c(tmp, tmp1)
      }
      epi_data_tmp <- epi_data %>% filter(traits %in% tmp)
    }
    
    tmp2[[j]] <- tmp 
  }

  for (j in 1:5) {
    writeXStringSet(seq[names(seq) %in% tmp2[[j]]], paste0("homo_sampling/seq_homo_sim",q,"_", j ,".fasta"), format="fasta")
    write.table(epi_data %>% filter(traits %in% tmp2[[j]]), paste0("homo_sampling/meta_homo_sim",q,"_", j ,".txt"), quote = FALSE, sep = "\t")
  }
  
  ##############################
  #2. Heterogeneous surveillance (varying genomic surveillance intensity in Mexico)
  ##############################
  rate <- c(0.0450, 0.0400, 0.0350, 0.0300, 0.0250, 0.0200, 0.0150, 0.0100, 0.0050, 0.0010) 
  tmp6 <- list()
  set.seed(99)
  for (z in c(1:length(rate))) {
    tmp5 <- list()
    for (i in 1:5) {
      if (z == 1) {
        data <- epi_data %>% filter(traits %in% tmp2[[i]])
      } else {
        data <- epi_data %>% filter(traits %in% tmp6[[(z-1)]][[i]])
      }
      
      data1 <- data %>% filter(regions != "Mexico")
      data2 <- data %>% filter(regions == "Mexico")
      
      if (z == 1) {
        data3 <- data2 %>% slice_sample(n = ceiling(nrow(data2)*(rate[z]/0.05)), replace = F)
      } else {
        data3 <- data2 %>% slice_sample(n = ceiling(nrow(data2)*(rate[z]/rate[(z-1)])), replace = F)
      }
      
      tmp5[[i]] <- c(data1$traits, data3$traits)
    }
    tmp6[[z]] <- tmp5
  }
  
  for (i in c(1:length(rate))) {
    for (j in 1:5) {
      writeXStringSet(seq[names(seq) %in% tmp6[[i]][[j]]], paste0("rate_sampling/seq_rate_",rate[i],"_sim",q,"_", j ,".fasta"), format="fasta")
      write.table(epi_data %>% filter(traits %in% tmp6[[i]][[j]]), paste0("rate_sampling/meta_rate_",rate[i],"_sim",q,"_", j ,".txt"), quote = FALSE, sep = "\t")
    }
  }
  
  ##############################
  #3. Heterogeneous surveillance (varying time window performing genomic surveillance in Mexico)
  ##############################
  window <- c(seq(as.Date("2009-04-05"), as.Date("2009-04-05")+63, 7))
  
  tmp8 <- list()
  set.seed(99)
  for (z in 1:length(window)) {
    tmp7 <- list()
    for (j in 1:5) {
      data <- epi_data %>% filter(traits %in% tmp2[[j]])
      data1 <- data %>% filter(regions != "Mexico")
      data2 <- data %>% filter(regions == "Mexico")
      data3 <- data2 %>% filter(samp.dates1 >= window[z])
      tmp7[[j]] <- c(data1$traits, data3$traits)
    }
    tmp8[[z]] <- tmp7
  }
  
  for (i in 1:length(window)) {
    for (j in 1:5) {
      writeXStringSet(seq[names(seq) %in% tmp8[[i]][[j]]], paste0("window_sampling/seq_window_",window[i],"_sim",q,"_", j ,".fasta"), format="fasta")
      write.table(epi_data %>% filter(traits %in% tmp8[[i]][[j]]), paste0("window_sampling/meta_window_",window[i],"_sim",q,"_", j ,".txt"), quote = FALSE, sep = "\t")
    }
  }

  ##############################
  #4. Heterogeneous surveillance (varying data sharing timeliness)
  ##############################
  set.seed(99)
  tmp12 <- c()
  for (j in 1:5) {
    print(j)
    data <- epi_data %>% filter(traits %in% tmp2[[j]])
    
    for (i in 1:nrow(data)) {
      if (data$regions[i] == "Mexico") {
        m = 63
      } else {
        m = 21
      }
      delay <- max(7, round(rgamma(1, shape = m/7, scale = 7), 0))
      data$sub.delay[i] <- delay
      data <- data %>% mutate(type = paste0("run",j))
    }
    
    tmp12 <- rbind(tmp12, data)
  }
  
  tmp12$sub.dates <- tmp12$samp.dates1 + tmp12$sub.delay

  dates <- c(seq(as.Date("2009-06-08"), as.Date("2009-06-08")+63, 7))

  tmp13 <- list()
  for (i in 1:length(dates)) {
    tmp14 <- list()
    for (j in 1:5) {
      data <- tmp12 %>% filter(type == paste0("run", j))
      data <- data %>% 
        filter(sub.dates <= dates[i])
      tmp14[[j]] <- data$traits
    }
    tmp13[[i]] <- tmp14
  }
  
  for (i in 1:length(dates)) {
    for (j in 1:5) {
      writeXStringSet(seq[names(seq) %in% tmp13[[i]][[j]]], paste0("timeliness_sampling/seq_timeliness_",dates[i],"_sim",q,"_", j ,".fasta"), format="fasta")
      write.table(epi_data %>% filter(traits %in% tmp13[[i]][[j]]), paste0("timeliness_sampling/meta_timeliness_",dates[i],"_sim",q,"_", j ,".txt"), quote = FALSE, sep = "\t")
    }
  }
  
  ##############################
  #=======5.combined scenarios=======
  ##############################
  rate1 <- c(0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02)
  window1 <- seq(as.Date("2009-04-05"), as.Date("2009-04-05")+42, 7)
  date1 <- c(as.Date("2009-07-10"), as.Date("2009-08-10"), as.Date("2009-09-10"))

  tmp17 <- c()
  for (i in 1:length(rate1)) {
    tmp16 <- list()
    for (j in 1:5) {
      if (i == 1) {
        data <- tmp12 %>% filter(type == paste0("run",j)) %>% filter(traits %in% tmp2[[j]])
      } else {
        data <- tmp12 %>% filter(type == paste0("run",j)) %>% filter(traits %in% tmp6[[i-1]][[j]])
      }
      
      tmp15 <- list()
      for (z in 1:length(window1)) {
        data1 <- data %>% filter(regions != "Mexico")
        data2 <- data %>% filter(regions == "Mexico")
        data3 <- data2 %>% filter(samp.dates1 >= window1[z])
        data4 <- rbind(data1, data3)
        
        tmp15[[z]] <- data4$traits
      }
      tmp16[[j]] <- tmp15
    }
    tmp17[[i]] <- tmp16
  }

  for (i in 1:length(rate1)) {
    for (j in 1:5) {
      for (z in 1:length(window1)) {
        data <- tmp12 %>% filter(type == paste0("run",j)) %>% filter(traits %in% tmp17[[i]][[j]][[z]])
        for (v in 1:length(date1)) {
          data1 <- data %>% filter(sub.dates <= date1[v])
          writeXStringSet(seq[names(seq) %in% data1$traits], paste0("combine_sampling/combine_",rate1[i],"_",window1[z],"_",date1[v],"_sim",q,"_", j ,".fasta"), format="fasta")
          write.table(epi_data %>% filter(traits %in% data1$traits), paste0("combine_sampling/combine_",rate1[i],"_",window1[z],"_",date1[v],"_sim",q,"_", j ,".txt"), quote = FALSE, sep = "\t")
        }
      }
    }
  }
}
