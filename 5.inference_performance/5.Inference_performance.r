#==Load libraries==
library(readr)
library(stringr)
library(dplyr)
library(scoringRules)
library(RevGadgets)
library(coda)
library(lubridate)
library(future.apply)
library(parallel)
library(ggtree)
library(treeio)
library(tidyr)
library(phangorn)

#==1. root location inference==
t_data <- read_log("phylogeographic_output/seq_homo_sim100_1_spa.geo.states.log", skip = 4)
t_data$source <- sapply(str_split(gsub('["]', '', t_data$X1), "\t"), function(x) x[2])
t_data <- t_data[round(nrow(t_data)*0.1, 0):nrow(t_data),]
root_location <- t_data %>% 
  group_by(source) %>% 
  summarise(n = n()) %>%
  mutate(total = sum(n)) %>%
  mutate(post_support = n/total) %>%
  filter(source == "region3") # region3 is Mexico

#==2. tMRCA inference==
summarise_hpd_lower <- function(x) {
  if(length(x) <= 1) {
    return(x[1]);
  }
  return(HPDinterval(as.mcmc(x),prob = 0.95)[1])
}

summarise_hpd_upper <- function(x) {
  if(length(x) <= 1) {
    return(x[1]);
  }
  return(HPDinterval(as.mcmc(x),prob = 0.95)[2])
}

tMRCA_data <- readTrace("phylogenetic_output/seq_homo_sim100_1.log", burnin = 0.1)[[1]]
tMRCA_data <- tMRCA_data[,6]

mean_value <- mean(tMRCA_data)
median_value <- median(tMRCA_data)
low <- summarise_hpd_lower(tMRCA_data)
upp <- summarise_hpd_upper(tMRCA_data)
min <- min(tMRCA_data)
max <- max(tMRCA_data)
es <- es_sample(y = 2009, dat = matrix(tMRCA_data, nrow = 1)) #Energy Score metric

#==3. Ancestral geographic state inference for pairs of sampled tips==
ncores <- detectCores() - 4
plan(multisession, workers = ncores)

#define function
build_phylo_from_transmission <- function(df) {
  all_nodes <- unique(c(df$child, df$parent))
  roots <- setdiff(df$parent, df$child)
  tips <- setdiff(df$child, df$parent)
  internal <- setdiff(all_nodes, c(tips, roots))
  internal <- c(roots[1], internal)
  n_tip <- length(tips)
  n_node <- length(internal)
  tip_index <- setNames(seq_along(tips), tips) 
  node_index <- setNames(seq_len(n_node) + n_tip, internal)
  
  edges <- t(apply(df, 1, function(row) {
    child <- row["child"]
    parent <- row["parent"]
    child_num <- if (child %in% tips) tip_index[child] else node_index[child]
    parent_num <- if (parent %in% tips) tip_index[parent] else node_index[parent]
    c(parent_num, child_num)
  }))
  colnames(edges) <- NULL
  
  phy <- list(
    edge = edges,
    tip.label = tips,
    node.label = internal,
    Nnode = n_node,
    edge.length = rep(1, nrow(edges)) 
  )
  class(phy) <- "phylo"

  phy
}

#converted the transmission chain into a transmission tree containing polytomies
for (i in c(100)) {
  # i = 100
  tmp <- read.delim(paste0("../2.seq_generation/sim",i,"_transmission_chain.txt")) %>% 
    select(c("case","source")) %>%
    mutate(source = paste0("region",source)) %>%
    rename("child" = "case", "parent" = "source")
  epi_tree_tmp <- build_phylo_from_transmission(tmp)
  write.tree(epi_tree_tmp, paste0("epi_tree_sim",i,".tre"))
}

#mcc files (with one simulation as example)
files2 <- c(list.files(path = "phylogeographic_output/mcc_tree", full.names = T))

pair_source <- c()
for (j in c(100)) {
  print(j)
  
  #1. Epi Data
  epi_data <- read.delim(paste0("../2.seq_generation/sim",j,"_transmission_chain.txt")) %>% 
    mutate(sampling.time = as.Date("2009-01-01") + infectious.time + 3)
  
  epi_tree <- read.tree(paste0("epi_tree_sim",j,".tre"))
  epi_p1 <- ggtree(epi_tree)
  epi_tmp <- epi_p1$data
  
  files3 <- files2[str_detect(files2, paste0("sim",j))]
  
  taxa_homo <- files3[str_detect(files3, "homo")]
  tree_homo <- read.beast(taxa_homo)
  pairs_df_homo <- as.data.frame(t(combn(tree_homo@phylo$tip.label,2))) %>%
    mutate(V3 = paste0( sapply(str_split(V1,"_"), function(x) x[1]) , "_", sapply(str_split(V1,"_"), function(x) x[2]) ), 
           V4 = paste0( sapply(str_split(V2,"_"), function(x) x[1]) , "_", sapply(str_split(V2,"_"), function(x) x[2]) ), 
           epi_source = NA) 
  pairs_df_homo$epi_source <- future_apply(pairs_df_homo[, 3:4], 1, function(x) MRCA(epi_tree, as.character(x)))
  # saveRDS(pairs_df_homo, paste0("pairs_true_sim",j,".RDS"))
  
  dat_tmp <- c()
  for (i in 1:length(files3)) {
    # i = 2
    file_name = strsplit(files3[i],"/")[[1]][length(strsplit(files3[i],"/")[[1]])]
    n_sim = rev(strsplit(file_name,"_")[[1]])[3]
    
    #2. MCC tree data
    tree <- read.beast(files3[i])
    date_max <- decimal2Date(max(as.numeric(sapply(str_split(tree@phylo$tip.label, "_"), function(x) x[3])))) + 3
    ggtree(tree, mrsd = date_max)+ theme_tree2() -> fig
    t_data <- as.data.frame(fig$data)
    
    pairs_df <- as.data.frame(t(combn(tree@phylo$tip.label,2)))
    pairs_df <- pairs_df %>%
      mutate(V3 = paste0( sapply(str_split(V1,"_"), function(x) x[1]) , "_", sapply(str_split(V1,"_"), function(x) x[2]) ), 
             V4 = paste0( sapply(str_split(V2,"_"), function(x) x[1]) , "_", sapply(str_split(V2,"_"), function(x) x[2]) )) 
    
    MRCA_infer <- mrca(tree@phylo)
    pairs_df$infer_source <- MRCA_infer[as.matrix(pairs_df[, c("V1", "V2")])]
    pairs_df <- left_join(pairs_df, pairs_df_homo)
    print(table(is.na(pairs_df)))
    
    pairs_df1 <- left_join(pairs_df, epi_tmp[,c(2,4)], by = c("epi_source" = "node")) %>%
      left_join(epi_data[,c(3,5,10,16)], by = c("label" = "case")) %>%
      mutate(infection.time = as.Date("2009-01-01") + infection.time,
             infectious.time = as.Date("2009-01-01") + infectious.time) %>%
      mutate(true_geo = sapply(str_split(label, "_"), function(x) x[1])) %>%
      left_join(t_data[,c(2,5,21)], by = c("infer_source" = "node")) %>%
      mutate(tmrca = as.Date(date_decimal(x))) %>%
      select(-"x") %>%
      mutate(sim = paste0("sim",j),
             type = str_remove_all(file_name,"_mcc.tre"))
    
    dat_tmp <- rbind(dat_tmp, pairs_df1)
  }
  # saveRDS(dat_tmp, paste0("pairs_output_sim",j,".RDS"))
  pair_source <- rbind(pair_source, dat_tmp)
}
pair_source <- rbind(readRDS("pairs_output_sim100.RDS"))

#define match types
pair_source$match_geo <- pair_source$true_geo == pair_source$geo
pair_source$pair_type <- substr(pair_source$V1,1,7) == substr(pair_source$V2,1,7)
pair_source$node_type <- "deep"
pair_source$node_type[pair_source$sampling.time >= as.Date("2009-04-01")] <- "recent"

#Deep and recent pairs of tips
pair_dis1 <- pair_source %>%
  group_by(sim, type, node_type, match_geo) %>%
  summarise(n = n()) 

mmp1 <- data.frame(type = rep(unique(pair_dis1$type), each = 4),
                   node_type = rep(rep(c("deep", "recent"), each = 2), length(unique(pair_dis1$type))),
                   match_geo = rep(rep(c(FALSE, TRUE), 2), length(unique(pair_dis1$type))))

pair_dis1 <- right_join(pair_dis1, mmp1)
pair_dis1 <- pair_dis1 %>%
  mutate(n = ifelse(is.na(n),0,n)) %>%
  mutate(sim = ifelse((is.na(sim) & str_detect(type, "homo") ), sapply(str_split(type,"\\_"), function(x) x[3]), sim)) %>%
  mutate(sim = ifelse((is.na(sim) & str_detect(type, "combine") ), sapply(str_split(type,"\\_"), function(x) x[5]), sim)) %>%
  mutate(sim = ifelse((is.na(sim) & str_detect(type, "seq") ), sapply(str_split(type,"\\_"), function(x) x[4]), sim)) %>%
  group_by(sim, type, node_type) %>%
  mutate(total = sum(n)) %>%
  filter(match_geo == T) %>%
  mutate(accuracy = n*100/total)

#Between-region and within-region pairs of tips
pair_source1 <- pair_source %>%
  group_by(sim, type, pair_type, match_geo) %>%
  summarise(n = n()) 

mmp <- data.frame(type = rep(unique(pair_source1$type), each = 4),
                  pair_type = rep(rep(c(FALSE, TRUE), each = 2), length(unique(pair_source1$type))),
                  match_geo = rep(rep(c(FALSE, TRUE), 2), length(unique(pair_source1$type))))

pair_source1 <- right_join(pair_source1, mmp)
pair_source1 <- pair_source1 %>%
  mutate(n = ifelse(is.na(n),0,n)) %>%
  mutate(sim = ifelse((is.na(sim) & str_detect(type, "sim75") ), "sim75", sim)) %>%
  mutate(sim = ifelse((is.na(sim) & str_detect(type, "sim51") ), "sim51", sim)) %>%
  group_by(sim, type, pair_type) %>%
  mutate(total = sum(n)) %>%
  filter(match_geo == T) %>%
  mutate(accuracy = n*100/total)
