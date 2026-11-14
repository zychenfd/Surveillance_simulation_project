#==1. Load packages==
library(stringr)
library(dplyr)
library(seqinr)
library(lubridate)

#==2. Read data==
#non-antigenic sites
ref_1 <- read.fasta("ref_non_antigenic.fasta", forceDNAtolower = FALSE, set.attributes = FALSE)[[1]]
mu1 <- 1e-2/365
nonantigenic_sites <- c((1*3-2):(86*3),
                        (93*3-2):(140*3),
                        (143*3-2):(153*3),
                        (160*3-2):(169*3),
                        (175*3-2):(175*3),
                        (182*3-2):(182*3),
                        (188*3-2):(200*3),
                        (213*3-2):(219*3),
                        (223*3-2):(237*3),
                        (240*3-2):(251*3),
                        (255*3-2):(566*3))

#antigenic sites
ref_2 <- read.fasta("ref_antigenic.fasta", forceDNAtolower = FALSE, set.attributes = FALSE)[[1]]
mu2 <- 1.5e-2/365
antigenic_sites <- c((87*3-2):(92*3),
                     (141*3-2):(142*3),
                     (154*3-2):(159*3),
                     (170*3-2):(174*3),
                     (176*3-2):(181*3),
                     (183*3-2):(187*3),
                     (201*3-2):(212*3),
                     (220*3-2):(222*3),
                     (238*3-2):(239*3),
                     (252*3-2):(254*3))

#other parameters
root.date <- as.Date("2009-01-01")
kappa <- 10
regions = c("Canada", "Caribbean", "Central America", "Mexico", "United States")
initPatch <- 4

#==3. Sequence simulation==
#Taking simulation 100 as an example
for (k in 100) {
  #==1. Read epi simulation data==
  sim <- readRDS(paste0("../1.epi_model/res_v",k,".RDS"))
  
  #==2. Check epidemic size==
  sim$size
  
  #==3. Convert cases.info from a list to a data frame==
  cases.info <- data.frame(do.call("rbind", sim$sourcing$data)) 
  names(cases.info) <- c("patch.source",
                         "ID.source",
                         "infection.time",
                         "latent.period",
                         "infectious.time", 
                         "infectious.period",
                         "serial.time",
                         "stop",
                         "hard.stop")
  
  cases.info$id <- row.names(cases.info)
  cases.info$source <- paste0(cases.info$patch.source, "_", cases.info$ID.source)

  #==4. Number of initial sequences==
  n.init <- length(which(cases.info$ID.source < 0))
  names.init <- rownames(cases.info[which(cases.info$ID.source < 0),])
  
  #==5. Information from the reference sequence==
  #==For non-antigenic sites==
  genome.length <- length(ref_1)
  baseComp <- as.numeric(table(ref_1)) 
  baseFreq <- baseComp / genome.length
  bases <- names(table(ref_1))
  names(baseComp) <- names(baseFreq) <- bases
  
  #====For antigenic sites====
  genome.length2 <- length(ref_2)
  baseComp2 <- as.numeric(table(ref_2)) 
  baseFreq2 <- baseComp2 / genome.length2
  bases2 <- names(table(ref_2))
  names(baseComp2) <- names(baseFreq2) <- bases2
  
  #==6. Create lists, with first case deposited== 
  snp <- vector("list", nrow(cases.info))
  names(snp) <- rownames(cases.info)
  snp[[names.init[1]]] <- list(loc1 = c(), 
                               loc2 = c(), 
                               nuc1 = ref_1,
                               nuc2 = ref_2,
                               date = decimal_date(as.Date(root.date)), 
                               nb.snp1 = 0,
                               nb.snp2 = 0,
                               d.dist1 = 0,
                               d.dist2 = 0)
  
  # Time max can exceed the time.max used for the simulation if the epidemic doesn't stop before the end of the simulation 
  # To speed up the evolution loop, we define only the time steps that should be treated
  time.vector <- sort(unique(cases.info$infection.time))
  time.vector <- time.vector[which(time.vector >= 0)] 
  
  #==7. EVOLUTIONARY PARAMETERS OF THE HKY MODEL==
  #==For non-antigenic sites==
  #7.1 HKY parameters
  beta <- mu1/(2*(baseFreq['a']*baseFreq['c']+baseFreq['a']*baseFreq['t']+baseFreq['c']*baseFreq['g']+baseFreq['g']*baseFreq['t'])+
                 2*kappa*(baseFreq['a']*baseFreq['g']+baseFreq['c']*baseFreq['t']))
  alpha <- kappa*beta
  #7.2 rate matrix Q
  Q <- matrix(NA, ncol = 4, nrow = 4)
  Q[upper.tri(Q)] <- c(beta*baseFreq['c'], alpha*baseFreq['g'], 
                       beta*baseFreq['g'], beta*baseFreq['t'], 
                       alpha*baseFreq['t'], beta*baseFreq['t']) 
  Q[lower.tri(Q)] <- c(beta*baseFreq['a'], alpha*baseFreq['a'], 
                       beta*baseFreq['a'], beta*baseFreq['c'], 
                       alpha*baseFreq['c'], beta*baseFreq['g'])
  diag(Q) <- -apply(Q, 1, sum, na.rm = TRUE)
  #7.3 Compute invertible matrix of Q and its eigenvalues
  E <- eigen(Q)$vectors  
  E_1 <- solve(E)
  L <- diag(eigen(Q)$values)
  
  #====For antigenic sites====
  #7.1 HKY parameters
  beta2 <- mu2/(2*(baseFreq2['a']*baseFreq2['c']+baseFreq2['a']*baseFreq2['t']+baseFreq2['c']*baseFreq2['g']+baseFreq2['g']*baseFreq2['t'])+
                  2*kappa*(baseFreq2['a']*baseFreq2['g']+baseFreq2['c']*baseFreq2['t']))
  alpha2 <- kappa*beta2
  #7.2 rate matrix Q
  Q2 <- matrix(NA, ncol = 4, nrow = 4)
  Q2[upper.tri(Q2)] <- c(beta2*baseFreq2['c'], alpha2*baseFreq2['g'], 
                         beta2*baseFreq2['g'], beta2*baseFreq2['t'], 
                         alpha2*baseFreq2['t'], beta2*baseFreq2['t']) 
  Q2[lower.tri(Q2)] <- c(beta2*baseFreq2['a'], alpha2*baseFreq2['a'], 
                         beta2*baseFreq2['a'], beta2*baseFreq2['c'], 
                         alpha2*baseFreq2['c'], beta2*baseFreq2['g'])
  diag(Q2) <- -apply(Q2, 1, sum, na.rm = TRUE)
  #7.3 Compute invertible matrix of Q and its eigenvalues
  E2 <- eigen(Q2)$vectors  
  E_1_2 <- solve(E2)
  L2 <- diag(eigen(Q2)$values)
  
  #==8. SEQUENCE EVOLUTION==
  evol <- function(x, t) {

    # Retrieve individual parameters
    vec <- cases.info[rownames(cases.info) == x, ]
    
    # Get the Serial interval 
    serial.time <- as.numeric(vec["serial.time"])
    
    # Compute the diagonal matrix for the current generation time
    L_curr <- L
    diag(L_curr) <- exp(diag(L)*serial.time)
    
    L_curr2 <- L2
    diag(L_curr2) <- exp(diag(L2)*serial.time)
    
    # Compute the transition probability matrix 
    p_t <- E %*% L_curr %*% solve(E)
    colnames(p_t) <- rownames(p_t) <- bases
    if (!all(apply(p_t, 1, sum) > rep(1-1e-10, length(bases)))) stop("Transition probabilities do not sum to 1!")
    
    p_t2 <- E2 %*% L_curr2 %*% solve(E2)
    colnames(p_t2) <- rownames(p_t2) <- bases2
    if (!all(apply(p_t2, 1, sum) > rep(1-1e-10, length(bases2)))) stop("Transition probabilities do not sum to 1!")
    
    # Update the composition of the sequence
    ID.infector <- paste(vec[c('patch.source', 'ID.source')], collapse = "_")
    
    if (ID.infector == names.init[1]){
      
      # Sample new sequence
      newSequence1 <- sapply(ref_1, function(x) sample(bases,
                                                       size = 1,
                                                       prob = p_t[x, ]))
      newSequence2 <- sapply(ref_2, function(x) sample(bases2,
                                                       size = 1,
                                                       prob = p_t2[x, ]))
      
      # Store the changes
      loc1 <- which(newSequence1 != ref_1)
      loc2 <- which(newSequence2 != ref_2)
      nuc1 <- newSequence1[loc1]
      nuc2 <- newSequence2[loc2]
      d.dist1 <- length(loc1)
      d.dist2 <- length(loc2)
      
    } else {
      
      # Sequence of the source
      source1 <- ref_1
      source2 <- ref_2
      if (length(snp[[ID.infector]]$loc1) > 0) source1[snp[[ID.infector]]$loc1] <- snp[[ID.infector]]$nuc1
      if (length(snp[[ID.infector]]$loc2) > 0) source2[snp[[ID.infector]]$loc2] <- snp[[ID.infector]]$nuc2
      
      # Sample new sequence
      newSequence1 <- sapply(source1, function(x) sample(bases,
                                                         size = 1,
                                                         prob = p_t[x, ]))
      newSequence2 <- sapply(source2, function(x) sample(bases2,
                                                         size = 1,
                                                         prob = p_t2[x, ]))
      
      # SNPs that differ from the reference genome are stored
      
      loc1 <- which(newSequence1 != ref_1) 
      loc2 <- which(newSequence2 != ref_2) 
      nuc1 <- newSequence1[loc1]
      nuc2 <- newSequence2[loc2]
      d.dist1 <- length(which(newSequence1 != source1))  # Compute the number of new SNPs compared to the reference sequence 
      d.dist2 <- length(which(newSequence2 != source2))  # Compute the number of new SNPs compared to the reference sequence 
    }
    
    # root.date should be expressed as "yyyy-mm-dd", ie "%Y-%m-%d" or "%Y/%m/%d"
    final.date <- decimal_date(as.Date(root.date) + as.numeric(vec['infectious.time']))
    final.date <- round(final.date, 6)
    
    # Output
    res <- list(nuc1 = nuc1, nuc2 = nuc2, loc1 = loc1, loc2 = loc2, date = final.date,
                nb.snp1 = length(loc1), nb.snp2 = length(loc2), d.dist1 = d.dist1, d.dist2 = d.dist2)
    return(res)
  }
  
  # Sequence evolution by time step
  set.seed(100)
  start <- Sys.time()
  for (t in time.vector) {
    print(t)
    id <- rownames(cases.info[cases.info$infection.time == t, ])
    print(id)
    snp[id] <- lapply(id, evol, t = t)
  }
  end <- Sys.time()
  
  #==9. WRITING SEQUENCES==
  #==Sampling time, name, location of sequences==
  samp.dates <- sapply(names(snp), function(x) snp[[x]]$date)
  mutation1 <- sapply(names(snp), function(x) snp[[x]]$nb.snp1)
  mutation2 <- sapply(names(snp), function(x) snp[[x]]$nb.snp2)
  table(mutation1)
  table(mutation2)
  long.names <- paste(names(snp), samp.dates, sep = "_")
  short.names <- paste0("region", names(snp))
  patch <- sapply(long.names, function(x) as.numeric(strsplit(x, "_")[[1]][1])) + 1
  long.names <- paste0("region", long.names)
  
  #==Write text file containing the spatial and temporal metadata==
  directory <- ""
  simulationNb <- k
  continuous_traits <- data.frame(traits = long.names, patch = patch, samp.dates = samp.dates, short.names = short.names)
  write.table(continuous_traits, paste0(directory, "sim", simulationNb, "_metadata.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  metadata <- continuous_traits 
  metadata$regions <- sapply(metadata$patch, function(x) regions[x])
  table(metadata$regions)
  index.case <- paste0(initPatch-1, "_0")
  names.init[1] == index.case
  
  #==10. FUNCTION TO WRITE SEQUENCES==
  
  write.seq <- function(case, index.case, snp, error.rate = 0) {
    # Write full sequence based on reference sequence
    out1 <- snp[[index.case]]$nuc1
    out1[snp[[case]]$loc1] <- unlist(snp[[case]]$nuc1)
    
    out2 <- snp[[index.case]]$nuc2
    out2[snp[[case]]$loc2] <- unlist(snp[[case]]$nuc2)
    
    out <- rep("N", 1698)
    out[nonantigenic_sites] <- out1
    out[antigenic_sites] <- out2
    print(table(out))
    
    return(out)
  }
  
  set.seed(100)
  sequences <- lapply(str_remove_all(metadata$short.names, "region"), 
                      write.seq, 
                      index.case = index.case, 
                      snp = snp)
  table(duplicated(sequences))
  write.fasta(sequences, 
              metadata$traits,
              paste0(directory, "sim", simulationNb, "_sequences", ".fasta"))
  
  traits <- metadata[, c("traits", "regions")]
  write.table(traits, 
              paste0(directory, "sim", simulationNb, "_traits", ".txt"),
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
  
  # Convert epidemic output to transmission chain in data frame
  transmission.chain <- data.frame(do.call("rbind", sim$sourcing$data)) 
  names(transmission.chain) <- c("patch.source",
                                 "ID.source",
                                 "infection.time",
                                 "latent.period",
                                 "infectious.time", 
                                 "infectious.period",
                                 "serial.time",
                                 "stop",
                                 "hard.stop")
  
  # Add the case and the source for each infection event
  transmission.chain$case <- rownames(transmission.chain)
  transmission.chain$source <- paste(transmission.chain$patch.source, transmission.chain$ID.source, sep = "_")
  
  # Modify the patch.source counting
  transmission.chain$patch.source <- sapply(transmission.chain$patch.source, function(x) 
    ifelse(x < 0, x, x + 1))
  transmission.chain <- left_join(transmission.chain %>% mutate(case = paste0("region",case)), metadata, by = c("case" = "short.names"))
  
  # Write the transmission chain
  write.table(transmission.chain, paste0(directory, "sim", simulationNb, "_transmission_chain.txt"), quote = FALSE, sep = "\t")
}
