library(Rcpp)

#==Load epidemic model==
sourceCpp("epi_model.cpp") # Epidemic model

#==Parameters==
initPatch <- 4
regions <- c("Canada", "Caribbean", "Central America", "Mexico", "United States")
human_pop <- c(33593920, 41142850, 42835202, 111049420, 304725917) # for each region
n.patch <- length(human_pop)
init.inf <- rep(0, n.patch)
init.inf[initPatch] <- 1
beta <- 0.411/human_pop
connection_strength <- 60
time.max <- as.integer(as.Date("2009-06-11")-as.Date("2009-01-01"))
min.cases <- 30000
period_infectious <- 3.38
period_latent <- 2.62

#==Simulation==
#Note: Demo matrices data of human movement data between these regions have be deposited in the folder "./demo_matrices_data/"
#The original air travel data are proprietary and confidential and cannot be transferred as per the Terms of Service.
set.seed(100)
for (i in 1:100) {
  res <- cpp_model(human_pop, beta, init.inf, time.max, period_infectious, period_latent, connection_strength)
  while (res$size < min.cases) {
    print(i)
    res <- cpp_model(human_pop, beta, init.inf, time.max, period_infectious, period_latent, connection_strength)
  }
  saveRDS(res, paste0("simulation",i,".RDS"))
}
