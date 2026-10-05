#E1

#optimization without constraints
#variables: (1) Pmain, (2) Pact, (3) Mkids, (4) Pcp, (5) Qdecrease


#Load dataset

server = FALSE

if (server == TRUE) {
  
  #libraries
  library(tidyr)
  library(dplyr)
  library(withr)
  library(parallel)

  dir = "/storage/coda1/p-pk50/0/awakabayashi3/cruise_ship/optimization/"
  
  #source functions
  source(paste0(dir, "run/inputs/functions_cs_part1.R"))
  source(paste0(dir,"run/inputs/functions_cs_part2.R"))
  source(paste0(dir,"run/inputs/cs_cal_function.R"))
  
  
  #load inputs
  load(paste0(dir,"run/inputs/opt_data.rdata"))
  
} else {
  
  #libraries
  library(tidyr)
  library(dplyr)
  library(withr)
  library(parallel)
  
  dir = "code/"
  
  #source functions
  source(paste0(dir,"functions_cs_part1.R"))
  source(paste0(dir,"functions_cs_part2.R"))
  source(paste0(dir,"calibration/cs_cal_function.R"))
  
  
  #load inputs
  load(paste0(dir,"calibration/optimization/inputs/opt_data.rdata"))
}


#objective function
objective_function <- function(parameters) {
  
  num_cores <- 6
  cl <- makeCluster(num_cores)
  
  # Define or source your function on each node
  clusterEvalQ(cl, {
    server = FALSE
    
    if (server == TRUE) {
      
      #libraries
      library(tidyr)
      library(dplyr)
      library(withr)
      library(parallel)
      
      dir = "/storage/coda1/p-pk50/0/awakabayashi3/cruise_ship/optimization/"
      
      #source functions
      source(paste0(dir, "run/inputs/functions_cs_part1.R"))
      source(paste0(dir,"run/inputs/functions_cs_part2.R"))
      source(paste0(dir,"run/inputs/cs_cal_function.R"))
      
      
      #load inputs
      load(paste0(dir,"run/inputs/opt_data.rdata"))
      
    } else {
      
      #libraries
      library(tidyr)
      library(dplyr)
      library(withr)
      library(parallel)
      
      dir = "code/"
      
      #source functions
      source(paste0(dir,"functions_cs_part1.R"))
      source(paste0(dir,"functions_cs_part2.R"))
      source(paste0(dir,"calibration/cs_cal_function.R"))
      
      
      #load inputs
      load(paste0(dir,"calibration/optimization/inputs/opt_data.rdata"))
    }
    
  })
  
  #clusterSetRNGStream(cl = cl, 2024)
  
  results <- parSapply(cl, 1:100, function(x) simulate_voyage_opt(net = std_net_cal, initial = initial_std_net_cal,
                                                                 cs_input_cal = cs_input_top,  covid_input = covid_input,
                                                                 network_input = network_input, observed_data = observed,
                                                                 parameters_opt = parameters))
  stopCluster(cl)
  
  results_mean <- mean(results)
  
  return(results_mean)
}


#run
result_E1 <- optim( par = initial_guess,  
                    fn = objective_function,
                    gr = NULL,
                    method = "Nelder-Mead", 
                    control = list(maxit = 100, abstol = 5, reltol = NULL,
                                   trace = TRUE))


if (server == TRUE) {
  dir = "/storage/coda1/p-pk50/0/awakabayashi3/cruise_ship/optimization/"
  
  save(result_E1, file = paste0(dir, "results/result_E1_cluster.rdata"))
  
  
} else {
  
  dir = "code/"
  
  save(result_E1, file = paste0(dir, "calibration/optimization/results/result_E1_local_100rep.rdata"))
  
}



#Run test

Pmain <- result_E1$par[1]
Pact <- result_E1$par[2]
Pkids <- result_E1$par[2]*(1+result_E1$par[3])
Pemp_cc <- result_E1$par[2]
Pemp_cp <-  result_E1$par[4]
Q_decrease <- result_E1$par[5]

cs_input_test <- cs_input_top %>%
  mutate("Crew Contagion Probability - In network" = Pmain,
         "Passenger Contagion Probability - In network"  = Pmain,
         "Pass-Pass Dining Contagion Probability" = Pact,
         "Pass-Pass Entertainment Contagion Probability" = Pact,
         "Pass-Pass Nightlife Contagion Probability" = Pact,
         "Pass-Pass Kids Contagion Probability" = Pkids,
         "Crew-Crew Food Contagion Probability" = Pemp_cc,
         "Crew-Crew Housekeeping Contagion Probability" = Pemp_cc,
         "Crew-Crew Entertainment Contagion Probability" = Pemp_cc,
         "Crew-Passenger Food Contagion Probability" = Pemp_cp,
         "Crew-Passenger Entertainment Contagion Probability" = Pemp_cp,
         "Crew-Passenger Housekeeping Contagion Probability" =  Pemp_cp,
         "Decrease in prob transmission" = Q_decrease)
  
test <- simulate_voyage_onerep(net = std_net_cal , initial =initial_std_net_cal,
                       cs_input = cs_input_test,
                       scenario = 1,
                       covid_input = covid_input,
                       network_input = network_input,
                       observed_data = observed,
                       plotNetwork = FALSE,
                       detailed_output = FALSE, stop_early = TRUE)
