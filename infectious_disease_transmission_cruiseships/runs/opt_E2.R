#2

#optimization without constraints
#variables: (1) Pmain, (2) Pact, (3) Mkids, (4) Pcp = Pact, (5) Qdecrease


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
  source(paste0(dir,"run/inputs/cs_cal_function_v2.R"))
  
  
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
  source(paste0(dir,"calibration/cs_cal_function_v2.R"))
  
  
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
      source(paste0(dir,"run/inputs/cs_cal_function_v2.R"))
      
      
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
      source(paste0(dir,"calibration/cs_cal_function_v2.R"))
      
      
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

initial_guess_E2 <- initial_guess[c(1,2,3,5)]
#run
result_E2 <- optim( par = initial_guess_E2,  
                    fn = objective_function,
                    gr = NULL,
                    method = "Nelder-Mead", 
                    control = list(maxit = 100, abstol = 5, reltol = NULL,
                                   trace = TRUE))


if (server == TRUE) {
  dir = "/storage/coda1/p-pk50/0/awakabayashi3/cruise_ship/optimization/"
  
  save(result_E2, file = paste0(dir, "results/result_E2_cluster_100rep.rdata"))
  
  
} else {
  
  dir = "code/"
  
  save(result_E2, file = paste0(dir, "calibration/optimization/results/result_E2_local_100rep.rdata"))
  
}



