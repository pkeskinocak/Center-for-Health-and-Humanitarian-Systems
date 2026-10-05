#optimization with penalties

#libraries
library(tidyr)
library(dplyr)
library(withr)
library(parallel)
library(stats)

#lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.2/"

dir = "/storage/coda1/p-pk50/0/awakabayashi3/cruise_ship/optimization/"

#functions
source(paste0(dir, "run/inputs/functions_cs_part1.R"))
source(paste0(dir, "run/inputs/functions_cs_part2.R"))
source(paste0(dir, "run/inputs/cs_cal_function.R"))

#inputs
load(paste0(dir, "run/inputs/opt_data.Rdata"))

#objective function
objective_function <- function(parameters) {
  
  num_cores <- 10
  cl <- makeCluster(num_cores)
  
  # Define or source your function on each node
  clusterEvalQ(cl, {
    dir = "/storage/coda1/p-pk50/0/awakabayashi3/cruise_ship/optimization/"
    #functions
    source(paste0(dir, "run/inputs/functions_cs_part1.R"))
    source(paste0(dir, "run/inputs/functions_cs_part2.R"))
    source(paste0(dir, "run/inputs/cs_cal_function.R"))
    #libraries
    library("dplyr")
    library("tidyr")
    library("parallel")
    #data
    load(paste0(dir, "run/inputs/opt_data.Rdata"))
  })
  
  #clusterSetRNGStream(cl = cl, 2024)
  
  results <- parSapply(cl, 1:10, function(x) simulate_voyage_opt(net = std_network_cal, initial = initial_std_network_cal,
                                                                cs_input_cal = cs_input_top,  covid_input = covid_input,
                                                                network_input = network_input, observed_data = observed_data,
                                                                parameters_opt = parameters))
  stopCluster(cl)
  
  results_mean <- mean(results)
  
  return(results_mean)
}

# Penalty function for constraints
penalty_function <- function(parameters) {
  
  penalty <- 0
  
  #x1 >= x2
  if (parameters[1] < parameters[2]) {
    penalty <- penalty + (parameters[2] - parameters[1])^2
  }
  
  #x1 >= x3
  if (parameters[1] < parameters[3]) {
    penalty <- penalty + (parameters[3] - parameters[1])^2
  }
  
  #x1 >= x4
  if (parameters[1] < parameters[4]) {
    penalty <- penalty + (parameters[4] - parameters[1])^2
  }
  
  #x3 >= x2
  if (parameters[3] < parameters[2]) {
    penalty <- penalty + (parameters[2] - parameters[3])^2
  }
  
  #x3 >= x4
  if (parameters[3] < parameters[4]) {
    penalty <- penalty + (parameters[4] - parameters[3])^2
  }
  
  return(penalty)
}


# Combined objective and penalty function
combined_function <- function(parameters) {
  return(objective_function(parameters) + 100 * penalty_function(parameters))  # Penalty factor 100
}


#run
result_A3 <- optim( par = initial_guess,  
                    fn = combined_function,
                    gr = NULL,
                    method = "Nelder-Mead", 
                    control = list(maxit = 2, abstol = 5, reltol = NULL,
                                   trace = TRUE))



save(result_A3, file = paste0(dir, "results/result_A3.rdata"))