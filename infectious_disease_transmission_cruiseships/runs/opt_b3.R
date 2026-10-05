#optimization with penalties

#libraries
library(tidyr)
library(dplyr)
library(withr)
library(parallel)
library(stats)

#lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.2/"

dir = "code/"

#functions
source(paste0(dir, "functions_cs_part1.R"))
source(paste0(dir, "functions_cs_part2.R"))
source(paste0(dir, "calibration/cs_cal_function_p2.R"))

#inputs
load(paste0(dir, "calibration/opt_data.Rdata"))

#objective function
objective_function <- function(parameters) {
  
  num_cores <- detectCores() -1
  cl <- makeCluster(num_cores)
  
  # Define or source your function on each node
  clusterEvalQ(cl, {
    dir = "code/"
    #functions
    source(paste0(dir, "functions_cs_part1.R"))
    source(paste0(dir, "functions_cs_part2.R"))
    source(paste0(dir, "calibration/cs_cal_function_p2.R"))
    #libraries
    library("dplyr")
    library("tidyr")
    library("parallel")
    #data
    load(paste0(dir, "calibration/opt_data.Rdata"))
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

  return(penalty)
}


# Combined objective and penalty function
combined_function <- function(parameters) {
  return(objective_function(parameters) + 100 * penalty_function(parameters))  # Penalty factor 100
}

initial_guess <- initial_guess[c(1,2,5)]


#run
initial <- Sys.time()
result_B3 <- optim( par = initial_guess,  
                    fn = combined_function,
                    gr = NULL,
                    method = "Nelder-Mead", 
                    control = list(maxit = 100, abstol = 5, reltol = NULL,
                                   trace = TRUE))

final <- Sys.time()


save(result_B3, file = paste0(dir, "calibration/optimization/result_B3.rdata"))

