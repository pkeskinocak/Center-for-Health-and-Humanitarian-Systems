#optimization with constraints 

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


#run
result_A1 <- optim( par = initial_guess,  
                    fn = objective_function,
                    gr = NULL,
                    method = "Nelder-Mead", 
                    control = list(maxit = 100, abstol = 5, reltol = NULL,
                                        trace = TRUE))



save(result_A1, file = paste0(dir, "results/result_A1.rdata"))