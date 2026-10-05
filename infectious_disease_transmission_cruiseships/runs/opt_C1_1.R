#optimization with constraints 

#libraries
library(tidyr)
library(dplyr)
library(withr)
library(parallel)
library(stats)

#lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.2/"

#local
dir = "code/"

#functions
source(paste0(dir, "functions_cs_part1.R"))
source(paste0(dir, "functions_cs_part2.R"))
source(paste0(dir, "calibration/cs_cal_function_p3.R"))

#inputs
load(paste0(dir, "calibration/opt_data.Rdata"))

#server
#dir = "/storage/coda1/p-pk50/0/awakabayashi3/cruise_ship/optimization/"
#source(paste0(dir, "run/inputs/functions_cs_part1.R"))
#source(paste0(dir, "run/inputs/functions_cs_part2.R"))
#source(paste0(dir, "run/inputs/cs_cal_function_p3.R"))

#inputs
#load(paste0(dir, "run/inputs/opt_data.Rdata"))

#objective function
objective_function <- function(parameters) {
  
  num_cores <- 6
  cl <- makeCluster(num_cores)
  
  # Define or source your function on each node
  clusterEvalQ(cl, {
    #server
    # dir = "/storage/coda1/p-pk50/0/awakabayashi3/cruise_ship/optimization/"
    # #functions
    # source(paste0(dir, "run/inputs/functions_cs_part1.R"))
    # source(paste0(dir, "run/inputs/functions_cs_part2.R"))
    # source(paste0(dir, "run/inputs/cs_cal_function.R"))
    # #libraries
    # library("dplyr")
    # library("tidyr")
    # library("parallel")
    # #data
    # load(paste0(dir, "run/inputs/opt_data.Rdata"))
    
    #local
    dir = "code/"
    #functions
    source(paste0(dir, "functions_cs_part1.R"))
    source(paste0(dir, "functions_cs_part2.R"))
    source(paste0(dir, "calibration/cs_cal_function_p3.R"))
    #libraries
    library("dplyr")
    library("tidyr")
    library("parallel")
    #inputs
    load(paste0(dir, "calibration/opt_data.Rdata"))
  })
  
  #clusterSetRNGStream(cl = cl, 2024)
  
  results <- parSapply(cl, 1:10, function(x) simulate_voyage_opt(net = std_network_cal, initial = initial_std_network_cal,
                                                                 cs_input_cal = cs_input_top,  covid_input = covid_input,
                                                                 network_input = network_input, observed_data = observed_data,
                                                                 parameters_opt = parameters))
  stopCluster(cl)
  
  
  results_mean <- mean(results*10)
  
  return(results_mean)
}

initial_guess_C1_1 <- initial_guess[c(1,2,3,5)]
initial_guess_C1_1[3] <- 0.1

#run
initial <- Sys.time()
result_C1_1 <- optim( par = initial_guess_C1_1,  
                    fn = objective_function,
                    gr = NULL,
                    method = "Nelder-Mead", 
                    control = list(maxit = 100, abstol = 5, reltol = NULL,
                                   trace = TRUE))

final <- Sys.time()

save(result_C1_1, file = paste0(dir, "results/result_C1_1.rdata"))


