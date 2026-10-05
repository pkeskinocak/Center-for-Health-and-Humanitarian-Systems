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



#run optimization
#with contraints
ui <- rbind(c(1,-1,0,0,0),
            c(1,0,-1,0,0),
            c(1,0,0,-1,0),
            c(0,-1,1,0,0),
            c(0,0,1,-1,0))

ci <- c(0,0,0,0,0)

initial_guess[3] <- initial_guess[3]  + 0.001


#run
result_A2 <- constrOptim(theta = initial_guess,  
                         f = function(parameters) objective_function(parameters), 
                         grad = NULL, ui = ui, ci = ci, 
                         outer.iterations = 100,
                         outer.eps = 5, 
                         method = "Nelder-Mead", 
                         control = list(maxit = 100, abstol = 5,reltol = NULL,
                                        trace = TRUE))


save(result_A2, file = paste0(dir, "results/result_A2.rdata"))