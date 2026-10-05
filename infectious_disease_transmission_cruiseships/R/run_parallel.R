library(parallel)

#inputs
#source functions
source("R/functions.R")

#load inputs
covid_parameters <- read_csv("Inputs/covid_parameters.csv") #disease-specific parameters
covid_input <- covid_parameters
scenarios <- read_csv("Inputs/scenarios_sample.csv",  col_types = readr::cols(`Testing Days` = readr::col_character(), .default = readr::col_guess())) #scenario settings to run
network_parameters <- read_csv("Inputs/network_parameters.csv") #network parameters including group sizes and min and max contacts for interactions 

load("Inputs/std_net.rdata") #contains all networks in the simulation, this script is ran with the standard ship network
load("Inputs/initial_std_net_I3.rdata") #identifies initial infected individuals, in this script this configuration is for a std cruise ship with 3 initial infections 

#set up

results = list() #save results
n_reps = 3 #number of replications per scenario
n_scenarios = nrow(scenarios) #number of scenarios to run

#run in parallel
num_cores <- 7
cl <-parallel::makeCluster(num_cores, setup_strategy = "sequential")


# Define or source your function on each node
clusterEvalQ(cl, {
  
  #libraries
  library(tidyr)
  library(dplyr)
  library(withr)
  library(parallel)
  library(readr)
  
  
  #source functions
  source("R/functions.R")
  
  #load inputs
  covid_parameters <- read_csv("Inputs/covid_parameters.csv") #disease-specific parameters
  covid_input <- covid_parameters
  scenarios <- read_csv("Inputs/scenarios_sample.csv",  col_types = readr::cols(`Testing Days` = readr::col_character(), .default = readr::col_guess())) #scenario settings to run
  network_parameters <- read_csv("Inputs/network_parameters.csv") #network parameters including group sizes and min and max contacts for interactions 
  
  load("Inputs/std_net.rdata") #contains all networks in the simulation, this script is ran with the standard ship network
  load("Inputs/initial_std_net_I3.rdata") #identifies initial infected individuals, in this script this configuration is for a std cruise ship with 3 initial infections 
  

})


# Run scenarios 
for (s in 1:n_scenarios) {
  
  type <-  scenarios[s,]$CruiseShip #get cruise ship type 
  initial <-  scenarios[s,]$Initial #get initial number of infections 
  index <- s #get index
  
  clusterExport(cl, varlist = c("s", "type", "initial", "index", "scenarios"))
  
  results[[index]] <- parSapply(cl, 1:n_reps,
                                function(x) simulate_voyage_onerep(net = get(type) , 
                                                                   initial = get(initial),
                                                                   cs_input = scenarios,
                                                                   scenario = s,
                                                                   covid_input = covid_parameters,
                                                                   network_input = network_parameters,
                                                                   observed_data = NULL,
                                                                   plotNetwork = FALSE,
                                                                   detailed_output = FALSE, stop_early = FALSE)) #change this as necessary 
  

  
  print(paste0("Scenario: ", s, "; ", "Cruise ship type: ", type, "; Initial set: ", initial))
  
}

stopCluster(cl)

#save results
results_sample_parallel <- results

save(results_sample_parallel, file = "Results/Simulation/results_sample_parallel.rdata")




