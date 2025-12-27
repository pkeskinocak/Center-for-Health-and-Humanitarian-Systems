#create demographics
library(readxl)
library(sna)

#functions
source("R/functions.R")

#1 Create ship demographic
#load group sizes and combinations per cruise ship type
std_group <- read_excel("Inputs/demographics_scenarios.xlsx", 
                         sheet = "standard")
fam_group <- read_excel("Inputs/demographics_scenarios.xlsx", 
                        sheet = "family")
senior_group <- read_excel("Inputs/demographics_scenarios.xlsx", 
                        sheet = "senior")


#load cruise ship configuration
cs_input <- read_excel("Inputs/cs_input_for_network_setup.xlsx", 
                      skip = 1)

#load info about networks
network_parameters <- read_csv("Inputs/network_parameters.csv") #network parameters including group sizes and min and max contacts for interactions 


#Create networks for the different type of ships
std_net <- makeDemographicNetwork(cs_input = cs_input, network_input = network_parameters, scenario = 1, dist = std_group)

fam_net <-  makeDemographicNetwork(cs_input = cs_input, network_input = network_parameters,scenario = 1, dist = fam_group)

sen_net <- makeDemographicNetwork(cs_input = cs_input, network_input = network_parameters,scenario = 1, dist = senior_group)
  

#check networks

#ages
table(std_net$passAge)/length(std_net$passAge)
table(fam_net$passAge)/length(fam_net$passAge)
table(sen_net$passAge)/length(sen_net$passAge)

#group sizes
table(table(std_net$passGroupID))/max(std_net$passGroupID)
table(table(fam_net$passGroupID))/max(fam_net$passGroupID)
table(table(sen_net$passGroupID))/max(sen_net$passGroupID)

#2 Create initial populations 

#load disease paramaters
covid_parameters <- read_csv("Inputs/covid_parameters.csv") #disease-specific parameters

#std - scenario 1
initial_std_1 <- simulation_initialize(net = std_net, cs_input = cs_input, scenario = 1, covid_input = covid_parameters, seed = 5, groupseed = 3,
                                       plotNetwork = TRUE)

#std - scenario 1
initial_std_2 <- simulation_initialize(net = std_net, cs_input = cs_input, scenario = 2, covid_input = covid_parameters, seed = 5, groupseed = 3,
                                       plotNetwork = TRUE)

#std - scenario 3
initial_std_3 <- simulation_initialize(net = std_net, cs_input = cs_input, scenario = 3, covid_input = covid_parameters, seed = 5, groupseed = 3,
                                       plotNetwork = TRUE)


#family - scenario 1
initial_fam_1 <- simulation_initialize(net = fam_net, cs_input = cs_input, scenario = 1, covid_input = covid_parameters, seed = 5, groupseed = 3,
                                       plotNetwork = TRUE)

#family - scenario 2
initial_fam_2 <- simulation_initialize(net = fam_net, cs_input = cs_input, scenario = 2, covid_input = covid_parameters, seed = 5, groupseed = 3,
                                       plotNetwork = TRUE)

#family - scenario 3
initial_fam_3 <- simulation_initialize(net = fam_net, cs_input = cs_input, scenario = 3, covid_input = covid_parameters, seed = 5, groupseed = 3,
                                       plotNetwork = TRUE)

#senior  scenario 1
initial_sen_1 <- simulation_initialize(net = sen_net, cs_input = cs_input, scenario = 1, covid_input = covid_parameters, seed = 5, groupseed = 3,
                                       plotNetwork = TRUE)
#senior  scenario 2
initial_sen_2 <- simulation_initialize(net = sen_net, cs_input = cs_input, scenario = 2, covid_input = covid_parameters, seed = 5, groupseed = 3,
                                       plotNetwork = TRUE)
#senior  scenario 3
initial_sen_3 <- simulation_initialize(net = sen_net, cs_input = cs_input, scenario = 3, covid_input = covid_parameters, seed = 5, groupseed = 3,
                                       plotNetwork = TRUE)
  

#save
save(std_net, fam_net, sen_net,
     initial_std_1, initial_std_2, initial_std_3,
     initial_fam_1, initial_fam_2, initial_fam_3,
     initial_sen_1, initial_sen_2, initial_sen_3,
     file = "Inputs/settings.rdata")

