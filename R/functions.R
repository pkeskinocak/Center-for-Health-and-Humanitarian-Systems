#Agent-based simulation functions


## A. Agent functions 

### 1. makeAgent: This function creates an agent with a biostate (COVID state), age (we use this to create agents in their initial states)

makeAgent <- function(biostate,age=30)
  
{
  
  return (list(biostate=biostate,
               age=age,
               nextbiostate=NA,
               biostatecountdown=NA,
               viralload = NA,
               nextduration = NA))
}


### 2. setAgentState: This function gives the agent a specified state and assigns time duration in that state, and what is their next state. (we use this once agents transition to exposed and beyond)


setAgentState <- function(agent, new_biostate, covid_input) {
  
  # 1 Susceptible, 2 Exposed, 3 Asymptomatic, 4 Pre-symptomatic, 5 Symptomatic, 6 Quarantined (infected), 7 Quarantined (not infected) 
  # 8 Hospitalized, 9 Dead, 10 Recovered, 11 Removed 
  
  states = 11  # 11 state compartments 
  
  #keeps track of the agents original state (before transitioning)
  previous_biostate <- agent$biostate
  #new biostate to be assigned 
  biostate = new_biostate
  #set agent to the new biostate 
  agent$biostate <- biostate
  
  #The following section creates the durations for the new assigned biostate and the next biostate in the transition for the individual
  
  #If the new biostate assigned is ... 
  #(2) Exposed state 
  if (biostate == 2) {
    
    #Define exposed times
    bioMean <- dplyr::filter(covid_input, Parameter == "Te")$Mean
    bioSD  <- dplyr::filter(covid_input, Parameter == "Te")$SD
    
    #Define transition probabilities 
    Psym <- dplyr::filter(covid_input, Parameter == "Psym" & AgeGroup == agent$age)$Mean
    bioTransition <-c(0,0, (1-Psym), Psym, 0,0,0,0,0,0,0) 
    
    #Get agent next state
    agent$nextbiostate <- sample(1:states, prob=bioTransition,size=1)
    
    #Get agent duration in the current state
    agent$biostatecountdown <- min(ceiling(rlnorm(1, bioMean, bioSD)), 10) 
    
  }
  
  #(3) Infected asymptomatic state (gets a viral load)
  else if (biostate == 3) {
    
    #Define infected asymptomatic times
    bioMean <- dplyr::filter(covid_input, Parameter == "Tasym_rec")$Mean
    bioSD <-  dplyr::filter(covid_input, Parameter == "Tasym_rec")$SD
    
    #Get agent duration in current state
    agent$biostatecountdown <- min(ceiling(rlnorm(1, bioMean, bioSD)), 20) #other states, keep exposed duration
    
    #Create viral load curve, for every day in the asymptomatic state (lenght of viral load is always equal to the biostatecountdown assigned)
    agent$viralload <- getViralLoad(type = "a", asymp_recovery = agent$biostatecount)
    
    #Next state 
    agent$nextbiostate <- 10 #determine next state based on transition matrix (always to asym to recovered)
    
  }
  
  #(4) infected pre-symptomatic (gets a viral load)
  else if (biostate == 4) {
    
    #4 Get pre-symptomatic times 
    bioMean <-  dplyr::filter(covid_input, Parameter == "Tpre-sym")$Mean    
    bioSD <- dplyr::filter(covid_input, Parameter == "Tpre-sym")$SD
    
    #Get how long in the presymptomatic 
    countdown <- min(ceiling(rlnorm(1, bioMean, bioSD)), 10) #pre-symptomatic
    agent$biostatecountdown <- countdown
    
    #5 Get symptomatic times
    bioMean <-  dplyr::filter(covid_input, Parameter == "Tsym_rec")$Mean    
    bioSD <- dplyr::filter(covid_input, Parameter == "Tsym_rec")$SD
    
    #Get how long in symptomatic state
    agent$nextduration <- min(ceiling(rlnorm(1, bioMean, bioSD)), 20) #symptomatic
    
    #Create viral load curve, for every day in the pre-symptomatic and symptomatic states
    agent$viralload <- getViralLoad(type = "s", presymp_time = agent$biostatecountdown,
                                    symp_recovery = agent$nextduration)
    
    #Get agent next state 
    agent$nextbiostate <- 5 #determine next state based on transition matrix (always pre-symp to symp)
    
    
    
  }
  
  #(5) infected symptomatic (continues the viral load from pre-symptomatic)
  else if (biostate == 5) {
    
    #duration in the symptomatic state
    agent$biostatecountdown = agent$nextduration #gets the duration that we got in the pre-symptomatic step
    
    #transition probabilities
    Phosp <- dplyr::filter(covid_input, Parameter == "Phosp" & AgeGroup == agent$age)$Mean
    bioTransition <- c(0,0,0,0,0,0,0,Phosp, 0,1-Phosp,0)
    
    #viral load keeps the same from previous state (previous state is always pre-symptomatic)
    
    #next state
    agent$nextbiostate <- sample(1:states, prob=bioTransition,size=1) #determine next state based on transition matrix 
    agent$nextduration <- NA
  }
  
  #(8) hospitalized (viral load becomes zero because patient is no longer infecting others)
  else if (biostate == 8) {
    
    agent$viraload = NA
    
    #get transition probabilites 
    Pdead <- dplyr::filter(covid_input, Parameter == "Pdead" & AgeGroup == agent$age)$Mean
    bioTransition <- c(0,0,0,0,0,0,0,0, Pdead,1-Pdead,0)
    
    #next biostate
    agent$nextbiostate <- sample(1:states, prob=bioTransition,size=1) #determine next state based on transition matrix 
    
    #duration depends if dead or recovery
    if (agent$nextbiostate == 9 ) {  #if next state is (9) dead
      bioMean <-  dplyr::filter(covid_input, Parameter == "Thosp_dead")$Mean    
      bioSD <- dplyr::filter(covid_input, Parameter == "Thosp_dead")$SD
    } else { #if next state if (10) recovered 
      bioMean <-  dplyr::filter(covid_input, Parameter == "Thosp_rec")$Mean    
      bioSD <- dplyr::filter(covid_input, Parameter == "Thosp_rec")$SD
    }
    
    #get agent duration in current state
    agent$biostatecountdown <- min(ceiling(rlnorm(1, bioMean, bioSD)), 30)
    
  }
  
  #terminal states 9 dead, 10 recovered, 11 removed, or 1 susceptible 
  else {
    
    agent$biostatecountdown <- NA
    agent$nextbiostate <- NA  
    
  }
  
  
  
  return(agent) 
  
}


### 3. updateAgent: This function updates counter in the biostate countdown, if counter is 0, transition agent to next biostate 

updateAgent <- function(agent, covid_input)
  
{
  #print(i)
  if(!is.na(agent$biostatecountdown)) #proceed if agent has an active biostate countdown 
  {
    agent$biostatecountdown <- agent$biostatecountdown -1 #lower countdown by one 
    
    if(agent$biostatecountdown <= 0)  ##transition too next state
    {
      #agent <- transitionAgent(agent, covid_input)
      agent <- setAgentState(agent, agent$nextbiostate, covid_input)
      
    }
  }
  return(agent)
}

## B. Individual viral loads

### 1. getViralLoad:  Creates viral load for an infected agent. Gives ct values for infected period

#input = symptomatic or asymptomatic
# if asymptomatic:
# recovery rate (given)
# (draw profileration time within function)

# if symptomatic 
# pre-symptomatic duration (given)
# recovery smptomatic duration (given)
# draw time of of viral peak with respect to the symptom on set (0,1,2)

#returns viral load in Ct 

getViralLoad <- function(type = c("s", "a"),  presymp_time = NA, symp_recovery =NA,
                         asymp_recovery = NA){
  
  #when agent is symptomatic, parameters to be updated based on pathogen 
  if (type == "s") {
    
    #get parameters
    peak_vl_mean = 22.34788
    peak_vl_sd = 3.561828 
    peak_time = sample(c(0,1,2), size = 1)
    
    #create durations 
    tp = presymp_time + peak_time
    tc = max(presymp_time + symp_recovery - tp,1)
    
    #peak viral load
    peak_vl = rnorm(n =1 , mean = peak_vl_mean, sd = peak_vl_sd)
    
    
  } else {
    
    peak_vl_mean = 22.3798
    peak_vl_sd = 4.712329 
    
    #create duration
    tp = round(rlnorm(n = 1, meanlog = 0.9972, sdlog = 0.7625))
    tc = max(asymp_recovery - tp,1)
    
    #peak viral load
    peak_vl = rnorm(n =1 , mean = peak_vl_mean, sd = peak_vl_sd)
    
    
  }
  
  #create daily viral load 
  
  #first piece, proliferation length from 40 to peak
  if (tp > 1) {
    days = seq(1,tp,1)
  } else {
    days = 1
  }
  m = ifelse(length(days) != 1, (peak_vl - 40)/(tp-1), 1)
  piece1  = m*(days-1) + 40
  
  
  #second piece, clearance length from peak to 40
  days = seq(tp+1,tp+tc, 1)
  m = (40 - peak_vl)/tc
  piece2 = m*(days - tp) + peak_vl
  
  #combine
  viral_load = c(piece1,piece2)
  
  #return viral load
  return(viral_load)
  
  
}

### 2. convert_Ct_logGEML: Transforms viral load in Ct to RNA concentration 

convert_Ct_logGEML <- function(Ct, m_conv=-3.609714286, b_conv=40.93733333) {
  out <- (Ct-b_conv)/m_conv * log10(10) + log10(250)
  return(out) 
}


### 3. getTransmissibilityFactor_new: Convert viral load to transmissibility factor for the individual 

#take agents viral load,  agent characteristics
#return trans factor for the day in infection

getTransmissibilityFactor_new <- function(agent, OR = 1.3) {
  
  #get agent characteristics 
  viralload = agent$viralload
  biostate = agent$biostate
  
  #get viral load in correct units (Ct to viral load)
  trans_viralload = convert_Ct_logGEML(viralload)
  
  #calculate average 
  averageViralLoad = mean(trans_viralload)
  
  trans <- trans_viralload
  
  #get correct index
  if (biostate == 3 ) {         #asymptomatic 
    index = length(trans) - agent$biostatecountdown + 1
  } else if (biostate == 4 ) { #pre-symptomatic
    index = length(trans) - (agent$biostatecountdown + agent$nextduration) + 1
  } else {
    index = length(trans) - agent$biostatecountdown  + 1
  }
  
  viralLoad = trans[index]
  
  return(list(viralLoad = viralLoad,
              avg_viralLoad = averageViralLoad))
  
}

## C. Intervention functions 

### 1. getTestIDs_S: Identifies pool of symptomatic agents to test 

getTestIDs_S <- function(numberTest = 1000, pool_pass, pool_crew) {
  
  #empty vector
  tested_pass <- vector()
  tested_crew <- vector()
  
  #get agents states
  states_pass <- sapply(pool_pass,FUN=function(x){x$biostate}) #get states of all agents
  states_crew <- sapply(pool_crew,FUN=function(x){x$biostate}) #get states of all agents
  
  #get passenger ages
  ages_pass <- sapply(pool_pass,FUN=function(x){x$age}) 
  
  #Select agent IDs to be tested based on capacity
  
  #1 Symptomatic 65+
  pool <- intersect(which(states_pass == 5),  which(ages_pass == "65+")) 
  
  
  if (length(pool) == 0) {
    
    tested_pass <- tested_pass #empty vector
    tested_crew <- tested_crew #empty vector 
    
    #if more tests than people -> test everyone
  } else if (length(pool) <= numberTest) {
    
    tested_pass <- pool    #vector of IDS
    numberTest = numberTest - length(pool) #calculate leftover tests
    
    #if we are here (more tests than symptomatic 65+, we go to symptomatic all ages)
    
    if (numberTest == 0) {
      return(list(tested_pass = tested_pass, 
                  tested_crew = tested_crew,
                  numberTest = numberTest)) #stop here 
      
    }
    
  } else { #we have more symptomatic senior individuals than test, we draw and then stop. 
    
    tested_pass <- sample(pool, size = numberTest)
    numberTest = numberTest - length(tested_pass)
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest)) #stop here
    
  }
  
  #2 Symptomatic all ages
  constant <- length(states_pass)
  pool_symp_pass <- intersect(which(states_pass == 5),  which(ages_pass != "65+")) #pool of pass
  pool_symp_crew <- which(states_crew == 5) + constant  #pool of crew, we add to the crew IDs so they are the same as the passenger IDs
  pool <- c(pool_symp_pass, pool_symp_crew) #total pool
  
  if (length(pool) == 0) {
    
    tested_pass <- tested_pass
    tested_crew <- tested_crew
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest)) #stop here
    
    #if more tests than people -> test everyone
  } else if (length(pool) <= numberTest) {
    
    tested_pass <- c(tested_pass, pool_symp_pass)
    tested_crew <- pool_symp_crew - constant
    numberTest = numberTest - length(pool) #calculate leftover tests
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest))
    
    
  } else { #we have more symptomatic individuals than tests, we draw and then stop. 
    
    draw <- sample(pool, size = numberTest) #draw frrom symptomatic
    pass <- draw[which(draw <= constant)] #get pass
    crew <- draw[which(draw > constant)] - constant #get crew
    
    tested_pass <- c(tested_pass, pass) #update tested IDS
    tested_crew <- crew  #update tested IDS
    
    numberTest = numberTest - length(draw)
    
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest))
    
    
  }
}

### 2. getTestIDs_S_CC: Identifies pool of symptomatic + close contacts agents to test 

getTestIDs_S_CC <- function(numberTest = 1000, pool_pass, pool_crew, net) {
  
  #empty vector
  tested_pass <- vector()
  tested_crew <- vector()
  
  #get agents states
  states_pass <- sapply(pool_pass,FUN=function(x){x$biostate}) #get states of all agents
  states_crew <- sapply(pool_crew,FUN=function(x){x$biostate}) #get states of all agents
  
  #get passenger ages
  ages_pass <- sapply(pool_pass,FUN=function(x){x$age}) 
  
  #Select agent IDs to be tested based on capacity
  
  #1 Symptomatic 65+
  pool <- intersect(which(states_pass == 5),  which(ages_pass == "65+")) 
  
  
  #if more tests than people -> test everyone
  if (length(pool) == 0) {
    
    tested_pass <- tested_pass
    tested_crew <- tested_crew
    
  } else if (length(pool) <= numberTest) {
    
    tested_pass <- pool    #vector of IDS
    
    numberTest = numberTest - length(pool) #calculate leftover tests
    
    if (numberTest == 0) {
      return(list(tested_pass = tested_pass, 
                  tested_crew = tested_crew,
                  numberTest = numberTest))
      
    }
    
  } else { #we have more symptomatic senior individuals than test, we draw and then stop. 
    
    tested_pass <- sample(pool, size = numberTest)
    numberTest = numberTest - length(tested_pass)
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest))
    
  }
  
  #2 Symptomatic all ages
  constant <- length(states_pass)
  pool_symp_pass <- intersect(which(states_pass == 5),  which(ages_pass != "65+")) #pool of pass
  pool_symp_crew <- which(states_crew == 5) + constant  #pool of crew, we add to the crew IDs so they are the same as the passenger IDs
  pool <- c(pool_symp_pass, pool_symp_crew) #total pool
  
  if (length(pool) == 0) {
    
    tested_pass <- tested_pass
    tested_crew <- tested_crew
    
    #if more tests than people -> test everyone
    
  } else if (length(pool) <= numberTest) {
    
    tested_pass <- c(tested_pass, pool_symp_pass)
    tested_crew <- pool_symp_crew - constant
    numberTest = numberTest - length(pool) #calculate leftover tests
    
    if (numberTest == 0) {
      return(list(tested_pass = tested_pass, 
                  tested_crew = tested_crew,
                  numberTest = numberTest))
    }
    
  } else { #we have more symptomatic individuals than tests, we draw and then stop. 
    
    draw <- sample(pool, size = numberTest) #draw frrom symptomatic
    pass <- draw[which(draw <= constant)] #get pass
    crew <- draw[which(draw > constant)] - constant #get crew
    
    tested_pass <- c(tested_pass, pass) #update tested IDS
    tested_crew <- crew  #update tested IDS
    
    numberTest = numberTest - length(draw)
    
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest))
    
    
  }
  
  # 3   Close contacts from travel main group and cabin main group 
  
  pool_pass_symptomatic <- which(states_pass == 5 )
  pool_crew_symptomatic <- which(states_crew == 5 )
  
  #get close contacts
  
  # PASSENGERS
  #pass groups
  cc_groups_pass <- unique(net$passGroupID[pool_pass_symptomatic])
  #passengers in groups
  cc_all_pass <-  which(net$passGroupID %in% cc_groups_pass) #all passenger in selected groups
  cc_pass <- setdiff(cc_all_pass,pool_pass_symptomatic) #get rid of those already tested
  
  # CREW 
  #crew groups
  cc_groups_crew <- unique(net$crewGroupID[pool_crew_symptomatic])
  #crew in groups
  cc_all_crew <-  which(net$crewGroupID %in% cc_groups_crew) #all crew in selected groups
  cc_crew <- setdiff(cc_all_crew,pool_crew_symptomatic) #get rid of those already tested
  
  # First passengers 
  if (length(cc_pass)  == 0 ) {
    
    tested_pass <- tested_pass
    
    #if more tests than people -> test everyone
  } else if (length(cc_pass) <= numberTest) {
    
    tested_pass <- c(tested_pass, cc_pass)    #vector of IDS
    numberTest = numberTest - length(cc_pass) #calculate leftover tests
    
    if (numberTest == 0) {
      return(list(tested_pass = tested_pass, 
                  tested_crew = tested_crew,
                  numberTest = numberTest))
    }
    
  } else { #we have more more close contacts than tests, we draw and then stop. 
    
    draw <- sample(cc_pass, size = numberTest)
    tested_pass <- c(tested_pass, draw)
    
    numberTest = numberTest - length(draw)
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest))
    
  }
  
  #continue with crew if test are remaining 
  # First passengers 
  if (length(cc_crew)  == 0 ) {
    
    tested_crew <- tested_crew
    
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest))
    
    #if more tests than people -> test everyone
  } else if (length(cc_crew) <= numberTest) {
    
    tested_crew <- c(tested_crew, cc_crew)    #vector of IDS
    numberTest = numberTest - length(cc_crew) #calculate leftover tests
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest))
    
    
  } else { #we have more more close contacts than tests, we draw and then stop. 
    
    draw <- sample(cc_crew, size = numberTest)
    tested_crew <- c(tested_pass, draw)
    
    numberTest = numberTest - length(draw)
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest))
    
  }
  
  
}

### 3. getTestIDs_asym: Identifies pool of asymptomatic agents to test

getTestIDs_asym <- function(numberTest = 1000, pool_pass, pool_crew, tested_pass_prev,
                            tested_crew_prev) {
  
  #empty vector
  tested_pass <- vector()
  tested_crew <- vector()
  
  #get agents states
  states_pass <- sapply(pool_pass,FUN=function(x){x$biostate}) #get states of all agents
  states_crew <- sapply(pool_crew,FUN=function(x){x$biostate}) #get states of all agents
  
  #get passenger ages
  ages_pass <- sapply(pool_pass,FUN=function(x){x$age}) 
  
  #No Symptomatic all ages
  
  #Passengers
  pool_pass <- which(states_pass %in% c(1,2,3,4)) #all ages
  #get rid of the ones previously tested
  pool_pass <- setdiff(pool_pass, tested_pass_prev)
  
  #Crew
  constant <- length(states_pass)
  pool_crew <- which(states_crew %in% c(1,2,3,4))
  #get rid of the ones previously tested
  pool_crew <- setdiff(pool_crew, tested_crew_prev)
  pool_crew <- pool_crew + constant  #pool of crew, we add to the crew IDs so they are the same as the passenger IDs
  
  
  pool <- c(pool_pass, pool_crew) #total pool
  
  #no tests
  if (length(pool) ==  0){
    tested_pass <- tested_pass
    tested_crew <- tested_crew
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew))
    
    #if more tests than people -> test everyone
  } else if (length(pool) <= numberTest) {
    
    tested_pass <- c(tested_pass, pool_pass)
    tested_crew <- c(tested_crew, pool_crew - constant)
    numberTest = numberTest - length(pool) #calculate leftover tests
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew))
    
    # we have more people than tests, we draw   
  } else {  
    
    draw <- sample(pool, size = numberTest) #draw froom no symptoms pool 
    pass <- draw[which(draw <= constant)] #get pass
    crew <- draw[which(draw > constant)] - constant #get crew
    
    tested_pass <- c(tested_pass, pass) #update tested IDS
    tested_crew <- c(tested_crew, crew)  #update tested IDS
    
    numberTest = numberTest - length(draw)
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew))
    
  }
}



### 4. getTestIDs_S_CC_D: Identifies pool of symptomatic agents to test with delay 

getTestIDs_S_CC_D <- function(numberTest = 1000, pool_pass, pool_crew, net) {
  
  #empty vector
  tested_pass <- vector()
  tested_crew <- vector()
  
  #get agents states
  states_pass <- sapply(pool_pass,FUN=function(x){x$biostate}) #get states of all agents
  states_crew <- sapply(pool_crew,FUN=function(x){x$biostate}) #get states of all agents
  
  #get passenger ages
  ages_pass <- sapply(pool_pass,FUN=function(x){x$age}) 
  
  #Select agent IDs to be tested based on capacity
  
  #1 Symptomatic 65+
  pool <- intersect(which(states_pass == 5),  which(ages_pass == "65+")) 
  
  
  #if more tests than people -> test everyone
  if (length(pool) == 0) {
    
    tested_pass <- tested_pass
    tested_crew <- tested_crew
    
  } else if (length(pool) <= numberTest) {
    
    tested_pass <- pool    #vector of IDS
    numberTest = numberTest - length(pool) #calculate leftover tests
    
    if (numberTest == 0) {
      return(list(tested_pass = tested_pass, 
                  tested_crew = tested_crew,
                  numberTest = numberTest))
      
    }
    
  } else { #we have more symptomatic senior individuals than test, we draw and then stop. 
    
    tested_pass <- sample(pool, size = numberTest)
    numberTest = numberTest - length(tested_pass)
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest))
    
  }
  
  #2 Symptomatic all ages
  constant <- length(states_pass)
  pool_symp_pass <- intersect(which(states_pass == 5),  which(ages_pass != "65+")) #pool of pass
  pool_symp_crew <- which(states_crew == 5) + constant  #pool of crew, we add to the crew IDs so they are the same as the passenger IDs
  pool <- c(pool_symp_pass, pool_symp_crew) #total pool
  
  if (length(pool) == 0) {
    
    tested_pass <- tested_pass
    tested_crew <- tested_crew
    
    #if more tests than people -> test everyone
    
  } else if (length(pool) <= numberTest) {
    
    tested_pass <- c(tested_pass, pool_symp_pass)
    tested_crew <- pool_symp_crew - constant
    numberTest = numberTest - length(pool) #calculate leftover tests
    
    if (numberTest == 0) {
      return(list(tested_pass = tested_pass, 
                  tested_crew = tested_crew,
                  numberTest = numberTest))
    }
    
  } else { #we have more symptomatic individuals than tests, we draw and then stop. 
    
    draw <- sample(pool, size = numberTest) #draw frrom symptomatic
    pass <- draw[which(draw <= constant)] #get pass
    crew <- draw[which(draw > constant)] - constant #get crew
    
    tested_pass <- c(tested_pass, pass) #update tested IDS
    tested_crew <- crew  #update tested IDS
    
    numberTest = numberTest - length(draw)
    
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                numberTest = numberTest))
    
    
  }
  
  # 3   Close contacts from travel main group and cabin main group 
  
  pool_pass_symptomatic <- which(states_pass == 5 )
  pool_crew_symptomatic <- which(states_crew == 5 )
  
  #get close contacts
  
  # PASSENGERS
  #pass groups
  cc_groups_pass <- unique(net$passGroupID[pool_pass_symptomatic])
  #passengers in groups
  cc_all_pass <-  which(net$passGroupID %in% cc_groups_pass) #all passenger in selected groups
  cc_pass <- setdiff(cc_all_pass,pool_pass_symptomatic) #get rid of those already tested
  
  # CREW 
  #crew groups
  cc_groups_crew <- unique(net$crewGroupID[pool_crew_symptomatic])
  #crew in groups
  cc_all_crew <-  which(net$crewGroupID %in% cc_groups_crew) #all crew in selected groups
  cc_crew <- setdiff(cc_all_crew,pool_crew_symptomatic) #get rid of those already tested
  
  # First passengers 
  if (length(cc_pass)  == 0 ) {
    
    nextday_pass <- NULL
    
    #if more tests than people -> test everyone
  } else if (length(cc_pass) <= numberTest) {
    
    nextday_pass <- cc_pass
    numberTest = numberTest - length(cc_pass) #calculate leftover tests
    
  } else { #we have more more close contacts than tests, we draw and then stop. 
    
    draw <- sample(cc_pass, size = numberTest)
    nextday_pass <- draw
    numberTest = numberTest - length(draw)
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                nextday_pass = nextday_pass,
                numberTest = numberTest))
    
  }
  
  #continue with crew if test are remaining 
  # First passengers 
  if (length(cc_crew)  == 0 ) {
    
    nextday_crew <- NULL
    
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                nextday_pass = nextday_pass,
                numberTest = numberTest))
    
    #if more tests than people -> test everyone
  } else if (length(cc_crew) <= numberTest) {
    
    nextday_crew <- cc_crew   #vector of IDS
    numberTest = numberTest - length(cc_crew) #calculate leftover tests
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                nextday_pass = nextday_pass,
                nextday_crew = nextday_crew,
                numberTest = numberTest))
    
    
  } else { #we have more more close contacts than tests, we draw and then stop. 
    
    draw <- sample(cc_crew, size = numberTest)
    nextday_crew- draw
    
    numberTest = numberTest - length(draw)
    
    return(list(tested_pass = tested_pass, 
                tested_crew = tested_crew,
                nextday_pass = nextday_pass,
                nextday_crew = nextday_crew,
                numberTest = numberTest))
    
  }
  
  
}

### 5. testDaily: Perform testing and identifies agents that go to isolation/quarantine 

#Testing procedures -> Tested individuals go to the quarantine states
# This function takes: 
# tested_pass and tested_crew obtained from getTestIDs , it has the crew and pass IDs
# true positive and false positive error
# Process: Sends tested people to quarantine infected and quarantine not infected
# Outputs: Pool of passengers + crew. Table with the states of the agents that went to quarantine
# pool_pass, pool_crew, state_pass_quarantined, state_crew_quarantined


testDaily <- function(tested_pass, tested_crew, pool_pass, pool_crew, fn_error_sym = 0.01, fp_error_sym = 0.01,
                      fn_error_asym = 0.01, fp_error_asym = 0.01) {
  
  #get states for passengers and crew 
  states_pass <- sapply(pool_pass,FUN=function(x){x$biostate}) #get states of all agents
  states_crew <- sapply(pool_crew,FUN=function(x){x$biostate}) #get states of all agents
  
  #(1a) Infected detectable - No symptoms: (3) Asymptomatic, (4) Pre-symptomatic 
  
  #Passenger
  pass_inf_detectable_all_nosym <- intersect(tested_pass, which(states_pass %in% c(3,4))) #get passenger IDs
  
  if (length(pass_inf_detectable_all_nosym) > 1) {
    pass_inf_detectable_nosym <- sample(pass_inf_detectable_all_nosym, size = round(length(pass_inf_detectable_all_nosym)*(1-fn_error_asym))) #sample to account for sampling error 
  } else {
    pass_inf_detectable_nosym <- pass_inf_detectable_all_nosym
  }
  
  #Crew 
  crew_inf_detectable_all_nosym <- intersect(tested_crew, which(states_crew %in% c(3,4))) #get crew IDs
  
  if (length(crew_inf_detectable_all_nosym) > 1) {
    crew_inf_detectable_nosym <- sample(crew_inf_detectable_all_nosym, size = round(length(crew_inf_detectable_all_nosym)*(1-fn_error_asym))) #sample to account for sampling error 
  } else {
    crew_inf_detectable_nosym <- crew_inf_detectable_all_nosym
  }
  
  #(1b) Infected detectable - Symptoms: (5) Symptomatic
  
  #Passenger
  pass_inf_detectable_all_sym <- intersect(tested_pass, which(states_pass %in% c(5))) #get passenger IDs
  
  if (length(pass_inf_detectable_all_sym) > 1) {
    pass_inf_detectable_sym <- sample(pass_inf_detectable_all_sym, size = round(length(pass_inf_detectable_all_sym)*(1-fn_error_sym))) #sample to account for sampling error 
  } else {
    pass_inf_detectable_sym <- pass_inf_detectable_all_sym
  }
  
  #Crew
  crew_inf_detectable_all_sym <- intersect(tested_crew, which(states_crew %in% c(5))) #get crew IDs
  
  if (length(crew_inf_detectable_all_sym) > 1) {
    crew_inf_detectable_sym <- sample(crew_inf_detectable_all_sym, size = round(length(crew_inf_detectable_all_sym)*(1-fn_error_sym))) #sample to account for sampling error 
  } else {
    crew_inf_detectable_sym <-  crew_inf_detectable_all_sym
  }
  
  #combine 1a and 1b
  pass_inf_detectable <- c(pass_inf_detectable_nosym,  pass_inf_detectable_sym)
  crew_inf_detectable <- c(crew_inf_detectable_nosym,  crew_inf_detectable_sym)
  
  #(2) Infected not detectable: (2) Exposed 
  
  #Passenger
  pass_inf_notdetectable_all <- intersect(tested_pass, which(states_pass == 2)) #get passenger IDs
  
  if (length(pass_inf_notdetectable_all) > 1) {
    pass_inf_notdetectable <- sample(pass_inf_notdetectable_all, size = round(length(pass_inf_notdetectable_all)*fp_error_asym)) #sample to account for sampling error 
  } else {
    pass_inf_notdetectable <-  pass_inf_notdetectable_all
  }
  
  #Crew
  crew_inf_notdetectable_all <- intersect(tested_crew, which(states_crew == 2)) #get crew IDs
  
  if (length(crew_inf_notdetectable_all) > 1) {
    crew_inf_notdetectable <- sample(crew_inf_notdetectable_all, size = round(length(crew_inf_notdetectable_all)*fp_error_asym)) #sample to account for sampling error 
  } else {
    crew_inf_notdetectable <-  crew_inf_notdetectable_all
  }
  
  #(3) Not infected: (1) Susceptible (assume always asymptomatic)
  
  #Passenger
  pass_notinf_all <- intersect(tested_pass, which(states_pass == 1)) #get passenger IDs
  
  if (length(pass_notinf_all) > 1) {
    pass_notinf <- sample(pass_notinf_all, size = round(length(pass_notinf_all)*fp_error_asym)) #sample to account for sampling error 
  } else {
    pass_notinf <- pass_notinf_all
  }
  
  #Crew 
  crew_notinf_all <- intersect(tested_crew, which(states_crew == 1)) #get crew IDs
  
  if (length(crew_notinf_all) >1 ) {
    crew_notinf <- sample(crew_notinf_all, size = round(length(crew_notinf_all)*fp_error_asym)) #sample to account for sampling error 
  } else {
    crew_notinf <- crew_notinf_all
  }
  
  #change status
  #passengers 
  #(1) quarantine infected
  for (i in c(pass_inf_detectable, pass_inf_notdetectable)) {
    pool_pass[[i]] <- setAgentState(pool_pass[[i]],6, covid_input) 
  }
  
  #(2) quarantine not infected
  for (i in pass_notinf) {
    pool_pass[[i]] <- setAgentState(pool_pass[[i]],7, covid_input) 
  }
  
  #crew
  #(1) quarantine infected
  for (i in c(crew_inf_detectable, crew_inf_notdetectable)) {
    pool_crew[[i]] <- setAgentState(pool_crew[[i]],6, covid_input) 
  }
  
  #(2) quarantine not infected
  for (i in crew_notinf) {
    pool_crew[[i]] <- setAgentState(pool_crew[[i]],7, covid_input) 
  }
  
  return(list(pool_pass = pool_pass, 
              pool_crew = pool_crew))
}


### 6. removePassenger: Remove passengers to remove from active pool

#takes pool of passengers  already present in the function
#return pool of passengers with updates
removePassengers <- function(numberRemoved = 100, pool_pass) {
  
  #get pax states
  states_pass <- sapply(pool_pass,FUN=function(x){x$biostate}) #get states of all agents
  
  #get agents belonging in the susceptible, exposed, asymptomatic, pre-symptomatic states
  eligible_pass <- which(states_pass %in% c(1,2,3,4))
  
  if (length(eligible_pass) > numberRemoved) {
    removedIDs <- sample(eligible_pass, size = numberRemoved)
  } else {
    removedIDs <- eligible_pass
  }
  
  #update state of those removed
  for (i in removedIDs) {
    pool_pass[[i]] <- setAgentState(pool_pass[[i]],11, covid_input) 
  }
  
  return(pool_pass = pool_pass)
  
  
}


### 7. masking_selection: Selection of agents for masking intervention

masking_selection <- function(network, selection = c("random", "structured"), coverage) {
  
  #cs_input <- cs_input_scenario %>% filter(ID == "5-A-1")
  
  #random selection: select x number of agents randomly
  if (selection == "random") {
    
    #passengers
    numPass <- length(network$passAge)
    numSelection <- round(coverage[1]/100*numPass)
    pass_selection <- sample(1:numPass, size = numSelection)
    
    #crew 
    numCrew <- length(network$crewAge)
    numSelection <- round(coverage[2]/100*numCrew)
    crew_selection_emp  <- sample(1:numCrew, size = numSelection)
    
  } else if (selection == "structured") {
    
    ##PASSENGERS 
    
    #seniors more likely to be masked
    #pick senior in one group: within that group, prob that senior is masked 80%, adult 70%, kid 60%
    #perform steps until mask coverage is reached
    
    #coverage
    numPass <- length(network$passAge)
    numSelection <- round(coverage[1]/100*numPass)
    
    #thresholds
    thresholds <- c("0-12" = 0.5, "13-17" = 0.6, "18-64" = 0.7, "65+" = 0.8)
    
    #create selection placeholder
    pass_selection <- vector()
    seniorChosen <- vector()
    
    while (length(pass_selection) < numSelection) {
      
      #target number
      numPass <- length(network$passAge)
      passTarget <- round(coverage[1]/100*numPass)
      
      #passenger ages
      passAges <- network$passAge
      
      #pick senior 
      seniorChoices <- setdiff(which(passAges == "65+"), seniorChosen)
      
      if (length(seniorChoices) > 1) {
        seniorPrelim <- sample(seniorChoices, size = 1)
      } else {
        choices <- setdiff(1:numPass,  pass_selection)
        seniorPrelim <- sample(choices, size = 1)  #once there's no senior to choose from, select randomly from eligible 
      }
      
      #get group
      passGroups <- network$passGroupID
      passPrelim <- which(passGroups == passGroups[seniorPrelim])
      agesPrelim <- passAges[c(passPrelim)]
      draw <- runif(length(agesPrelim))
      passSelected <- passPrelim[which((draw < thresholds[agesPrelim]) == TRUE)]
      
      #update list of seniors called
      seniorChosen <- unique(c(seniorChosen, seniorPrelim, passPrelim[which(agesPrelim == "65+")]))
      
      #add to running selection
      pass_selection <- unique(c(pass_selection, passSelected))
    }
    
    #crew 80% 
    numCrew <- length(network$crewAge)
    numSelection <- round(coverage[2]/100*numCrew)
    crew_selection_emp  <- sample(1:numCrew, size = numSelection)
    
    
  }
  
  pass_selection <- sort(pass_selection)
  #crew_selection_main <- sort(crew_selection_main)
  crew_selection_emp <- sort(crew_selection_emp)
  
  return(list(pass_selection, crew_selection_emp))
}


### 8. vaccination_selection_implemention: Selection of agents for vaccination internvention and assigns them to recovered stage depending on effectiveness

vaccination_selection_implemention <- function(network, pool_pass, pool_crew, selection = c("random", "structured"), coverage,
                                               effectiveness) {
  
  coverage_kids = coverage[1]
  coverage_teens = coverage[2]
  coverage_adults = coverage[3]
  coverage_seniors = coverage[4]
  
  eff_kids = effectiveness[1]
  eff_teens = effectiveness[2]
  eff_adults = effectiveness[3]
  eff_seniors = effectiveness[4]
  
  
  #random selection: select x number of agents randomly
  if (selection == "random") {
    
    # Create variables
    pass_kids_selection <- vector()
    pass_teens_selection <- vector()
    pass_adults_selection <- vector()
    pass_seniors_selection <- vector()
    
    #passengers kids
    numPass_kids <- length(which(network$passAge == "0-12")) 
    
    if (numPass_kids > 0 ) {
      numSelection <- round(coverage_kids/100*numPass_kids)
      pass_kids_selection <- sample(which(network$passAge %in% c("0-12")), size = numSelection)
    }
    
    #passenger teens
    
    numPass_teens <- length(which(network$passAge == "13-17")) 
    
    if (numPass_teens  > 0 ) {
      numSelection <- round(coverage_teens/100*numPass_teens)
      pass_teens_selection <- sample(which(network$passAge %in% c("13-17")), size = numSelection)
    }
    
    #passenger adults
    numPass_adults <- length(which(network$passAge == "18-64")) 
    numSelection <- round(coverage_adults/100*numPass_adults)
    pass_adults_selection <- sample(which(network$passAge == "18-64"), size = numSelection)
    
    #passenger senior
    numPass_seniors <- length(which(network$passAge == "65+")) 
    numSelection <- round(coverage_seniors/100*numPass_seniors)
    pass_seniors_selection <- sample(which(network$passAge == "65+"), size = numSelection)
    
    #crew 
    numCrew <- length(network$crewAge) 
    numSelection <- round(90/100*numCrew)
    crew_selection <- sample(1:numCrew, size = numSelection)
    
  } else if (selection == "structured") {
    
    ##PASSENGERS 
    
    #seniors more likely to be vaccinate
    #pick senior in one group: within that group, prob that senior is vaccinated 90%, adult 80%, kid 70%
    #pick adult in one group: within that group, prob that senior is vaccinated 90%, adult 80%, kid 70%
    #intercalating senior and adult ensures that group with no senior are vaccinated as well. 
    #perform steps until vaccination coverage is reached
    
    # (1) Create thresholds
    #thresholds <- c("0-12" = 0.95, "13-17" = 0.95, "18-64" = 0.8, "65+" = 0.9)
    #thresholds <- c("0-12" = 0, "13-17" = 0, "18-64" = 0, "65+" = 0)
    
    #if adults in the group are vaccinated, pick kids until threshold is met 
    
    passAges <- network$passAge
    
    # (2) Create variables
    pass_kids_selection <- vector()
    pass_teens_selection <- vector()
    pass_adults_selection <- vector()
    pass_seniors_selection <- vector()
    
    # (3) Coverage by age
    numPass_kids <- length(which(passAges== "0-12")) 
    numSelection_kids <- round(coverage_kids/100*numPass_kids)
    
    numPass_teens <- length(which(passAges== "13-17")) 
    numSelection_teens <- round(coverage_teens/100*numPass_teens)
    
    numPass_adults <- length(which(passAges == "18-64")) 
    numSelection_adults <- round(coverage_adults/100*numPass_adults)
    
    numPass_seniors <- length(which(passAges == "65+")) 
    numSelection_seniors <- round(coverage_seniors/100*numPass_seniors)
    
    # (4) Start kids/teens selection since this is the group which is most sparse
    
    #group IDs
    passGroups <- network$passGroupID
    
    # Teens
    if (numSelection_teens  > 0 ) {
      
      #select sample 
      pass_teens_selection <- sample(which(network$passAge %in% c("13-17")), size = numSelection_teens)
      
      # Get adults in their groups
      passPrelim <- which(passGroups %in% passGroups[pass_teens_selection])
      agesPrelim <- passAges[passPrelim]
      
      #Choose adults in their group and other kids as well 
      thresholds_teens =  c("0-12" = 0.9, "13-17" = 1, "18-64" = 0.95, "65+" = 0.95)
      
      # Draw passengers from the same groups of selected teens
      df <- data.frame("ID" = passPrelim, "ages" = agesPrelim)
      
      #random draw
      draw <- runif(length(agesPrelim))
      passSelected <- passPrelim[which((draw < thresholds_teens[agesPrelim]) == TRUE)]
      df <- df %>%
        filter(ID %in% passSelected)
      
      #seniors 
      pass_seniors_selection <- c(pass_seniors_selection, 
                                  df %>% filter(ages == "65+") %>% pull(ID))
      
      #adults 
      pass_adults_selection <- c(pass_adults_selection, 
                                 df %>% filter(ages == "18-64") %>% pull(ID))
      
      #kids
      pass_kids_selection <- c(pass_kids_selection, 
                               df %>% filter(ages == "0-12") %>% pull(ID))
      
    }
    
    # Continue with kids 
    
    if (numSelection_kids  > 0 & length(pass_kids_selection) < numSelection_kids) {
      
      remaining <- numSelection_kids - length(pass_kids_selection)
      eligible <- setdiff(which(network$passAge %in% c("0-12")), pass_kids_selection)
      pass_selection <- sample(eligible, size = remaining)
      
      # Get adults in their groups
      passPrelim <- which(passGroups %in% passGroups[pass_selection])
      agesPrelim <- passAges[passPrelim]
      
      thresholds_kids =  c("0-12" = 1, "13-17" = 0, "18-64" = 0.95, "65+" = 0.95)
      
      # Draw passengers from the same groups of selected teens
      df <- data.frame("ID" = passPrelim, "ages" = agesPrelim)
      
      #random draw
      draw <- runif(length(agesPrelim))
      passSelected <- passPrelim[which((draw < thresholds_kids[agesPrelim]) == TRUE)]
      df <- df %>%
        filter(ID %in% passSelected)
      
      #seniors 
      pass_seniors_selection <- c(pass_seniors_selection, 
                                  df %>% filter(ages == "65+") %>% pull(ID))
      pass_seniors_selection <- unique(pass_seniors_selection)
      
      #adults 
      pass_adults_selection <- c(pass_adults_selection, 
                                 df %>% filter(ages == "18-64") %>% pull(ID))
      pass_adults_selection <- unique(pass_adults_selection )
      #kids
      pass_kids_selection <- c(pass_kids_selection, pass_selection)
      pass_kids_selection  <- unique(pass_kids_selection)
    }
    
    
    # (4) Selection of adults until coverage is reached 
    
    seniorsChosen <- vector()
    adultsChosen <- vector()
    
    passGroups <- network$passGroupID
    
    
    
    while (length(pass_seniors_selection) < numSelection_seniors ||
           length(pass_adults_selection) < numSelection_adults ||
           length(pass_kids_selection) < numSelection_kids   ||
           length(pass_teens_selection) < numSelection_teens ) {
      
      thresholds <- c("0-12" = 0, "13-17" = 0, "18-64" = 0.8, "65+" = 0.9)
      
      remaining <- data.frame("ages" = c("65+", "18-64","13-17", "0-12"),
                              "remaining" = c(numSelection_seniors - length(pass_seniors_selection),
                                              numSelection_adults - length(pass_adults_selection),
                                              numSelection_teens - length(pass_teens_selection),
                                              numSelection_kids - length(pass_kids_selection)))
      
      
      
      
      # 4.1. Pick senior 
      if (remaining[1,2] > 0) {
        
        seniorChoices <- setdiff(which(passAges == "65+"), pass_seniors_selection)
        
        if (length(seniorChoices) > 1) {
          seniorPrelim <- sample(seniorChoices, size = 1)
        } else {
          seniorPrelim <- sample(which(passAges == "65+"), size = 1)
        }
        
        # 4.1.1 Pick others within group 
        
        #get group
        passPrelim <- which(passGroups == passGroups[seniorPrelim])
        agesPrelim <- passAges[c(passPrelim)]
        df <- data.frame("ID" = passPrelim, "ages" = agesPrelim)
        
        #only keep 
        #random draw
        draw <- runif(length(agesPrelim))
        passSelected <- passPrelim[which((draw < thresholds[agesPrelim]) == TRUE)]
        df <- df %>%
          filter(ID %in% passSelected)
        
        
        #check if reached capacity or not
        df <- df %>%
          left_join(remaining, by = "ages") %>%
          mutate(select = if_else(remaining > 0, "selected", "not selected")) %>%
          filter(select == "selected")
        
        #put items in corresponding vector
        
        #seniors 
        pass_seniors_selection <- c(pass_seniors_selection, 
                                    df %>% filter(ages == "65+") %>% pull(ID))
        pass_seniors_selection <- unique(pass_seniors_selection)
        
        #adults 
        pass_adults_selection <- c(pass_adults_selection, 
                                   df %>% filter(ages == "18-64") %>% pull(ID))
        pass_adults_selection <- unique(  pass_adults_selection)
        
        
      }
      
      
      # 5.1. Pick adults
      
      remaining <- data.frame("ages" = c("65+", "18-64","13-17", "0-12"),
                              "remaining" = c(numSelection_seniors - length(pass_seniors_selection),
                                              numSelection_adults - length(pass_adults_selection),
                                              numSelection_teens - length(pass_teens_selection),
                                              numSelection_kids - length(pass_kids_selection)))
      
      
      if (remaining[2,2] > 0) {
        
        adultsChoices <- setdiff(which(passAges == "18-64"), pass_adults_selection)
        
        if (length(adultsChoices) > 1) {
          adultPrelim <- sample(adultsChoices, size = 1)
        } else {
          adultPrelim <- sample(which(passAges == "18-64"), size = 1)
        }
        
        # 4.1.1 Pick others within group 
        
        #get group
        passPrelim <- which(passGroups == passGroups[adultPrelim])
        agesPrelim <- passAges[c(passPrelim)]
        df <- data.frame("ID" = passPrelim, "ages" = agesPrelim)
        
        #random draw
        draw <- runif(length(agesPrelim))
        passSelected <- passPrelim[which((draw < thresholds[agesPrelim]) == TRUE)]
        df <- df %>%
          filter(ID %in% passSelected)
        
        #check if reached capacity or not
        df <- df %>%
          left_join(remaining, by = "ages") %>%
          mutate(select = if_else(remaining > 0, "selected", "not selected")) %>%
          filter(select == "selected")
        
        #put items in corresponding vector
        
        #seniors 
        pass_seniors_selection <- c(pass_seniors_selection, 
                                    df %>% filter(ages == "65+") %>% pull(ID))
        pass_seniors_selection <- unique(pass_seniors_selection) 
        #adults 
        pass_adults_selection <- c(pass_adults_selection, 
                                   df %>% filter(ages == "18-64") %>% pull(ID))
        pass_adults_selection <- unique( pass_adults_selection)
        
      }
      
      
      
      
    }
    
    #crew 
    numCrew <- length(network$crewAge) 
    numSelection <- round(90/100*numCrew)
    crew_selection <- sample(1:numCrew, size = numSelection)
    
  }
  
  #Implementation part
  
  #Pass Kids selected
  for (i in pass_kids_selection) {
    if (!pool_pass[[i]]$biostate %in% c(3,4,5)) {
      if(runif(1) < eff_kids/100) { #if random < vaccine effectiveness, goes to recovered
        pool_pass[[i]] <- setAgentState(pool_pass[[i]],10, covid_input) 
      }
    }
  }
  
  #Pass teens selected
  for (i in pass_teens_selection) {
    if (!pool_pass[[i]]$biostate %in% c(3,4,5)) {
      if(runif(1) < eff_teens/100) { #if random < vaccine effectiveness, goes to recovered
        pool_pass[[i]] <- setAgentState(pool_pass[[i]],10, covid_input) 
      }
    }
  }
  
  
  #Pass Adults selected
  for (i in pass_adults_selection) {
    if (!pool_pass[[i]]$biostate %in% c(3,4,5)) {
      if(runif(1) < eff_adults/100) { #if random < vaccine effectiveness, goes to recovered
        pool_pass[[i]] <- setAgentState(pool_pass[[i]],10, covid_input) 
      }
    }
  }
  
  #Pass Seniors selected
  for (i in pass_seniors_selection) {
    if (!pool_pass[[i]]$biostate %in% c(3,4,5)) {
      if(runif(1) < eff_seniors/100) { #if random < vaccine effectiveness, goes to recovered
        pool_pass[[i]] <- setAgentState(pool_pass[[i]],10, covid_input) 
      }
    }
  }
  
  #Crew selected
  for (i in crew_selection) {
    if (!pool_crew[[i]]$biostate %in% c(3,4,5)) {
      if(runif(1) < eff_adults/100) { #if random < vaccine effectiveness, goes to recovered
        pool_crew[[i]] <- setAgentState(pool_crew[[i]],10, covid_input) 
      }
    }
  }
  
  pass_selection <- sort(c(pass_seniors_selection, pass_adults_selection, pass_teens_selection,  pass_kids_selection ))
  crew_selection <- sort(crew_selection)
  
  #table(passAgesNew[pass_selection])
  
  return(list(pool_pass = pool_pass,
              pool_crew = pool_crew,
              pass_selection = pass_selection,
              crew_selection = crew_selection))
}



### 9. interventions_plot: Plot selected agents for the interventions selected (masking or vaccination)

interventions_plot <- function(seed, coord, network, selection, main="", labels = c("selected", "not selected")) {
  
  # The range of indices you want to check
  range <- 1:dim(network)[1]
  
  # Create the new vector
  intervention_flag <- ifelse(range %in% selection, 1, 2)
  
  set.seed(seed)
  
  #create coordinates, if they are not specified in the function using the network matrix 
  if(is.null(coord))
  {
    coord  <- gplot.layout.fruchtermanreingold(network,layout.par=list(niter=500))
  }
  
  newmin <- mean(coord[,2]) - (-min(coord[,2]) + mean(coord[,2])) * 1.4
  palette =c("#fdb462", "#D3D3D3")
  
  #set plot layout
  plot(coord,col="black",bty="n",pch=16,cex=2.7,xaxt="n",yaxt="n",main=main,xlab="",ylab="",axes=F,
       ylim=c(newmin,max(coord[,2])),type="n")
  
  
  #plot segments 
  for(i in 1:nrow(network)) {
    segments(coord[i,1],
             coord[i,2],
             coord[network[i,]==1,1],
             coord[network[i,]==1,2],col="grey40")
    
  }
  
  
  #add circles
  points(coord,pch=16,cex=1,col="black", lwd = 1)
  
  #add colors based on states or ages
  
  points(coord,pch=16,cex=1 ,col = palette[intervention_flag], lwd = 1)
  
  legend(mean(coord[,1]),min(coord[,2]),bty='n',y.intersp=.7,cex=.8,
         labels, pch=16,col=palette)
  
  figure <- recordPlot()
  
  
  return(figure)
}

### 10. interventions_plot_CC: Plot selected agents for the testing intervention

interventions_plot_CC <- function(seed, coord, network, selection1, selection2, main="", labels = c("Tested - S", "Tested - CC", "Not tested")) {
  
  # The range of indices you want to check
  range <- 1:dim(network)[1]
  
  # Create the new vector
  intervention_flag <- case_when(range %in% selection1 ~ 1,
                                 range %in% selection2 ~ 2,
                                 TRUE ~ 3)
  
  set.seed(seed)
  
  #create coordinates, if they are not specified in the function using the network matrix 
  if(is.null(coord))
  {
    coord  <- gplot.layout.fruchtermanreingold(network,layout.par=list(niter=500))
  }
  
  newmin <- mean(coord[,2]) - (-min(coord[,2]) + mean(coord[,2])) * 1.4
  palette =c("#fdb462", "lightgreen", "#D3D3D3")
  
  #set plot layout
  plot(coord,col="black",bty="n",pch=16,cex=2.7,xaxt="n",yaxt="n",main=main,xlab="",ylab="",axes=F,
       ylim=c(newmin,max(coord[,2])),type="n")
  
  
  #plot segments 
  for(i in 1:nrow(network)) {
    segments(coord[i,1],
             coord[i,2],
             coord[network[i,]==1,1],
             coord[network[i,]==1,2],col="grey40")
    
  }
  
  
  #add circles
  points(coord,pch=16,cex=1,col="black", lwd = 1)
  
  #add colors based on states or ages
  
  points(coord,pch=16,cex=1 ,col = palette[intervention_flag], lwd = 1)
  
  legend(mean(coord[,1]),min(coord[,2]),bty='n',y.intersp=.7,cex=.8,
         labels, pch=16,col=palette)
  
  figure <- recordPlot()
  
  return(figure)
}



# #random
# selection_random <- masking_selection(std_net, selection = "random", coverage  = 50)
# selection_structured <- masking_selection(std_net, selection = "structured", coverage  = 5)
# 
# interventions_plot(seed = 123, coord = NULL, network = std_net$passNetwork, 
#                    selection = selection_random[[1]], main = "random")
# interventions_plot(seed = 123, coord = NULL, network = std_net$passNetwork, 
#                    selection = selection_structured[[1]], main = "structured")
# 
# #vax
# selection_random <- vaccination_selection(std_net, selection = "random", coverage  = c(2,5,10))
# selection_structured <-  vaccination_selection(std_net, selection = "structured", coverage  = c(2,5,10))
# interventions_plot(seed = 123, coord = NULL, network = std_net$passNetwork, 
#                     selection = random$pass_selection, main = "random")
# interventions_plot(seed = 123, coord = NULL, network = std_net$passNetwork, 
#                     selection = structured$pass_selection, main = "structured")


## D. Network functions 

### 1. makeDemographicNetwork: This function creates the demographic network of the ship. This includes creating passenger travel groups, crew roommate groups and assigning ages. 

makeDemographicNetwork <- function(cs_input, network_input , scenario, dist) { #takes demographic input from input spreadsheet 
  
  dist = dist 
  
  #get cruise ship settings 
  cs_input <- cs_input %>%
    filter(Scenario == scenario) 
  
  #get total number of passengers and crew 
  numPass = cs_input$`Number of passengers`
  numCrew = cs_input$`Number of crew`
  
  #Create crew cabin networks
  
  #crews have groups from 1 to 4
  crewGroup = seq(1,4,1)
  
  #probability of having groups of 1 to 4 for crew (this prob mean the % of crew that are in groups of 1, 2, 3, 4)
  crewGroupProb = c(cs_input$`Crew Groups of 1`, cs_input$`Crew Groups of 2`, cs_input$`Crew Groups of 3`, 
                    cs_input$`Crew Groups of 4`)
  
  #number of crew members in groups of 1,2,3,4
  crewGroupSize = crewGroupProb*numCrew
  
  #number of groups of each size (example: number of crew in groups of 2/group size), this tell you how many groups of each size there are
  #we round down in case the answer is not a whole number (example: 50 people in groups of 3/3 = 16.7 = 16 groups of 3 )
  crewNumGroup = floor(crewGroupSize/crewGroup)
  
  #if they are missing crew members (due to rounding down in the previous step), add them to groups of 2, if one is remaining, add them to groups of 1
  if (sum(crewNumGroup*crewGroup) != numCrew) {
    quotient =  (numCrew - sum(crewNumGroup*crewGroup))%/%2 #how many groups of 2 
    remainder = round((numCrew - sum(crewNumGroup*crewGroup))%%2) #how many groups of 1
    crewNumGroup[2] = crewNumGroup[2] + quotient #add new groups of 2
    crewNumGroup[1] = crewNumGroup[1] + remainder #add new groups of 3
  }
  
  
  #create vector to store group IDs
  crewGroupID = vector()
  
  #the following block create the agents in groups using the number of groups of each size
  
  #groups of 1
  #if there are groups of 1
  if (crewNumGroup[1] > 0) {
    crewGroup1 = seq(1, crewNumGroup[1],1) #create IDs
    crewGroupID = append(crewGroupID, crewGroup1) #append to running vector 
  }
  
  #groups of 2 
  #if there are groups of 2 
  if (crewNumGroup[2] > 0) {
    index = ifelse(length(crewGroupID)!= 0 , max(crewGroupID), 0) #get index from running crewGroupID vector
    crewGroup2 = rep(seq(index+1, index + crewNumGroup[2],1) ,2) #create IDs
    crewGroupID = append(crewGroupID, crewGroup2) #append to running vector 
    
  }
  
  #groups of 3 
  #if there are groups of 3 
  if (crewNumGroup[3] > 0) {
    index = ifelse(length(crewGroupID)!= 0, max(crewGroupID), 0) #get index from running crewGroupID vector
    crewGroup3 = rep(seq(index+1, index + crewNumGroup[3],1) ,3) #create IDs
    crewGroupID = append(crewGroupID, crewGroup3) #append to running vector 
    
  }
  
  #groups of 4
  #if there are groups of 4 
  if (crewNumGroup[4] > 0) {
    index = ifelse(length(crewGroupID)!= 0, max(crewGroupID), 0)  #get index from running crewGroupID vector
    crewGroup4 = rep(seq(index+1, index + crewNumGroup[4],1) ,4) #create IDs
    crewGroupID = append(crewGroupID, crewGroup4) #append to running vector 
    
  }
  
  #assign crew age (all crew members are in 18-64 age group)
  crewAge = rep("18-64", numCrew) #assign age 
  
  #check that number of IDs match numCrew
  if (length(crewGroupID) != numCrew) {
    stop("Check crew group assignemnts ") 
  }
  
  #sort IDS
  crewGroupID <- sort(crewGroupID)
  
  #get number of people per employment group
  num_FoodBack <- network_input[network_input[2] == "Dining Back - Crew",]
  num_FoodFront <- network_input[network_input[2] == "Dining Front - Crew",]
  num_Enter <- network_input[network_input[2] == "Entertainment - Crew",]
  num_House <- network_input[network_input[2] == "Housekeeping - Crew",]
  num_Other <- network_input[network_input[2] == "Other - Crew",]
  
  
  #food back network
  numberGroup = num_FoodBack$`Max groups`
  prop = num_FoodBack$`Prop included`
  
  crewFoodBackGroups = (ceiling(prop*numCrew/numberGroup)) #number of groups of 5 for food (back group)
  numFoodBack = crewFoodBackGroups*numberGroup
  crewFoodBackID <- sample(rep(c(1:crewFoodBackGroups), times = numberGroup)) #order groups ID randomly, this means that people in the same cabin
  #group are not necessarily in the same food group, but people in the same cabin group are likely to be in the same occupational role
  #gives zero to crew members that are not in this group
  crewFoodBackID <- c(crewFoodBackID, rep(0, times = numCrew - numFoodBack))
  
  #create food back matrix 
  crewFoodBackMat <- (outer(crewFoodBackID ,crewFoodBackID ,"==") * outer( crewFoodBackID , crewFoodBackID ,"*"))>0
  
  #create food front matrix
  numberGroup = num_FoodFront$`Max groups`
  prop = num_FoodFront$`Prop included`
  
  crewFoodFrontGroups = (ceiling(prop*numCrew/numberGroup)) #number of groups of 5 for food (front group)
  numFoodFront = crewFoodFrontGroups*numberGroup
  crewFoodFrontID <- sample(rep(c(1:crewFoodFrontGroups), times = numberGroup)) #order groups ID randomly, this means that people in the same cabin
  #group are not necessarily in the same food group, but people in the same cabin group are likely to be in the same occupational role
  #gives zero to crew members that are not in this group (include 0s in the begining for those who are food backcrew)
  crewFoodFrontID <- c(rep(0, times = numFoodBack), crewFoodFrontID, rep(0, times = numCrew - numFoodBack - numFoodFront))
  
  #create food back matrix 
  crewFoodFrontMat <- (outer(crewFoodFrontID  , crewFoodFrontID  ,"==") * outer( crewFoodFrontID  , crewFoodFrontID  ,"*"))>0
  
  
  #create housekeeping
  numberGroup = num_House$`Max groups`
  prop = num_House$`Prop included`
  
  housekeepingGroups = ceiling(prop*numCrew/numberGroup)
  numHousekeeping =   housekeepingGroups*numberGroup
  housekeepingID <- sample(rep(c(1: housekeepingGroups), times = numberGroup)) #order groups ID randomly 
  #gives zero to crew members that are not in this group (include 0s in the begining for those who are food back/front crew)
  housekeepingID <- c(rep(0, times = numFoodBack + numFoodFront),housekeepingID, rep(0, times = numCrew - numFoodBack - numFoodFront -numHousekeeping))
  
  #create housekeep matrix
  housekeepingMat <- (outer(housekeepingID ,housekeepingID ,"==") * outer( housekeepingID , housekeepingID ,"*"))>0
  
  #create entertainment group
  numberGroup = num_Enter$`Max groups`
  prop = num_Enter$`Prop included`
  
  entretainmentGroups = ceiling(prop*numCrew/numberGroup)
  numEntretainment =   entretainmentGroups*numberGroup
  entretainmentID <- sample(rep(c(1: entretainmentGroups), times = numberGroup)) #order groups ID randomly 
  #gives zero to crew members that are not in this group (include 0s in the beginning for those who are food back/front and housekeeping crew)
  entretainmentID <- c(rep(0, times = numFoodBack + numFoodFront + numHousekeeping),entretainmentID, rep(0, times = numCrew - numFoodBack - numFoodFront -numHousekeeping - numEntretainment))
  
  #create housekeep matrix
  entretainmentMat <- (outer(entretainmentID ,entretainmentID ,"==") * outer( entretainmentID , entretainmentID ,"*"))>0
  
  
  #create "other" group
  numberGroup = num_Other$`Max groups`
  
  #number of people in groups
  total = numFoodBack + numFoodFront + numEntretainment + numHousekeeping
  remaining = numCrew - total
  
  otherGroups = floor(remaining/numberGroup)
  numOther =   otherGroups*numberGroup
  otherID <- sample(rep(c(1:otherGroups), times = numberGroup)) #order groups ID randomly 
  #gives zero to crew members that are not in this group (include 0s in the beginning for those who are food back/front and housekeeping crew)
  otherID <- c(rep(0, times = numFoodBack + numFoodFront + numHousekeeping + numEntretainment),otherID, rep(0, times = max(0,numCrew - numFoodBack - numFoodFront -numHousekeeping - numEntretainment - numOther)))
  
  #create housekeep matrix
  otherMat <- (outer(otherID ,otherID  ,"==") * outer( otherID  , otherID  ,"*"))>0
  
  
  
  #assing roles based on proportions
  crewRole <- c(rep("food_back", numFoodBack),
                rep("food_front", numFoodFront),
                rep("housekeep", numHousekeeping),
                rep("entertainment", numEntretainment),
                rep("other", remaining))
  
  
  #Create passenger network
  
  ## groups of 1
  numGroups1 = dist$freq[dist$groupSize == 1]
  
  #get IDS
  passGroup1 = seq(1, sum(numGroups1),1)
  
  #get ages
  passGroup1_ages <- dist$comb[dist$groupSize == 1] #get age group distribution info for groups of 1 
  
  passAge_Group1 <- vector()
  
  #get passenger ages from the distribution dataframe 
  for (i in 1:length(passGroup1_ages)) {
    passAge_Group1 =  c(passAge_Group1, rep(passGroup1_ages[i], numGroups1[i]))
  }
  
  #groups of 2
  
  #get IDS
  index = max(passGroup1)  #get index from previous group
  numGroups2 = dist$freq[dist$groupSize == 2]  # number of groups of size 2
  passGroup2 = sort(rep(seq(index+1 , index+ sum(numGroups2),1),2)) #create IDS
  
  #get ages
  passGroup2_ages <- dist$comb[dist$groupSize == 2]
  
  passAge_Group2 <- vector()
  
  #assign ages
  for (i in 1:length(passGroup2_ages)) {
    comb = trimws(unlist(strsplit(passGroup2_ages[i], split=",")))
    passAge_Group2 =  c(passAge_Group2, rep(comb, numGroups2[i]))
  }
  
  
  #groups of 3 
  
  #get IDS
  index = max(passGroup2)  #get index from previous group
  numGroups3 = dist$freq[dist$groupSize == 3] # number of groups of size 3
  passGroup3 = sort(rep(seq(index+1 , index+ sum(numGroups3),1),3)) #create IDS
  
  
  #get ages
  passGroup3_ages <- dist$comb[dist$groupSize == 3]
  
  passAge_Group3 <- vector()
  
  #assign ages
  for (i in 1:length(passGroup3_ages)) {
    comb = trimws(unlist(strsplit(passGroup3_ages[i], split=",")))
    passAge_Group3 =  c(passAge_Group3, rep(comb, numGroups3[i]))
  }
  
  #groups of 4 
  
  #get IDS
  index = max(passGroup3)  #get index from previous group
  numGroups4 = dist$freq[dist$groupSize == 4]  # number of groups of size 4
  passGroup4 = sort(rep(seq(index+1 , index+ sum(numGroups4),1),4)) #get IDs
  
  
  #get ages
  passGroup4_ages <- dist$comb[dist$groupSize == 4]
  
  passAge_Group4 <- vector()
  
  #assign ages
  for (i in 1:length(passGroup4_ages)) {
    comb = trimws(unlist(strsplit(passGroup4_ages[i], split=",")))
    passAge_Group4 =  c(passAge_Group4, rep(comb, numGroups4[i]))
  }
  
  
  #groups of 5 
  
  #get IDS
  index = max(passGroup4)  #get index from previous group
  numGroups5 = dist$freq[dist$groupSize == 5] # number of groups of size 5
  passGroup5 = sort(rep(seq(index+1 , index+ sum(numGroups5),1),5)) #get IDS
  
  
  #get ages
  passGroup5_ages <- dist$comb[dist$groupSize == 5]
  
  passAge_Group5 <- vector()
  
  for (i in 1:length(passGroup5_ages)) {
    comb = trimws(unlist(strsplit(passGroup5_ages[i], split=",")))
    passAge_Group5 =  c(passAge_Group5, rep(comb, numGroups5[i]))
  }
  
  #add all
  passGroupID <- c(passGroup1,passGroup2, passGroup3, passGroup4,passGroup5)
  passAge <- c(passAge_Group1,passAge_Group2, passAge_Group3, passAge_Group4, passAge_Group5)
  
  #check that number of IDs match numPass
  
  if (length(passGroupID) != numPass) {
    stop("Check pass group assignemnts ") 
  }
  
  #check that lenght of IDs match age vector
  if (length(passGroupID) != length(passAge)) {
    stop("Check pass group age assignemnts ") 
  }
  
  
  ## lay out the social network
  passMat <- ( outer(passGroupID,passGroupID,"=="))>0 #create matrix size numCrew x numCrew , assign TRUE if crew members below to the same group
  crewMat <- ( outer(crewGroupID,crewGroupID,"=="))>0  #create matrix size numPass x numPass , assign TRUE if passengers  below to the same group
  
  return (list(
    passNetwork =passMat,
    passGroupID = passGroupID,
    crewNetwork = crewMat,
    crewGroupID = crewGroupID, 
    crewFoodBackNetwork = crewFoodBackMat,
    crewFoodFrontNetwork = crewFoodFrontMat,
    crewHousekeepingNetwork =  housekeepingMat,
    crewEntertainmentNetwork =  entretainmentMat,
    crewOtherNetwork = otherMat,
    passAge = passAge,
    crewAge = crewAge,
    crewRole = crewRole))
}

### 2. mygplot: This function plots passenger and crew network

mygplot <- function(seed, coord, network,age, states, employment,  plot = c("states", "ages", "employment"), main="") {
  
  plot = match.arg(plot)
  state_names <- c("susceptible", "exposed", "infected_asym", "infected_presymp" ,"infected_sym", "quarantine_infected", 
                   "quarantine_notinfected", "hospitalized", "dead", "recovered", "removed")
  
  state_labels <- c("Susceptible", "Exposed", "Infected \nasymptomatic", "Infected \npre-symptomatic", "Infected \nsymptomatic", 
                    "In isolation \ninfected", "In isolation \nnot \ninfected", "Hospitalized", "Dead", "Recovered", "Removed")
  
  set.seed(seed)
  
  #create coordinates, if they are not specified in the function using the network matrix 
  if(is.null(coord))
  {
    coord  <- gplot.layout.fruchtermanreingold(network,layout.par=list(niter=500))
  }
  
  newmin <- mean(coord[,2]) - (-min(coord[,2]) + mean(coord[,2])) * 1.4
  palette_state=c("#D3D3D3","#ffff99","#cab2d6", "#fdb462","#fb8072", "#1f78b4", "#a6cee3", "#ff7f00", "#e31a1c", "#33a02c","#6a3d9a" )
  palette_age = c("#7fc97f", "#beaed4", "#fdc086","#386cb0" )
  palette_employment = c("#7fc97f", "#beaed4", "#fdc086","#386cb0","#D3D3D3" )
  
  
  #set plot layout
  plot(coord,col="black",bty="n",pch=16,cex=2.7,xaxt="n",yaxt="n",main=main,xlab="",ylab="",axes=F,
       ylim=c(newmin,max(coord[,2])),type="n")
  
  
  #plot segments 
  for(i in 1:nrow(network)) {
    segments(coord[i,1],
             coord[i,2],
             coord[network[i,]==1,1],
             coord[network[i,]==1,2],col="grey40")
    
  }
  
  
  #add circles
  points(coord,pch=16,cex=1,col="black", lwd = 1)
  
  #add colors based on states or ages
  if (plot == "states") {
    points(coord,pch=16,cex=1,col= palette_state[states])
    
    legend(mean(coord[,1]),min(coord[,2]),bty='n',y.intersp=.7,cex=.8,
           state_names, pch=16,col=palette_state)
  } else if (plot == "ages") {
    age_group <- ifelse(age == "0-12", 1, ifelse(age == "13-17", 2, ifelse(age == "18-64", 3, 4)))
    age_labels <- c("0-12", "13-17", "18-64", "65+")
    points(coord,pch=16,cex=1,col= palette_age[age_group])
    
    legend(mean(coord[,1]),min(coord[,2]),bty='n',y.intersp=.7,cex=.8,
           age_labels, pch=16, col = palette_age) 
    
  } else {
    employment_group <- ifelse(employment == "food_back", 1, ifelse(employment == "food_front" , 2, 
                                                                    ifelse(employment == "housekeep" ,3, ifelse(employment == "entretainment", 4, 5))))
    employment_labels <- c("Food (Back)", "Food (Front)", "Housekeeping", "Entretainment", "Other")
    
    points(coord,pch=16,cex=1,col= palette_employment[employment_group])
    
    legend(mean(coord[,1]),min(coord[,2]),bty='n',y.intersp=.7,cex=.8,
           employment_labels, pch=16, col = palette_employment) 
    
  }
  
  figure <- recordPlot()
  
  
  return(figure)
}

### 3. activity_matrix_group_allages: Creates activity groups where travel groups mix with travel groups 

#takes: (1) Travel group matrix, (2) Information about the mixing: minimum and maximum number of groups mixing & probability of each,
# (3) Percentage of groups included

activity_matrix_group_allages <- function(net = net, min_groups = 1,
                                          max_groups = 4, proportions = proportions,
                                          prop_included = 1) {
  
  #number of travel groups (and account for losses)
  totalTravelGroups = length(unique(net$passGroupID)) #total
  travelGroups = floor(totalTravelGroups*prop_included) #accounting for only proportion that participates
  
  #get group sizes from min and max
  GroupsSize = seq(min_groups, max_groups)
  
  #group sizes (calculate how many groups of 1,2,...max_groups will be )
  multiplier = floor(travelGroups/sum(GroupsSize*proportions))
  numGroups = floor(multiplier*proportions) #number of groups of the sizes specified 
  
  #number of travel groups in each group size
  numTravelGroupsperGroup = numGroups*GroupsSize
  
  #check that travelGroups = sum of numTravelGroups_DiningGroup 
  diff = travelGroups - sum(numTravelGroupsperGroup)
  
  #add difference (we create groups of 1 with the difference)
  numGroups = numGroups + c(diff, rep(0, length(proportions)-1))
  numTravelGroupsperGroup = numGroups*GroupsSize
  
  #create activity group IDS based the number of groups per size
  #if else statement so we create an empty vector if proportion for that group size is 0 
  
  if (numTravelGroupsperGroup[1] != 0) {
    group1 = seq(1, length = numGroups[1])
  } else {
    group1 = NULL
  }
  
  if (numTravelGroupsperGroup[2]  != 0) {
    group2 = rep(seq(max(group1,0)+1, length = numGroups[2]),2)
  } else {
    group2 = NULL
  }
  
  if (numTravelGroupsperGroup[3]  != 0) {
    group3 = rep(seq(max(group1, group2,0)+1, length = numGroups[3]),3)
  } else {
    group3 = NULL
  }
  
  if (numTravelGroupsperGroup[4]  != 0) {
    group4 = rep(seq(max(group1, group2, group3,0)+1, length = numGroups[4]),4)
  } else {
    group4 = NULL
  }
  
  
  
  #group ID and travel ID cross-reference 
  cross = data.frame("passGroupID" = sample(1:totalTravelGroups, travelGroups), "GroupID" = c(group1,group2,group3, group4))
  
  #combine with original ID (pass groud by pass ID )
  comb = data.frame("passGroupID" = net$passGroupID) #create df with all travel group IDs
  comb = left_join(comb, cross,by = "passGroupID") #match travel group IDS with activity group IDs
  
  #fill with 0 those that did not participate
  comb <- replace(comb, is.na(comb), 0)
  
  #create matrix 
  matrix <- (outer(comb$GroupID, comb$GroupID ,"=="))*(outer(comb$GroupID,comb$GroupID ,"*"))>0
  
  return(matrix)
  
}


### 4. activity_matrix_ind_adults: Creates activity groups where adults with other adults 

#takes: (1) Travel group matrix, (2) Information about the mixing: minimum and maximum number of individuals mixing (assume uniform distribution)
# (3) Percentage of adults of 18-64 and 65+ included

activity_matrix_ind_adults <- function(net = net, min_contacts = 3,
                                       max_contacts = 10, 
                                       prop_18_64_included = 0.8,
                                       prop_65_p_included = 0.8) {
  
  #get all adults of hthe different ages
  totaladults_18_ID =  which(net$passAge %in% c("18-64"))
  totaladults_65_ID =  which(net$passAge %in% c("65+"))
  
  #get participating adults
  length_adults_18 = ceiling(length(totaladults_18_ID)*prop_18_64_included)
  length_adults_65 = ceiling(length(totaladults_65_ID)*prop_65_p_included)
  
  #sample participating adults
  part_adults_18_ID = sample(totaladults_18_ID, size = length_adults_18, replace = FALSE)
  part_adults_65_ID = sample(totaladults_65_ID, size = length_adults_65, replace = FALSE)
  
  numPass <- length(net$passAge)
  matrix_18 <- matrix(data = NA, nrow = numPass, ncol = numPass)
  matrix_65 <- matrix(data = NA, nrow = numPass, ncol = numPass)
  
  #complete matrix
  for (i in 1:numPass) {
    
    #if passenger row is part of the participant adults (18+), sample from participating adults
    if (i %in% part_adults_18_ID) {
      
      effective_max_contacts <- min(max_contacts, length(part_adults_18_ID))
      contacts <- sample(part_adults_18_ID, sample(min_contacts:effective_max_contacts, size = 1), replace = F)
      matrix_18[i,contacts] <- 1
      
    }
    
    #if passenger row is part of the participant adults (18+), sample from participating adults
    if (i %in% part_adults_65_ID) {
      
      effective_max_contacts <- min(max_contacts, length(part_adults_65_ID))
      contacts <- sample(part_adults_65_ID, sample(min_contacts:effective_max_contacts, size = 1), replace = F)
      matrix_65[i,contacts] <- 1
      
    }
  }
  
  matrix_18 <- replace( matrix_18 , is.na( matrix_18 ), 0)
  matrix_65 <- replace( matrix_65 , is.na( matrix_65 ), 0)
  
  matrix <- matrix_18 + matrix_65
  
  matrix  <- matrix == TRUE
  
  return(matrix)
  
  
}

### 5. activity_matrix_group_kids: Creates activity groups where kids mixed with other kids
#takes: (1) Travel group matrix, (2) Information about the mixing: minimum and maximum number of groups mixing & probability of each,
# (3) Percentage of groups included

activity_matrix_group_kids <- function(net = net, min_groups = 1,
                                       max_groups = 4, proportions = proportions,
                                       prop_included = 0.5) {
  
  #get travel groups that have kids by age group
  travelGroupsID_kids0 = unique(net$passGroupID[which(net$passAge == "0-12")])
  travelGroupsID_kids13 = unique(net$passGroupID[which(net$passAge == "13-17")])
  
  #number of travel groups with kids
  totalTravelGroups_kids0 = length(travelGroupsID_kids0)
  totalTravelGroups_kids13 = length(travelGroupsID_kids13)
  
  #participating travel groups with kids
  travelGroups_kids0 = floor(totalTravelGroups_kids0*prop_included)
  travelGroups_kids13 = floor(totalTravelGroups_kids13*prop_included)
  
  #get group sizes
  GroupsSize = seq(min_groups, max_groups)
  
  #group sizes (calculate how many groups of 1,2,...max_groups will be )
  multiplier_kids0 = floor(travelGroups_kids0/sum(GroupsSize*proportions))
  numGroups_kids0  = floor(multiplier_kids0*proportions)
  
  multiplier_kids13 = floor(travelGroups_kids13/sum(GroupsSize*proportions))
  numGroups_kids13  = floor(multiplier_kids13*proportions)
  
  #number of travel groups in each group size
  numTravelGroupsperGroup_kids0 = numGroups_kids0*GroupsSize
  numTravelGroupsperGroup_kids13 = numGroups_kids13*GroupsSize
  
  #check that travelGroups = numTravelGroups_DiningGroup 
  diff_kids0 = travelGroups_kids0 - sum(numTravelGroupsperGroup_kids0)
  diff_kids13 = travelGroups_kids13 - sum(numTravelGroupsperGroup_kids13)
  
  #add difference 
  numGroups_kids0 = numGroups_kids0 + c(diff_kids0, rep(0, length(proportions)-1))
  numTravelGroupsperGroup_kids0 = numGroups_kids0*GroupsSize
  
  numGroups_kids13 = numGroups_kids13 + c(diff_kids13, rep(0, length(proportions)-1))
  numTravelGroupsperGroup_kids13 = numGroups_kids13*GroupsSize
  
  
  #create group IDS
  
  #kids 0 -12 
  if (numTravelGroupsperGroup_kids0[1] != 0) {
    group1_kids0 = seq(1, length = numGroups_kids0[1])
  } else {
    group1_kids0 = NULL
  }
  
  if (numTravelGroupsperGroup_kids0[2]  != 0) {
    group2_kids0 = rep(seq(max(group1_kids0,0)+1, length = numGroups_kids0[2]),2)
  } else {
    group2_kids0 = NULL
  }
  
  
  if (numTravelGroupsperGroup_kids0[3]  != 0) {
    group3_kids0 = rep(seq(max(group1_kids0, group2_kids0,0)+1, length = numGroups_kids0[3]),3)
  } else {
    group3_kids0 = NULL
  }
  
  if (numTravelGroupsperGroup_kids0[4]  != 0) {
    group4_kids0 = rep(seq(max(group1_kids0, group2_kids0, group3_kids0,0)+1, length = numGroups_kids0[4]),4)
  } else {
    group4_kids0 = NULL
  }
  
  
  
  #kids 13-17
  if (numTravelGroupsperGroup_kids13[1] != 0) {
    group1_kids13 = seq(1, length = numGroups_kids13[1])
  } else {
    group1_kids13 = NULL
  }
  
  if (numTravelGroupsperGroup_kids13[2]  != 0) {
    group2_kids13 = rep(seq(max(group1_kids13,0)+1, length = numGroups_kids13[2]),2)
  } else {
    group2_kids13 = NULL
  }
  
  
  if (numTravelGroupsperGroup_kids13[3]  != 0) {
    group3_kids13 = rep(seq(max(group1_kids13, group2_kids13,0)+1, length = numGroups_kids13[3]),3)
  } else {
    group3_kids13 = NULL
  }
  
  if (numTravelGroupsperGroup_kids13[4]  != 0) {
    group4_kids13 = rep(seq(max(group1_kids13, group2_kids13,group3_kids13,0)+1, length = numGroups_kids13[4]),4)
  } else {
    group4_kids13 = NULL
  }
  
  
  
  #group ID and travel ID cross-reference 
  #sample passenger IDS from all kids using the number of participating kids
  #shuffle group IDS
  cross_kids0 = data.frame("passGroupID" = sample(travelGroupsID_kids0, size = travelGroups_kids0, replace = FALSE ), "GroupID" = c(group1_kids0,group2_kids0,group3_kids0, group4_kids0))
  cross_kids13 = data.frame("passGroupID" = sample(travelGroupsID_kids13, size = travelGroups_kids13, replace = FALSE ), "GroupID" = c(group1_kids13, group2_kids13,group3_kids13, group4_kids13))
  
  #combine with original ID (pass groud by pass ID )
  comb_kids0 = data.frame("passGroupID" = net$passGroupID, "passAge" = net$passAge)
  comb_kids0 = left_join(comb_kids0, cross_kids0,by = "passGroupID")
  
  comb_kids13 = data.frame("passGroupID" = net$passGroupID, "passAge" = net$passAge)
  comb_kids13 = left_join(comb_kids13, cross_kids13,by = "passGroupID")
  
  #replace with 0, individuals of different age that are in the group with kids 
  comb_kids0 <- comb_kids0 %>%
    mutate(GroupID = if_else(passAge == "0-12", GroupID, NA_real_))
  
  comb_kids13 <- comb_kids13 %>%
    mutate(GroupID = if_else(passAge == "13-17", GroupID, NA_real_))
  
  #fill with 0 those that did not participate
  comb_kids0 <- replace(comb_kids0, is.na(comb_kids0), 0)
  comb_kids13 <- replace(comb_kids13, is.na(comb_kids13), 0)
  
  #create matrix 
  matrix_kids0 <- (outer(comb_kids0$GroupID, comb_kids0$GroupID ,"=="))*(outer(comb_kids0$GroupID,comb_kids0$GroupID ,"*"))>0
  matrix_kids13 <- (outer(comb_kids13$GroupID, comb_kids13$GroupID ,"=="))*(outer(comb_kids13$GroupID,comb_kids13$GroupID ,"*"))>0
  
  #add two age groups 
  matrix <- 1*matrix_kids0 + 1*matrix_kids13
  
  matrix  <- matrix == TRUE
  
  return(matrix)
  
}

## E. Interaction functions

### 1. agents_interaction_sameclass: Model interactions of agents of the same class 

# Takes: (1) Pool of agents (crew or passengers), (2) Selected agent 1 and agent 2 based on the interactions, (3) contagion probability 
# Now it accounts for the individuals transmissibility 

agents_interaction_sameclass <- function(overall_pool, agent1_pool, agent2_pool, contagion_prob,
                                         masking_selected_IDs = NULL) {
  
  counter = 0 
  effective_contagion_prob = NA
  
  for(i in 1:length(agent1_pool)) {
    
    agent1 <- overall_pool[[agent1_pool[i]]] #get information of agent 1
    
    
    for (j in 1:length(agent2_pool[[i]])) {
      
      agent2 <-overall_pool[[agent2_pool[[i]][j]]] #get information of agent 1's interaction
      
      ##this constitutes the rules of infection.
      #if agent 1 is infected AND agent 2 is susceptible, contagion occurs based on contagion prob (we only consider when agent 1 is infected, since the agent 2 will eventually become agent 1 as the simulation progresses ),
      #if infection occurs, agent 2 becomes susceptible 
      if ((agent1$biostate %in% c(3,4,5)) & agent2$biostate == 1) {
        
        #get transmissibility factor
        #tf = getTransmissibilityFactor(agent1)
        tf = getTransmissibilityFactor_new(agent1)
        #tf = 1 
        
        #get reduction IF mask wearer 
        if (agent1_pool[i] %in% masking_selected_IDs) {
          masking_multiplier = 1-0.44
        } else {
          masking_multiplier = 1
        }
        
        #effective contagion probability for that day 
        #old
        #effective_contagion_prob = contagion_prob*tf*masking_multiplier
        
        #new
        scale = contagion_prob/tf$avg_viralLoad
        effective_contagion_prob = tf$viralLoad*scale*masking_multiplier
        #print(paste0("i= ", i," initial ", contagion_prob, " w/ multiplier final: ", contagion_prob*tf*masking_multiplier))
        
        #if random draw is below our prob 
        if (runif(1) < effective_contagion_prob) {
          
          ##infect agent 2 , becomes exposed
          overall_pool[[agent2_pool[[i]][j]]] <- setAgentState(agent2,2, covid_input) 
          
          #update counters
          counter = counter + 1
          
        }
        
      }
    }
  }
  
  #counter
  return(list(overall_pool_new = overall_pool,
              counter = counter,
              effective_contagion_prob =  effective_contagion_prob))
  
}


### 2. agents_interaction_twoclasses: Model interactions of agents of different classes

agents_interaction_twoclasses <- function(passenger_pool, crew_pool,  agent1_crew_pool, agent2_pass_pool, contagion_prob,
                                          masking_selected_pass = NULL,
                                          masking_selected_crew = NULL) {
  
  counter_crew = 0 
  counter_pass= 0 
  effective_contagion_prob = NA
  
  #crew to passenger direction 
  for(i in 1:length(agent1_crew_pool)) {
    
    agent1 <- crew_pool[[agent1_crew_pool[i]]] #get information of agent 1, crew 
    
    for (j in 1:length(agent2_pass_pool[[i]])) {
      
      agent2 <-passenger_pool[[agent2_pass_pool[[i]][j]]] #get information of agent's 1 interaction
      
      ##agent 1 infects agent 2 (crew to passenger)
      #if agent 1 is infected AND agent 2 is suspectible, contagion occurs based on contagion prob (we only consider when agent 1 is infected, since the agent 2 will eventually become agent 1 as the simulation progresses ),
      #if infection occurs, agent 2 becomes susceptible 
      if ( (agent1$biostate %in% c(3,4,5)) & agent2$biostate==1 ) {
        
        #get transmissibility factor
        #tf = getTransmissibilityFactor(agent1)
        #tf = 1 
        tf = getTransmissibilityFactor_new(agent1)
        
        #get reduction IF mask wearer 
        if (agent1_crew_pool[i] %in% masking_selected_crew) {
          masking_multiplier =  1-0.44
        } else {
          masking_multiplier = 1
        }
        
        #old
        #effective_contagion_prob = contagion_prob*tf*masking_multiplier
        #new
        scale = contagion_prob/tf$avg_viralLoad
        effective_contagion_prob = tf$viralLoad*scale*masking_multiplier
        
        #if random draw is below our prob 
        if (runif(1) <    effective_contagion_prob) {
          
          ##infect agent 2 , becomes exposed
          passenger_pool[[agent2_pass_pool[[i]][j]]] <- setAgentState(agent2,2, covid_input) ##infect agent 2 (passenger) , becomes exposed
          
          #update counters
          counter_pass = counter_pass + 1
          
        }
        
      } else if ((agent2$biostate %in% c(3,4,5)) & agent1$biostate==1) {
        # agent 2 infect agent 1 (passenger to crew)
        
        #get transmissibility factor
        #tf = getTransmissibilityFactor(agent2)
        #tf = 1 
        tf = getTransmissibilityFactor_new(agent2)
        
        #get reduction IF mask wearer 
        if (agent2_pass_pool[i] %in% masking_selected_pass) {
          masking_multiplier = 1-0.44
        } else {
          masking_multiplier = 1
        }
        
        #old
        #effective_contagion_prob = contagion_prob*tf*masking_multiplier
        #new
        scale = contagion_prob/tf$avg_viralLoad
        effective_contagion_prob = tf$viralLoad*scale*masking_multiplier
        
        #if random draw is below our prob 
        if (runif(1) < effective_contagion_prob) {
          
          ##infect agent 2 , becomes exposed
          crew_pool[[agent1_crew_pool[[i]]]] <- setAgentState(agent1,2, covid_input) ##infect agent 2 (crew) , becomes exposed
          
          #update counters
          counter_crew = counter_crew + 1
          
        }
        
        
      }
      
      
    }
  }
  
  return(list(passenger_pool_new = passenger_pool,
              crew_pool_new = crew_pool,
              counter_crew = counter_crew,
              counter_pass = counter_pass,
              effective_contagion_prob =  effective_contagion_prob))
  
}

## F. Simulation runs functions

## 1. simulation_initialize: initialize_simulation

# simulation_initialize: This function initializes simulation: Takes existing network and assigns infection state

#seed = fixed, means that the initial person(s) who are infected are fixed (will be the same if we use this function
#multiple times)
#groupseed = if it has number we fixed the group size of the initial infections 

# (A) seed = NULL, groupseed = NULL : Group size of infected individuals will be random as well as the IDs
# (B) seed = NULL, groupseed = 1,2,4,5 : Group size of infected individual will = groupseed, ID selceted will be random
# (C) seed = number, groupseed = NULL : The group size and the ID will follow the seed specified, the group size 
# of the infected individual won't change run to run, but could be 1-5
# (D) seed = number, groupseed = 1,2,3 4, 5 : Group size of infected individual will = groupseed, ID selected will follow 
# specified seed 

simulation_initialize <- function(net, cs_input, scenario = 1, covid_input, seed = 123 , groupseed = 2, plotNetwork) {
  
  #set seeds
  if (length(seed) > 0 ) {
    seed_num = seed
  } 
  
  if (length(seed) == 0 ) {
    seed_num = NULL
  }
  
  set.seed(seed_num)
  
  #filter correct scenario 
  cs_input <- cs_input %>%
    filter(Scenario == scenario )
  
  #get number of passengers and crew 
  numPass = cs_input$`Number of passengers`
  numCrew = cs_input$`Number of crew`
  
  #create list of agents
  pool_crew <- list()
  pool_pass <- list()
  
  #create crew pool (initially all are susceptible)
  for (i in 1:numCrew) {
    pool_crew[[i]] <- makeAgent(biostate= 1, 
                                age = net$crewAge[i])
  }
  
  #create passenger pool (initiall all are susceptible)
  for (i in 1:numPass) {
    pool_pass[[i]] <- makeAgent(biostate= 1, 
                                age = net$passAge[i])
  }
  
  #get the number of people who are initialy recovered/immune
  
  numRecovered_passenger <- cs_input$`Passenger - Recovered`
  numRecovered_crew <- cs_input$`Crew - Recovered`
  
  #if there is at least one recovered passenger
  if (numRecovered_passenger > 0) {
    
    #sample ID of those recovered
    sample_rec_passenger <- data.frame("ids" = sample(numPass, numRecovered_passenger, replace = FALSE),
                                       "state" = rep(10, numRecovered_passenger))
    
    #assign recovered state to ID 
    for (i in 1:nrow(sample_rec_passenger)) {
      
      pool_pass[[sample_rec_passenger[i, "ids"]]] <- setAgentState(pool_pass[[sample_rec_passenger[i, "ids"]]], sample_rec_passenger[i, "state"], covid_input)
      
    }
    
  }
  
  #if there is at least one recovered crew member 
  if (numRecovered_crew > 0) {
    
    
    #sample ID of those recovered
    sample_rec_crew <- data.frame("ids" = sample(numCrew, numRecovered_crew, replace = FALSE),
                                  "state" = rep(10, numRecovered_crew))
    
    #assign recovered state to ID 
    for (i in 1:nrow(sample_rec_crew)) {
      
      pool_crew[[sample_rec_crew[i, "ids"]]] <- setAgentState(pool_crew[[sample_rec_crew[i, "ids"]]], sample_rec_crew[i, "state"], covid_input)
      
    }
  }
  
  
  
  #initial infected 
  
  #passengers
  numInfected_passenger <- c(cs_input$`Passenger - Initial exposed`, 
                             cs_input$`Passenger - Initial infected asymptomatic`,
                             cs_input$`Passenger - Initial infected symptomatic`)
  
  numInfected_crew <- c(cs_input$`Crew - Initial exposed`, 
                        cs_input$`Crew - Initial infected asymptomatic`,
                        cs_input$`Crew - Initial infected symptomatic`)
  
  
  #if there's no groupseed (random)
  if (is.null(groupseed)) {
    
    #create eligible sample of passengers/crew who will be exposed, infected asymptomatic, and infected symptomatic
    
    #passengers
    if (numRecovered_passenger > 0) {
      eligible_pass <- setdiff(1:numPass,sample_rec_passenger$ids)
    } else {
      eligible_pass <- numPass
    }
    
    #crew
    if (numRecovered_crew > 0) {
      eligible_crew <- setdiff(1:numCrew,sample_rec_crew$ids)
    } else {
      eligible_crew <- numCrew
    }
    
    sample_passenger <- data.frame("ids" = sample(eligible_pass, sum(numInfected_passenger)),
                                   "state" = rep(c(2,3,4), numInfected_passenger)) #even if it's symptomatic, need to do pre-symp so we get complete viral load 
    
    sample_crew <- data.frame("ids" = sample(eligible_crew, sum(numInfected_crew)),
                              "state" = rep(c(2,3,4), numInfected_crew)) #even if it's symptomatic, need to do pre-symp so we get complete viral load 
    
    #passengers, set new states 
    if (nrow(sample_passenger) >= 1) {
      for (i in 1:nrow(sample_passenger)) {
        
        pool_pass[[sample_passenger[i, "ids"]]] <- setAgentState(pool_pass[[sample_passenger[i, "ids"]]], sample_passenger[i, "state"], covid_input)
        
      }
    }
    
    #crew, set new states 
    if (nrow(sample_crew) >= 1) {
      for (i in 1:nrow(sample_crew)) {
        
        pool_crew[[sample_crew[i, "ids"]]] <- setAgentState(pool_crew[[sample_crew[i, "ids"]]], sample_crew[i, "state"], 
                                                            covid_input)
      }
    }
    
    
    #if any were given pre-symptomatic, convert to symptomatic
    index_pass <- which(sapply(pool_pass,FUN=function(x){x$biostate}) == 4)
    index_crew <- which(sapply(pool_crew,FUN=function(x){x$biostate}) == 4)
    
    if (length(index_pass) > 0) {
      for (i in index_pass) {
        pool_pass[[i]] <- setAgentState(pool_pass[[i]] , 5, covid_input)
      }
    }
    
    if (length(index_crew) > 0) {
      for (i in index_crew) {
        pool_crew[[i]] <- setAgentState(pool_crew[[i]] , 5, covid_input)
      }
    }
    
  } else { #when groupseed is specified 
    
    #passengers
    #get a passenger ID where group of size = group seed , if more than one infected passenger, it might gives us different groups of passengers (i.e., infected)
    #passengers are in different groups 
    index <- sample(which(rowSums(net$passNetwork) == groupseed),sum(numInfected_passenger))
    
    #create dataframe with ID and state 
    sample_passenger <- data.frame("ids" = index, 
                                   "state" = rep(c(2,3,4), numInfected_passenger))
    
    #crew 
    #get a crew ID where group of size = 2 , if more than one infected passenger, it might gives us different groups
    #group seed is set to 2 for crew (as opposed to function argument for passenger)
    index <- sample(which(rowSums(net$crewNetwork) ==  2),sum(numInfected_crew))
    
    #get additional crew if needed
    sample_crew <- data.frame("ids" = index, 
                              "state" = rep(c(2,3,4), numInfected_crew))
    
    #assignd infections/exposed
    
    #passengers
    if (nrow(sample_passenger) >= 1) {
      for (i in 1:nrow(sample_passenger)) {
        
        pool_pass[[sample_passenger[i, "ids"]]] <- setAgentState(pool_pass[[sample_passenger[i, "ids"]]], sample_passenger[i, "state"], covid_input)
        
      }
      
    }
    
    #crew 
    if (nrow(sample_crew) >= 1) {
      for (i in 1:nrow(sample_crew)) {
        
        pool_crew[[sample_crew[i, "ids"]]] <- setAgentState(pool_crew[[sample_crew[i, "ids"]]], sample_crew[i, "state"], 
                                                            covid_input)
      } 
    }
    
    #if any were given pre-symptomatic, convert to symptomatic
    index_pass <- which(sapply(pool_pass,FUN=function(x){x$biostate}) == 4)
    index_crew <- which(sapply(pool_crew,FUN=function(x){x$biostate}) == 4)
    
    if (length(index_pass) > 0) {
      for (i in index_pass) {
        pool_pass[[i]] <- setAgentState(pool_pass[[i]] , 5, covid_input)
      }
    }
    
    if (length(index_crew) > 0) {
      for (i in index_crew) {
        pool_crew[[i]] <- setAgentState(pool_crew[[i]] , 5, covid_input)
      }
    }
    
  }
  
  
  
  
  
  #plot initial
  
  if(plotNetwork) {
    png(file = paste0("Results/Configuration/", scenario, "_passenger_day_", 0, ".png"), width = 1500, height = 1000)
    pass_initial <-mygplot(seed = 200, coord=NULL,net$passNetwork, net$passAge, sapply(pool_pass,FUN=function(x){x$biostate}) , plot = "states", main="Passenger - Day 0")
    dev.off() 
    
    png(file = paste0("Results/Configuration/", scenario, "_crew_day_", 0, ".png"), width = 1500, height =1000)
    crew_initial <-mygplot(seed = 200, coord=NULL,net$crewNetwork, net$crewAge, sapply(pool_crew,FUN=function(x){x$biostate}) , plot = "states", main="Crew - Day 0")
    dev.off() 
  }  
  
  plot_list <- list()
  
  return(list(pool_pass = pool_pass, pool_crew = pool_crew, sample_passenger = sample_passenger,
              sample_crew = sample_crew))
}

### 2. simulation_voyage_onerep: Runs one replication of the simulation

# Takes (1) network (crew and passenger), (2) initial set up , (3) cruise ship parameters, (4) covid parameters, (5) activitiy network parameters
# Selections: disembark = TRUE or FALSE , if people leave the ship some days
#observed_data = NULL or dataset, if we want to compute errors (so we can stop running if the error is too big)
# plotNetwork = TRUE or FALSE, if we want to plot the networks every day the simulation is running
# detailed_output = TRUE or FALSE, TRUE if we want to return dataset of each agent state per day (individual level)
# stop early = TRUE , stops at day 20 if error > 100

simulate_voyage_onerep <- function(net, initial, cs_input, scenario = 1, covid_input, network_input, observed_data = NULL,  plotNetwork = FALSE,
                                   detailed_output = FALSE, stop_early = FALSE) {
  
  #set.seed(seeds[1])
  
  #  -- * -- (1) Data Processing  -- * -- 
  
  #filter correct scenario 
  cs_input <- cs_input %>%
    filter(Scenario == scenario )
  
  network_input <- network_input %>%
    filter(Scenario == 1)
  
  #cruise ship characteristics
  numDays <- cs_input$`Trip duration (days)`
  numPass <-  cs_input$`Number of passengers`
  numCrew <- cs_input$`Number of crew`
  
  #netwok parameters
  dining_pass_parameters <- network_input[network_input[2] == "Dining - Pass",]
  enter_pass_parameters <- network_input[network_input[2] == "Entertainment - Pass",]
  night_pass_parameters <-   network_input[network_input[2] == "Nightlife - Pass",]
  kids_pass_parameters <-   network_input[network_input[2] == "Kids - Pass",]
  enter_pass_crew_parameters <- network_input[network_input[2] == "Entertainment - Pass and Crew",]
  dining_pass_crew_parameters <- network_input[network_input[2] == "Dining - Pass and Crew",]
  housekeeping_pass_crew_parameters <- network_input[network_input[2] == "Housekeeping - Pass and Crew",]
  
  #number of states
  states = 11
  
  #get initial states for passengers and crews 
  pool_pass = initial$pool_pass
  pool_crew = initial$pool_crew
  
  
  ## Interventions ##
  #set.seed(seeds[2])
  ## 1: Quarantine 
  
  quarantine = cs_input$Quarantine 
  
  if (quarantine == TRUE) {
    q_day <- cs_input$`Quarantine day`
    q_decrease_p <- cs_input$`Decrease in prob transmission`
    q_trigger <- cs_input$`Quarantine trigger`
  } 
  
  ## 2: Vaccination
  
  vaccination = cs_input$Vaccination
  
  if (vaccination == TRUE) {
    vax_selection <- cs_input$`Vaccination Selection`
    
    vax_coverage_kids <- cs_input$`Vax Coverage Kids`
    vax_coverage_teens <- cs_input$`Vax Coverage Teens`
    vax_coverage_adults <- cs_input$`Vax Coverage Adults`
    vax_coverage_seniors <- cs_input$`Vax Coverage Seniors`
    
    vax_effectiveness_kids <- cs_input$`Vax Effectiveness Kids`
    vax_effectiveness_teens <- cs_input$`Vax Effectiveness Teens`
    vax_effectiveness_adults <- cs_input$`Vax Effectiveness Adults`
    vax_effectiveness_seniors  <- cs_input$`Vax Effectiveness Seniors`
    
    #get vaccinated group
    vax_procedure <- vaccination_selection_implemention(network = net, 
                                                        pool_pass = pool_pass,
                                                        pool_crew = pool_crew, 
                                                        selection = vax_selection,
                                                        coverage = c(vax_coverage_kids, vax_coverage_teens,  vax_coverage_adults, vax_coverage_seniors),
                                                        effectiveness = c(vax_effectiveness_kids, vax_effectiveness_teens, vax_effectiveness_adults,
                                                                          vax_effectiveness_seniors))
    
    #update pool
    pool_pass = vax_procedure$pool_pass
    pool_crew = vax_procedure$pool_crew
    
    
  }
  
  ## 3: Testing
  
  testing = cs_input$Testing
  
  if (testing == TRUE) {
    test_protocol = cs_input$`Testing Protocol`
    test_capacity = cs_input$`Testing Capacity`
    test_days = cs_input$`Testing Days`
    test_fn_sym <- 1- cs_input$`Test TPR - Sym`
    test_fp_sym <- 1 - cs_input$`Test TNR - Sym`
    test_fn_asym <- 1 - cs_input$`Test TPR - Asym`
    test_fp_asym <- 1 - cs_input$`Test TNR - Asym`
    
  }
  
  ## 4: Masking
  
  masking = cs_input$Masking
  
  if (masking == TRUE) {
    
    #get inputs
    mask_selection = cs_input$`Masking Selection`
    mask_coverage_pass = cs_input$`Masking Coverage Passegers`
    mask_coverage_crew = cs_input$`Masking Coverage Crew`
    
    #people who use masks 
    masking_ID <- masking_selection(network = net, selection = mask_selection,
                                    coverage = c(mask_coverage_pass, mask_coverage_crew))
    
    passengers_masking <- masking_ID[[1]]
    crew_masking_emp <- masking_ID[[2]]
    
    
    
  } else {
    
    passengers_masking = NULL
    crew_masking_emp = NULL
  }
  
  
  #set.seed(seeds[3])
  #contation prob running 
  contagion_prob_running_main_pass <- vector()
  contagion_prob_running_main_crew <- vector()
  contagion_prob_running_act_pass <- vector()
  contagion_prob_running_act_kids <- vector()
  contagion_prob_running_emp_crew <- vector()
  contagion_prob_running_cp <- vector()
  
  #Contagion probability parameters from parameter sheet
  
  #Main group
  #Travel group
  contagionProb_pass_travel <- cs_input$`Passenger Contagion Probability - In network`
  #Crew cabin group 
  contagionProb_crew_cabin <- cs_input$`Crew Contagion Probability - In network`
  
  #Passenger activities 
  #Dining - Passenger network
  contagionProb_pass_dining <- cs_input$`Pass-Pass Dining Contagion Probability`
  #Entertainment - Passenger network
  contagionProb_pass_enter <- cs_input$`Pass-Pass Entertainment Contagion Probability`
  #Nightlife - Passenger network
  contagionProb_pass_night <- cs_input$`Pass-Pass Nightlife Contagion Probability`
  #Kids playgroups - Passenger network
  contagionProb_pass_kids <- cs_input$`Pass-Pass Kids Contagion Probability`
  
  #Crew employment network
  #Food - Crew network
  contagionProb_crew_food <- cs_input$`Crew-Crew Food Contagion Probability`
  #Housekeeping - Crew network
  contagionProb_crew_house <- cs_input$`Crew-Crew Housekeeping Contagion Probability`
  #Entertainment - Crew network
  contagionProb_crew_enter <- cs_input$`Crew-Crew Entertainment Contagion Probability`
  #Other - Crew network
  contagionProb_crew_other <- cs_input$`Crew-Crew Entertainment Contagion Probability`
  
  
  #Passenger-crew interactions based on employment networks 
  #Food - Pass-crew 
  contagionProb_crew_pass_food <- cs_input$`Crew-Passenger Food Contagion Probability`
  #Entertainment - Pass-crew
  contagionProb_crew_pass_enter <- cs_input$`Crew-Passenger Entertainment Contagion Probability`
  #Housekeeping - Pass-crew
  contagionProb_crew_pass_house <- cs_input$`Crew-Passenger Housekeeping Contagion Probability`
  
  # -- * -- 
  
  #  -- * -- (2) Set up variables to store results   -- * -- 
  
  # Create history of states (keeps track of the number of agents in each state daily)
  
  #get states for each agent at day 0 
  states_pass <- sapply(pool_pass,FUN=function(x){x$biostate}) #get states of all agents
  states_crew <- sapply(pool_crew,FUN=function(x){x$biostate}) #get states of all agents
  
  #create matrix that will store daily counts
  summaryhistory_pass <- matrix(c(0, table(factor(states_pass, levels = 1:11))), ncol = 12, byrow = TRUE)  #matrix of states x days 
  summaryhistory_crew <- matrix(c(0, table(factor(states_crew, levels = 1:11))), ncol = 12, byrow = TRUE)  #matrix of states x days 
  
  #create history of direction infections (who infected who)
  summaryhistory_direction <- matrix(NA, ncol = 4, nrow =numDays)
  
  #create history of in what network infections occur
  summaryhistory_network <- matrix(NA, ncol = 13, nrow =numDays)
  
  #create history number of agents exposed each day
  summaryhistory_exposed <- matrix(NA, ncol = 2, nrow =numDays)
  
  #create history of each agent per day (keep track of the agent state history throughout the voyage)
  indhistory_pass <- matrix(ncol = 3) #day, ID, state
  indhistory_crew <- matrix(ncol = 3) #day, ID, state
  
  #number of tests
  number_tests <- matrix(ncol = 2) #day, number of tests
  
  #keep track of operations 
  operations_track <- matrix(ncol = 2 )
  
  # -- * -- 
  
  #  -- * -- (2) Create interactions (agent 1 and agent 2) in the networks that do not change daily   -- * -- 
  
  # Networks that do not change daily are: 1.1 travel group, 1.2. cabin group
  # crew employment groups: 1.3 Food back, 1.4 Food front. 1.5. Housekeeping, 1.6. Entertainemnt 
  
  # 1.1 Travel networks
  
  #get matrix 
  matrixPass <- net$passNetwork
  diag(matrixPass) <- 0 
  
  #create sneezers (those that have at least 2 members if each group, at least one person in each matrix row)
  agent1_pass_travel <- which(rowSums(matrixPass) > 0)
  
  #create sneezed ons based on the travel groups matrix 
  agent2_pass_travel <- list() #total number of interactions (number of people x number of interactions per person)
  
  #each element on the list has the contacts for each sneezers 
  for(i in 1:length(agent1_pass_travel)) {
    agent2_pass_travel[[i]] <- which(matrixPass[agent1_pass_travel[i],] == 1)
  }
  
  
  # 1.2 Crew cabin networks 
  
  #get matrix 
  matrixCrew <- net$crewNetwork
  diag(matrixCrew) <- 0 
  
  #create sneezers (those that have at least 2 members if each group, at least one person in each matrix row)
  agent1_crew_cabin <- which(rowSums(matrixCrew) > 0)
  
  #create sneezed ons
  agent2_crew_cabin <- list() #total number of interactions (number of people x number of interactions per person)
  
  #each element on the list has the contacts for each sneezers
  for(i in 1:length(agent1_crew_cabin)) {
    agent2_crew_cabin[[i]] <- which(matrixCrew[agent1_crew_cabin[i],] == 1)
  }
  
  
  # 1.3 Food (back) crew-crew
  
  #get matrix
  matrixCrew <- net$crewFoodBackNetwork
  diag(matrixCrew) <- 0 
  
  #create sneezers (those that have at least 2 members if each group)
  agent1_crew_foodback <- which(rowSums(matrixCrew) > 0)
  
  #create sneezed ons
  agent2_crew_foodback<- list() 
  
  for(i in 1:length(agent1_crew_foodback)) {
    agent2_crew_foodback[[i]] <- which(matrixCrew[agent1_crew_foodback[i],] == 1)
  }
  
  # 1.4 Food (front) crew-crew
  
  #get matrix
  matrixCrew <- net$crewFoodFrontNetwork
  diag(matrixCrew) <- 0 
  
  #create sneezers (those that have at least 2 members if each group)
  agent1_crew_foodfront <- which(rowSums(matrixCrew) > 0)
  
  #create sneezed ons
  agent2_crew_foodfront = list() 
  
  for(i in 1:length(agent1_crew_foodfront)) {
    agent2_crew_foodfront[[i]] <- which(matrixCrew[agent1_crew_foodfront[i],] == 1)
  }
  
  #1.6 Entertainment (crew-to-crew)
  
  #assign interactions
  matrixCrew <- net$crewEntertainmentNetwork
  diag(matrixCrew) <- 0 
  
  #create sneezers (those that have at least 2 members if each group)
  agent1_crew_enter <- which(rowSums(matrixCrew) > 0)
  
  #create sneezed ons
  agent2_crew_enter <- list() 
  
  for(i in 1:length(agent1_crew_enter)) {
    agent2_crew_enter[[i]] <- which(matrixCrew[agent1_crew_enter[i],] == 1)
  }
  
  #1.5 Housekeeping crew-crew
  
  #get matrix
  matrixCrew <- net$crewHousekeeping
  diag(matrixCrew) <- 0 
  
  #create sneezers (those that have at least 2 members if each group)
  agent1_crew_house <- which(rowSums(matrixCrew) > 0)
  
  #create sneezed ons
  agent2_crew_house <- list() 
  
  for(i in 1:length(agent1_crew_house)) {
    agent2_crew_house[[i]] <- which(matrixCrew[agent1_crew_house[i],] == 1)
  }
  
  #1.6 Other crew-crew
  #get matrix
  matrixCrew <- net$crewOtherNetwork
  diag(matrixCrew) <- 0 
  
  #create sneezers (those that have at least 2 members if each group)
  agent1_crew_other <- which(rowSums(matrixCrew) > 0)
  
  #create sneezed ons
  agent2_crew_other <- list() 
  
  for(i in 1:length(agent1_crew_other)) {
    agent2_crew_other[[i]] <- which(matrixCrew[agent1_crew_other[i],] == 1)
  }
  
  
  
  #Identify housekeeping groups (for crew-passenger interactions)
  
  #identify groups of housekeeper 
  matrixCrew <- net$crewHousekeepingNetwork
  
  #create groups of 2 for crew passenger interactions
  numberGroup = 2
  numberHouse = length(which(net$crewRole == "housekeep"))
  
  before = length(which(net$crewRole %in% c("food_back", "food_front")))
  
  housekeepingGroups = floor(numberHouse/2)
  numHousekeeping =   housekeepingGroups*2
  housekeepingID <- sample(rep(c(1: housekeepingGroups), times = 2)) #order groups ID randomly 
  #gives zero to crew members that are not in this group (include 0s in the begining for those who are food back/front crew)
  housekeepingID <- c(rep(0, times = before),housekeepingID, rep(0, times = numCrew - (before + numHousekeeping)))
  
  #create housekeep matrix
  housekeepingMat_cp <- (outer(housekeepingID ,housekeepingID ,"==") * outer( housekeepingID , housekeepingID ,"*"))>0
  matrixCrew <- housekeepingMat_cp
  
  diag(matrixCrew) <- 0
  houseGroups <- as.data.frame(which(matrixCrew==TRUE, arr.ind = T)) #gives row and columns that = TRUE
  houseGroups$sorted_values <- apply(houseGroups , 1, function(row) paste(sort(row), collapse = ','))
  
  #the houseGroups dataframe shows the housecrew id and id of the crew member belonging to that group 
  houseGroups <- houseGroups %>%
    mutate(houseGroupID = match(sorted_values, unique(sorted_values))) %>%
    select(c(row, houseGroupID)) %>%
    rename(crewID = row)
  
  
  #Identify passenger groups
  #the travelGroups dataframe shows the travel group id and id of the passenger belonging to that group 
  travelGroups <- data.frame("passID" = 1:numPass, "travelGroupID" = net$passGroupID)
  
  # -- * -- 
  
  # Start day by day simulations
  #for(day in  28) {
  for(day in 1:numDays) {
    
    #print(day)
    #stop condition message
    stop_condition = "no error"
    
    #only if observed data exist
    if (!is.null(observed_data)) {
      observed_data_day <- observed_data %>%
        filter(Day == day)
    }
    
    
    #ONLY PERFORM KEEP DOING SIMULATION IF INFECTED IS MORE THAN 0 
    states_pass <- sapply(pool_pass,FUN=function(x){x$biostate}) #get states of all agents
    states_crew <- sapply(pool_crew,FUN=function(x){x$biostate}) #get states of all agents
    
    number_infectious <- length(which(c(states_pass, states_crew) %in% c(2,3,4,5)))
    
    if (number_infectious == 0) {
      stop_condition <- "Number of infectious = 0"
      break
    }
    
    #print(day) (comment in if we want to know in which day the simulation running is )
    
    #  -- * -- (3) Create networks that change daily and create interactions within each network  -- * -- 
    
    # These networks are 1.7. passenger dining, 1.8. entertainment, 1.9. nightlife,
    # 1.10 kids playgroups, 1.11 pass-crew food interactions, 1.12 pass-crew enter 
    # 1.13 housekeeping
    
    #get number of positive tests from end of day 
    number_positive <- length(which(c(states_pass, states_crew) %in% c(6,7)))
    
    #create quarantine flag
    if (quarantine == FALSE) {
      active_q = "inactive"
    } else if (quarantine == TRUE) {
      if (!is.na(q_day)) {
        active_q = if_else(day >= q_day , "active", "inactive")
      } else if (!is.na(q_trigger)) {
        active_q = if_else(number_positive >= ceiling(numPass*q_trigger/100) , "active", "inactive")
        
      }
    }
    
    #activate masking during quarantine
    if (!is.null(observed_data)) {
      if (day >= q_day) {
        
        #get inputs
        mask_selection = cs_input$`Masking Selection`
        mask_coverage_pass = cs_input$`Masking Coverage Passegers`
        mask_coverage_crew = cs_input$`Masking Coverage Crew`
        
        #people who use masks 
        masking_ID <- masking_selection(network = net, selection = mask_selection,
                                        coverage = c(mask_coverage_pass, mask_coverage_crew))
        
        passengers_masking <- masking_ID[[1]]
        crew_masking_emp <- masking_ID[[2]]
        
        
        
      } else {
        
        passengers_masking = NULL
        crew_masking_emp = NULL
      }
    }
    
    #track operations daily
    operations_track <- rbind(operations_track, c(day, active_q))
    
    # only create the networks + interactions if quarantine has not started (or not implemented)
    if (active_q == "inactive") {
      
      #print(paste0("Day ", day, ": Quarantine is not activated"))
      
      # A. Passenger networks
      
      # 1.7 Dining network
      
      #create network
      
      diningPassNetwork <- activity_matrix_group_allages(net = net, 
                                                         min_groups = dining_pass_parameters$`Min groups`,
                                                         max_groups = dining_pass_parameters$`Max groups`,
                                                         proportions = as.numeric(unlist(strsplit(dining_pass_parameters$`Group proportions`,","))),
                                                         prop_included = dining_pass_parameters$`Prop included`)
      
      #assign interactions 
      #get matrix
      matrixPass <- diningPassNetwork
      diag(matrixPass) <- 0
      
      #create sneezers (those that have at least 2 members if each group)
      agent1_pass_dining <- which(rowSums(matrixPass) > 0)
      
      #create sneezed ons
      agent2_pass_dining <- list() #total number of interactions (number of people x number of interactions per person)
      
      for(i in 1:length(agent1_pass_dining)) {
        agent2_pass_dining[[i]] <- which(matrixPass[agent1_pass_dining[i],] == 1)
      }
      
      # 1.8 Entertainment network 
      #create network
      enterPassNetwork <- activity_matrix_group_allages(net = net, 
                                                        min_groups = enter_pass_parameters$`Min groups`,
                                                        max_groups = enter_pass_parameters$`Max groups`,
                                                        proportions = as.numeric(unlist(strsplit(enter_pass_parameters$`Group proportions`,","))),
                                                        prop_included = enter_pass_parameters$`Prop included`)
      
      #assign interactions 
      #get matrix
      matrixPass <- enterPassNetwork
      diag(matrixPass) <- 0
      
      #create sneezers (those that have at least 2 members if each group)
      agent1_pass_enter <- which(rowSums(matrixPass) > 0)
      
      #create sneezed ons
      agent2_pass_enter <- list() #total number of interactions (number of people x number of interactions per person)
      
      for(i in 1:length(agent1_pass_enter)) {
        agent2_pass_enter[[i]] <- which(matrixPass[agent1_pass_enter[i],] == 1)
      }
      
      
      # 1.9 Nightlife 
      
      #create network
      nightPassNetwork <- activity_matrix_ind_adults(net = net,
                                                     min_contacts = night_pass_parameters$`Min contacts`,
                                                     max_contacts = night_pass_parameters$`Max contacts`,
                                                     prop_18_64_included = night_pass_parameters$`Prop included`,
                                                     prop_65_p_included = night_pass_parameters$`Prop included`)
      
      #assign interactions
      matrixPass <- nightPassNetwork
      diag(matrixPass) <- 0
      
      #create sneezers (those that have at least 2 members if each group)
      agent1_pass_night <- which(rowSums(matrixPass) > 0)
      
      #create sneezed ons
      agent2_pass_night <- list() #total number of interactions (number of people x number of interactions per person)
      
      for(i in 1:length(agent1_pass_night)) {
        agent2_pass_night[[i]] <- sample(agent1_pass_night, size = sample(night_pass_parameters$`Min contacts`:night_pass_parameters$`Max contacts`, size = 1), replace = FALSE) #draw from
        #the sneezers (all people in the network) from 2 to 10 (min and max contacts)
      }
      
      # 1.10. Kids group
      
      #ONLY IF THERE ARE KIDS 
      
      num_kids = length(which(net$passAge == "0-12" ))
      
      #create network
      if (num_kids > 0) {
        kidsPassNetwork <- activity_matrix_group_kids(net = net, 
                                                      min_groups = kids_pass_parameters$`Min groups`,
                                                      max_groups = kids_pass_parameters$`Max groups`,
                                                      proportions = as.numeric(unlist(strsplit(kids_pass_parameters$`Group proportions`,","))),
                                                      prop_included = kids_pass_parameters$`Prop included`)
        
        #assign interactions
        #get matrix 
        matrixPass <- kidsPassNetwork
        diag(matrixPass) <- 0
        
        #create sneezers (those that have at least 2 members if each group)
        agent1_pass_kids <- which(rowSums(matrixPass) > 0)
        
        #create sneezed ons
        agent2_pass_kids <- list() #total number of interactions (number of people x number of interactions per person)
        
        for(i in 1:length(agent1_pass_kids)) {
          agent2_pass_kids[[i]] <- which(matrixPass[agent1_pass_kids[i],] == 1)
          
        }
      }
      # 1.12 Entertainment 
      
      agent1_crew_pass_enter <- which(net$crewRole == "entertainment") #crew
      
      agent2_pass_crew_enter  <- list()
      
      for (i in 1:length(agent1_crew_pass_enter) ) {
        agent2_pass_crew_enter[[i]] <- sample(which(rowSums(enterPassNetwork) > 0), size = enter_pass_crew_parameters$`Max contacts`) #get 8 random passengers who participate in entertainment 
      }
      
    }
    
    
    # create all the networks + interaction that change every day and still happen when quarantine is in place 
    
    
    # 1.11 Food (front)
    
    agent1_crew_pass_foodfront <- which(net$crewRole == "food_front") #crew
    
    agent2_pass_crew_foodfront  <- list()
    
    for (i in 1:length(agent1_crew_pass_foodfront) ) {
      agent2_pass_crew_foodfront[[i]] <- sample(1:numPass, size = dining_pass_crew_parameters$`Max contacts`) #gets 8 (depending on the network ouptup file)
      
    }
    
    
    # 1.13 Housekeeping
    
    agent1_crew_pass_house <- houseGroups$crewID #crew
    
    #pass groups 
    travelGroupsID <- unique(travelGroups$travelGroupID)
    
    #house groups
    houseGroupsID <- unique(houseGroups$houseGroupID)
    
    ratio = floor(length(travelGroupsID)/length(houseGroupsID)) #number of travel groups per housekeeping group;8-9 per travel groups per pair of housekeepers
    vector <- sample(c(rep(houseGroupsID, ratio), seq(1, length(travelGroupsID) -  ratio*length(houseGroupsID)))) #create house to travel group assingment
    
    assignment <- data.frame("travelGroupsID" = travelGroupsID,
                             "houseGroupsID" = vector) #gets housekeeping group with assigned travel group
    
    #get 8 random from travel group
    
    selection_houseGroupID  <- list()
    
    #for each house group get 8 random people from the assigned travelgroups 
    for (i in houseGroupsID) {
      
      #get travel groups per houseID 
      travelGroups_selected = which(assignment$houseGroupsID == i)
      
      #number in travel group
      num = length(which(travelGroups$travelGroupID %in% travelGroups_selected))
      #get random  passengers
      selection <- sample(which(travelGroups$travelGroupID %in% travelGroups_selected), size =  housekeeping_pass_crew_parameters$`Max contacts`) #gets 4 or 8 depending if quarantine is activated
      
      
      #assigns selected passengers to house group
      selection_houseGroupID[[i]] <- selection
      
    }
    
    #assign passengers to each housekeeper 
    agent2_pass_crew_house <- list()
    index = 1
    
    for (i in agent1_crew_pass_house) {
      
      #get house group of agent 
      houseGroup_agent <- houseGroups[which(houseGroups$crewID == i), "houseGroupID"]
      
      #get passengers assinged to house group 
      selected <- selection_houseGroupID[[houseGroup_agent]]
      
      agent2_pass_crew_house[[index]] <- selected
      
      index = index + 1
      
    }
    
    # -- * -- 
    
    #  -- * -- (4) Create counters to keep track of infections  -- * -- 
    
    #counters 
    
    #direction
    counter_crew_to_pass <- 0
    counter_crew_to_crew <- 0
    counter_pass_to_crew <- 0
    counter_pass_to_pass <- 0
    
    #networks
    
    counter_pass_travel <- 0 #1.1
    counter_crew_cabin <- 0 #1.2 
    
    #crew roles networks
    counter_crew_dining_back <- 0 #1.3
    counter_crew_dining_front <- 0 #1.4
    counter_crew_house <- 0 #1.5
    counter_crew_enter <- 0 #1.6
    counter_crew_other <- 0
    
    #passenger activities
    counter_pass_dining <- 0 #1.7
    counter_pass_enter <- 0 #1.8
    counter_pass_night <- 0  #1.9
    counter_pass_kids <- 0 #1.10
    
    #pass-crew networks 
    counter_cp_dining <- 0  #1.11
    counter_cp_enter <- 0  #1.12
    counter_cp_house <- 0  #1.13
    
    
    #number of exposed 
    counter_exposed_pass <- 0
    counter_exposed_crew <- 0 
    
    #get vector that keeps track who has gotten tested throughout the cruise 
    list_tested_pass <- vector()
    list_tested_crew <- vector()
    
    #  -- * -- (5) Simulate interactions in the different networks  -- * -- 
    
    
    # SIMULATE INTERACTION 
    
    #for some networks the probability of transmission changes during quarantine, create a variable for this
    
    # no quarantine or quarantine has not started
    if (active_q == "inactive") { 
      prob_decrease = 0
      #print(paste0("Day ", day, ": Prob decrease is ", prob_decrease))
    } else if (active_q == "active") {
      prob_decrease = q_decrease_p
      #print(paste0("Day ", day, ": Prob decrease is ", prob_decrease))
      
      #quarantine active, masking 
      
    }
    
    ## Networks that are present any time ## 
    
    # 1.1 Passenger travel network
    
    #run interactions
    int_pass_travel <- agents_interaction_sameclass(pool_pass, agent1_pass_travel, agent2_pass_travel, contagionProb_pass_travel,
                                                    masking_selected_IDs = NULL)
    
    contagion_prob_running_main_pass<- c(contagion_prob_running_main_pass,int_pass_travel$effective_contagion_prob)
    
    #update passenger pool and counters if counter changed
    if (int_pass_travel$counter > 0) {
      pool_pass <- int_pass_travel$overall_pool_new
      counter_pass_travel = counter_pass_travel + int_pass_travel$counter
      counter_exposed_pass = counter_exposed_pass + int_pass_travel$counter
      counter_pass_to_pass = counter_pass_to_pass + int_pass_travel$counter
    }
    
    # 1.2 Crew cabin network
    
    #run interactions
    int_crew_cabin <- agents_interaction_sameclass(pool_crew, agent1_crew_cabin, agent2_crew_cabin, contagionProb_crew_cabin,
                                                   masking_selected_IDs = NULL)
    
    
    contagion_prob_running_main_crew <- c(contagion_prob_running_main_crew, int_crew_cabin$effective_contagion_prob)
    
    #update crew cabin pool and counters if counter changed
    if (int_crew_cabin$counter > 0) {
      pool_crew <- int_crew_cabin$overall_pool_new
      counter_crew_cabin = counter_crew_cabin + int_crew_cabin$counter
      counter_exposed_crew = counter_exposed_crew + int_crew_cabin$counter
      counter_crew_to_crew = counter_crew_to_crew + int_crew_cabin$counter
    }
    
    # 1.3 Food front, crew to crew (changes p during quarantine)
    
    #run interactions
    int_crew_foodfront <- agents_interaction_sameclass(pool_crew, agent1_crew_foodfront, agent2_crew_foodfront, contagionProb_crew_food*(1-prob_decrease),
                                                       crew_masking_emp)
    
    contagion_prob_running_emp_crew <- c( contagion_prob_running_emp_crew  , int_crew_foodfront$effective_contagion_prob)
    
    #update crew cabin pool and counters if counter changed
    if (int_crew_foodfront$counter > 0) {
      pool_crew <- int_crew_foodfront$overall_pool_new
      counter_crew_dining_front = counter_crew_dining_front + int_crew_foodfront$counter
      counter_exposed_crew = counter_exposed_crew + int_crew_foodfront$counter
      counter_crew_to_crew = counter_crew_to_crew + int_crew_foodfront$counter
    }
    
    # 1.4. Food back, crew to crew  (changes p during quarantine)
    
    #run interactions
    int_crew_foodback <- agents_interaction_sameclass(pool_crew, agent1_crew_foodback, agent2_crew_foodback, contagionProb_crew_food*(1-prob_decrease),
                                                      crew_masking_emp)
    
    contagion_prob_running_emp_crew <- c(contagion_prob_running_emp_crew, int_crew_foodback$effective_contagion_prob)
    
    #update crew cabin pool and counters if counter changed
    if (int_crew_foodback$counter > 0) {
      pool_crew <- int_crew_foodback$overall_pool_new
      counter_crew_dining_back = counter_crew_dining_back + int_crew_foodback$counter
      counter_exposed_crew = counter_exposed_crew + int_crew_foodback$counter
      counter_crew_to_crew = counter_crew_to_crew + int_crew_foodback$counter
    }
    
    
    # 1.5 Housekeeping : crew to crew (changes p during quarantine)
    
    #run interactions
    int_crew_house <- agents_interaction_sameclass(pool_crew, agent1_crew_house, agent2_crew_house, contagionProb_crew_house*(1-prob_decrease),
                                                   crew_masking_emp)
    
    contagion_prob_running_emp_crew <- c(contagion_prob_running_emp_crew, int_crew_house$effective_contagion_prob)
    
    #update crew cabin pool and counters if counter changed
    if (int_crew_house$counter > 0) {
      pool_crew <-  int_crew_house$overall_pool_new
      counter_crew_house =  counter_crew_house + int_crew_house$counter
      counter_exposed_crew = counter_exposed_crew + int_crew_house$counter
      counter_crew_to_crew = counter_crew_to_crew + int_crew_house$counter
    }
    
    # 1.6 Entertainment, crew to crew (changes p during quarantine)
    
    #run interactions
    int_crew_enter <- agents_interaction_sameclass(pool_crew, agent1_crew_enter, agent2_crew_enter, contagionProb_crew_enter*(1-prob_decrease),
                                                   crew_masking_emp)
    
    contagion_prob_running_emp_crew <- c(contagion_prob_running_emp_crew, int_crew_enter$effective_contagion_prob)
    
    
    #update crew cabin pool and counters if counter changed
    if (int_crew_enter$counter > 0) {
      pool_crew <-  int_crew_enter$overall_pool_new
      counter_crew_enter =  counter_crew_enter + int_crew_enter$counter
      counter_exposed_crew = counter_exposed_crew + int_crew_enter$counter
      counter_crew_to_crew = counter_crew_to_crew + int_crew_enter$counter
    }
    
    # Other, crew to crew (changes p during quarantine)
    
    #run interactions
    int_crew_other <- agents_interaction_sameclass(pool_crew, agent1_crew_other, agent2_crew_other, contagionProb_crew_other*(1-prob_decrease),
                                                   crew_masking_emp)
    
    contagion_prob_running_emp_crew <- c(contagion_prob_running_emp_crew, int_crew_other$effective_contagion_prob)
    
    
    #update crew cabin pool and counters if counter changed
    if (int_crew_other$counter > 0) {
      pool_crew <-  int_crew_other$overall_pool_new
      counter_crew_other =  counter_crew_other + int_crew_other$counter
      counter_exposed_crew = counter_exposed_crew + int_crew_other$counter
      counter_crew_to_crew = counter_crew_to_crew + int_crew_other$counter
    }
    
    
    
    #1.11 Front food, crew to pass (changes p during quarantine)
    
    int_crew_pass_food_front <- agents_interaction_twoclasses(pool_pass, pool_crew, agent1_crew_pass_foodfront, agent2_pass_crew_foodfront,
                                                              contagionProb_crew_pass_food*(1-prob_decrease),
                                                              masking_selected_pass = passengers_masking,
                                                              masking_selected_crew = crew_masking_emp)
    
    contagion_prob_running_cp <- c(contagion_prob_running_cp , int_crew_pass_food_front$effective_contagion_prob)
    
    
    #update crew and pass pool and counters if counter changed
    if (int_crew_pass_food_front$counter_crew + int_crew_pass_food_front$counter_pass > 0) {
      #update pools
      pool_pass <- int_crew_pass_food_front$passenger_pool_new
      pool_crew <- int_crew_pass_food_front$crew_pool_new
      
      #update counters
      counter_cp_dining <- counter_cp_dining + int_crew_pass_food_front$counter_crew + int_crew_pass_food_front$counter_pass
      counter_exposed_pass <-  counter_exposed_pass  + int_crew_pass_food_front$counter_pass
      counter_exposed_crew <-  counter_exposed_crew  + int_crew_pass_food_front$counter_crew
      counter_crew_to_pass <- counter_crew_to_pass + int_crew_pass_food_front$counter_pass
      counter_pass_to_crew <- counter_pass_to_crew + int_crew_pass_food_front$counter_crew
    }
    
    # 1.13 Housekeeping, crew to pass (changes p during quarantine)
    
    int_crew_pass_house <- agents_interaction_twoclasses(pool_pass, pool_crew, agent1_crew_pass_house, agent2_pass_crew_house,
                                                         contagionProb_crew_pass_house*(1-prob_decrease),
                                                         masking_selected_pass = passengers_masking,
                                                         masking_selected_crew = crew_masking_emp)
    
    contagion_prob_running_cp <- c(contagion_prob_running_cp , int_crew_pass_house$effective_contagion_prob)
    
    
    #update crew and pass pool and counters if counter changed
    if (int_crew_pass_house$counter_crew + int_crew_pass_house$counter_pass > 0) {
      
      pool_pass <- int_crew_pass_house$passenger_pool_new
      pool_crew <- int_crew_pass_house$crew_pool_new
      
      #update counters
      counter_cp_house <- counter_cp_house + int_crew_pass_house$counter_crew + int_crew_pass_house$counter_pass
      counter_exposed_pass <-  counter_exposed_pass  + int_crew_pass_house$counter_pass
      counter_exposed_crew <-  counter_exposed_crew  + int_crew_pass_house$counter_crew
      counter_crew_to_pass <- counter_crew_to_pass + int_crew_pass_house$counter_pass
      counter_pass_to_crew <- counter_pass_to_crew + int_crew_pass_house$counter_crew
      
    }
    
    
    ## Network that are only present before the quarantine ## 
    # only simulate if the quarantine has not started (or not implemented)
    if (active_q == "inactive") { 
      
      #1.7 Dining, pass to pass
      
      int_pass_dining <- agents_interaction_sameclass(pool_pass, agent1_pass_dining, agent2_pass_dining, contagionProb_pass_dining,
                                                      passengers_masking)
      
      contagion_prob_running_act_pass <- c(contagion_prob_running_act_pass , int_pass_dining$effective_contagion_prob)
      
      
      #update pass pool counter only if counter is updated 
      if (int_pass_dining$counter > 0) {
        #update passanger pool
        pool_pass <-   int_pass_dining$overall_pool_new
        
        #update counters
        counter_pass_dining  = counter_pass_dining  +  int_pass_dining$counter
        counter_exposed_pass = counter_exposed_pass +  int_pass_dining$counter
        counter_pass_to_pass = counter_pass_to_pass +  int_pass_dining$counter
      }
      
      #1.8 Entertainment , pass to pass
      
      int_pass_enter <- agents_interaction_sameclass(pool_pass, agent1_pass_enter, agent2_pass_enter, contagionProb_pass_enter,
                                                     passengers_masking)
      
      contagion_prob_running_act_pass  <- c(contagion_prob_running_act_pass , int_pass_enter$effective_contagion_prob)
      
      #update pass pool counter only if counter is updated 
      if (int_pass_enter$counter > 0) {
        #update passanger pool
        pool_pass <-   int_pass_enter$overall_pool_new
        
        #update counters
        counter_pass_enter  = counter_pass_enter  +  int_pass_enter$counter
        counter_exposed_pass = counter_exposed_pass +  int_pass_enter$counter
        counter_pass_to_pass = counter_pass_to_pass +  int_pass_enter$counter
      }
      
      #1.9 Nightlife , pass to pass 
      
      int_pass_night <- agents_interaction_sameclass(pool_pass, agent1_pass_night, agent2_pass_night, contagionProb_pass_night,
                                                     passengers_masking)
      
      
      contagion_prob_running_act_pass  <- c(contagion_prob_running_act_pass , int_pass_night$effective_contagion_prob)
      
      #update pass pool counter only if counter is updated 
      if (int_pass_night$counter > 0) {
        
        #update passanger pool
        pool_pass <-   int_pass_night$overall_pool_new
        
        #update counters
        counter_pass_night  = counter_pass_night  +  int_pass_night$counter
        counter_exposed_pass = counter_exposed_pass +  int_pass_night$counter
        counter_pass_to_pass = counter_pass_to_pass +  int_pass_night$counter
      }
      
      
      #1.10 Kids, pass to pass 
      
      if (num_kids > 0) {
        int_pass_kids <- agents_interaction_sameclass(pool_pass, agent1_pass_kids, agent2_pass_kids, contagionProb_pass_kids,
                                                      passengers_masking)
        contagion_prob_running_act_kids <- c(contagion_prob_running_act_kids, int_pass_kids$effective_contagion_prob)
        
        
        #update pass pool counter only if counter is updated 
        if (int_pass_kids$counter > 0) {
          #update passanger pool
          pool_pass <-   int_pass_kids$overall_pool_new
          
          #update counters
          counter_pass_kids = counter_pass_kids  +  int_pass_kids$counter
          counter_exposed_pass = counter_exposed_pass +  int_pass_kids$counter
          counter_pass_to_pass = counter_pass_to_pass +  int_pass_kids$counter
        }
      }
      #1.12 Entertainment, crew to pass
      
      int_crew_pass_enter <- agents_interaction_twoclasses(pool_pass, pool_crew, agent1_crew_pass_enter, agent2_pass_crew_enter,
                                                           contagionProb_crew_pass_enter,
                                                           masking_selected_pass = passengers_masking,
                                                           masking_selected_crew = crew_masking_emp)
      
      contagion_prob_running_cp <- c(contagion_prob_running_cp, int_crew_pass_enter$effective_contagion_prob)
      
      #update pass pool counter only if counter is updated 
      if (int_crew_pass_enter$counter_crew  + int_crew_pass_enter$counter_pass > 0) {
        #update pools
        pool_pass <- int_crew_pass_enter$passenger_pool_new
        pool_crew <- int_crew_pass_enter$crew_pool_new
        
        #update counters
        counter_cp_enter <- counter_cp_enter + int_crew_pass_enter$counter_crew + int_crew_pass_enter$counter_pass
        counter_exposed_pass <-  counter_exposed_pass  + int_crew_pass_enter$counter_pass
        counter_exposed_crew <-  counter_exposed_crew  + int_crew_pass_enter$counter_crew
        counter_crew_to_pass <- counter_crew_to_pass + int_crew_pass_enter$counter_pass
        counter_pass_to_crew <- counter_pass_to_crew + int_crew_pass_enter$counter_crew
      }
      
    }
    
    
    #  -- * -- (6) Implement daily testing  -- * -- 
    
    
    #if getting data from the observed data set 
    if (!is.null(observed_data)) {
      
      #daily testing is on on/after quarantine starts 
      if (active_q == "active") {
        
        tests = observed_data_day$`Number of tests`
        
        #only test is not na
        if (!is.na(tests)) {
          
          # First, we identify IDs of the passengers/crew that will get tested based on the priority 
          #testIDs <- getTestIDs_observed(numberTest = tests, pool_pass, pool_crew)
          
          #First test S+CC
          testIDs <- getTestIDs_S_CC(numberTest = tests, pool_pass = pool_pass, pool_crew = pool_crew, net = net)
          
          test_left <- testIDs$numberTest
          
          #Anybody else
          if (test_left > 0) {
            add_testIDs <- getTestIDs_asym(numberTest = test_left, pool_pass = pool_pass, pool_crew = pool_crew, 
                                           tested_pass_prev = testIDs$tested_pass,
                                           tested_crew_prev = testIDs$tested_crew)
            
            #add
            testIDs$tested_pass <- c(testIDs$tested_pass, add_testIDs$tested_pass)
            testIDs$tested_crew <- c(testIDs$tested_crew, add_testIDs$tested_crew)
          }
          # Second, using the IDs selected, we send to quarantine those with the positive test results 
          testResults <- testDaily(tested_pass = testIDs$tested_pass, tested_crew = testIDs$tested_crew, pool_pass = pool_pass,
                                   pool_crew = pool_crew) 
          
          # Update pool of passengers and crew
          pool_pass <- testResults$pool_pass
          pool_crew <- testResults$pool_crew
          
        }
      }
      
    } else if (testing == TRUE) { #no observed data, but testing is happening 
      
      #get days
      test_days_vector <- as.numeric(strsplit(test_days, ",")[[1]])
      
      #if it is day of testing
      if (day %in% test_days_vector ) {
        
        tests = test_capacity
        
        #only symptomatic 
        if (test_protocol == "S") {
          
          # First, we identify IDs of the passengers/crew that will get tested based on the priority 
          testIDs <- getTestIDs_S(numberTest = tests, pool_pass, pool_crew)
          
          #Remove people who have been tested before
          pass_tested <- setdiff(testIDs$tested_pass, list_tested_pass)
          crew_tested <- setdiff(testIDs$tested_crew, list_tested_crew)
          
          #update
          testIDs$tested_pass <- pass_tested
          testIDs$tested_crew <- crew_tested
          
          # Second, using the IDs selected, we send to quarantine those with the positive test results 
          testResults <- testDaily(tested_pass = testIDs$tested_pass, tested_crew = testIDs$tested_crew, pool_pass = pool_pass,
                                   pool_crew = pool_crew, fn_error_sym= test_fn_sym, fp_error_sym = test_fp_sym,
                                   fn_error_asym= test_fn_asym, fp_error_asym = test_fp_asym)
          
          # Update pool of passengers and crew
          pool_pass <- testResults$pool_pass
          pool_crew <- testResults$pool_crew
          
          # Update list of tested
          list_tested_pass <- c(list_tested_pass, testIDs$tested_pass)
          list_tested_crew <- c(list_tested_crew, testIDs$tested_crew)
          
          ntest = length(testIDs$tested_crew) + length(testIDs$tested_pass)
          
          #Save tests
          number_tests <- rbind(number_tests, c(day, ntest))
          
        } else if (test_protocol == "S+CC") {
          
          # First, we identify IDs of the passengers/crew that will get tested based on the priority 
          testIDs <- getTestIDs_S_CC(numberTest = tests, pool_pass, pool_crew, net = net)
          
          #Remove people who have been tested before
          pass_tested <- setdiff(testIDs$tested_pass, list_tested_pass)
          crew_tested <- setdiff(testIDs$tested_crew, list_tested_crew)
          
          #update
          testIDs$tested_pass <- pass_tested
          testIDs$tested_crew <- crew_tested
          
          # Second, using the IDs selected, we send to quarantine those with the positive test results 
          testResults <- testDaily(tested_pass = testIDs$tested_pass, tested_crew = testIDs$tested_crew, pool_pass = pool_pass,
                                   pool_crew = pool_crew,  fn_error_sym= test_fn_sym, fp_error_sym = test_fp_sym,
                                   fn_error_asym= test_fn_asym, fp_error_asym = test_fp_asym)
          
          
          # Update pool of passengers and crew
          pool_pass <- testResults$pool_pass
          pool_crew <- testResults$pool_crew
          
          # Update list of tested
          list_tested_pass <- c(list_tested_pass, testIDs$tested_pass)
          list_tested_crew <- c(list_tested_crew, testIDs$tested_crew)
          
          
          ntest = length(testIDs$tested_crew) + length(testIDs$tested_pass)
          
          #Save tests
          number_tests <- rbind(number_tests, c(day, ntest))
          
          
        } else if (test_protocol == "S+CC_D") {
          
          ntest = 0 
          
          #only from second day of testing 
          if (day != test_days_vector[1] & exists("nextday_pass")) {
            
            testResults <- testDaily(tested_pass = nextday_pass, tested_crew = nextday_crew, pool_pass = pool_pass,
                                     pool_crew = pool_crew, fn_error = test_fn, fp_error = test_fp)
            
            # Update pool of passengers and crew
            pool_pass <- testResults$pool_pass
            pool_crew <- testResults$pool_crew
            
            ntest = length(nextday_pass) + length(nextday_crew)
          }
          
          
          #symptomatic 
          
          # First, we identify IDs of the passengers/crew that will get tested based on the priority 
          testIDs <- getTestIDs_S_CC_D(numberTest = tests, pool_pass, pool_crew, net = net)
          
          # Second, using the IDs selected, we send to quarantine those with the positive test results 
          testResults <- testDaily(tested_pass = testIDs$tested_pass, tested_crew = testIDs$tested_crew, pool_pass = pool_pass,
                                   pool_crew = pool_crew, fn_error = test_fn, fp_error = test_fp)
          
          # Update pool of passengers and crew
          pool_pass <- testResults$pool_pass
          pool_crew <- testResults$pool_crew
          
          #Close contacts 
          testResults <- testDaily(tested_pass = testIDs$nextday_pass, tested_crew = testIDs$nextday_crew, pool_pass = pool_pass,
                                   pool_crew = pool_crew, fn_error = test_fn, fp_error = test_fp)
          
          # Update pool of passengers and crew
          pool_pass <- testResults$pool_pass
          pool_crew <- testResults$pool_crew
          
          #Vector for next time
          nextday_pass <- setdiff(testIDs$nextday_pass,  which(sapply(pool_pass,FUN=function(x){x$biostate}) %in% c(6,7)))
          nextday_crew <- setdiff(testIDs$nextday_crew,  which(sapply(pool_crew,FUN=function(x){x$biostate}) %in% c(6,7)))
          
          #Save tests
          
          ntest = ntest + length( testIDs$tested_pass) + length(testIDs$tested_crew) + length(testIDs$nextday_pass) + length(testIDs$nextday_crew)
          number_tests <- rbind(number_tests, c(day, ntest))
          
          
        }
        
      }
      
      
      
    }
    
    
    
    #  -- * -- (7) Remove passengers (if indicated)  -- * -- 
    
    
    #if getting data from the observed data set 
    if (!is.null(observed_data)) {
      if (!is.na(observed_data_day$Disembark)) {
        #number of passengers removed 
        removed = observed_data_day$Disembark
        #removal process
        removedPax <- removePassengers(numberRemoved = removed, pool_pass = pool_pass)
        #update pool of passengers
        pool_pass <- removedPax
        
      }
    }
    
    #  -- * -- (8) Update state counters at the end of the day  -- * -- 
    
    # Increment agent one day for state counters
    # Transition to next state if counter is 0
    for(i in 1:numPass) {
      #print(paste0("day ", day, " ", i))
      pool_pass[[i]] <- updateAgent(pool_pass[[i]], covid_input)
    }
    
    for(i in 1:numCrew) {
      pool_crew[[i]] <- updateAgent(pool_crew[[i]], covid_input)
    }
    
    #get states for each agent
    states_pass <- sapply(pool_pass,FUN=function(x){x$biostate}) #get states of all agents
    states_crew <- sapply(pool_crew,FUN=function(x){x$biostate}) #get states of all agents
    
    #get daily new exposed, new infections
    
    #identify
    selected_pass <- which(states_pass %in% c(2,3,4,5))
    selected_crew <- which(states_crew %in% c(2,3,4,5))
    
    #Save individual passenger state by day (Exposed, Infected symtomatic, Infected asymtomatic)
    if (length(selected_pass) > 0) {
      indhistory_pass <- rbind(indhistory_pass, matrix(c(rep(day, length(selected_pass)),  selected_pass, states_pass[selected_pass]), byrow = FALSE, nrow = length(selected_pass)))
    }
    
    if (length(selected_crew) > 0) {
      indhistory_crew <- rbind(indhistory_crew, matrix(c(rep(day, length(selected_crew)),  selected_crew, states_crew[selected_crew]), byrow = FALSE, nrow = length(selected_crew)))
    }
    
    #Summary per day: Get daily counts in each state 
    distrib_pass <- table(factor(states_pass,levels=1:11)) #count how many in each state 
    distrib_crew <- table(factor(states_crew,levels=1:11)) #count how many in each state 
    
    #Add new day
    summaryhistory_pass <- rbind(summaryhistory_pass, c(day, distrib_pass)) #get distribtution hisotry for the day 
    summaryhistory_crew <- rbind(summaryhistory_crew, c(day, distrib_crew)) #get distribtution hisotry for the day 
    
    if (detailed_output == TRUE) {
      
      
      #counter per day
      summaryhistory_direction[day,] <- c(counter_crew_to_pass, counter_crew_to_crew, 
                                          counter_pass_to_crew, counter_pass_to_pass)
      
      summaryhistory_network[day,] <- c(counter_pass_travel, #1
                                        counter_crew_cabin,  #2 
                                        counter_crew_dining_back, #3 
                                        counter_crew_dining_front, #4
                                        counter_crew_house, #5
                                        counter_crew_enter, #6
                                        counter_pass_dining, #7
                                        counter_pass_enter, #8
                                        counter_pass_night, #9
                                        counter_pass_kids, #10
                                        counter_cp_dining, #11
                                        counter_cp_enter, #12
                                        counter_cp_house) #13
      
      #counter new exposed
      summaryhistory_exposed[day,] <- c(counter_exposed_crew, counter_exposed_pass)
      
    }
    
    
    
    #  -- * -- (9) Plot networks  -- * -- 
    
    
    if(plotNetwork){
      
      png(file = paste0("Results/Simulation/", scenario,"_passenger_day_", day, ".png"), width = 1000, height = 800)
      mygplot(seed = 200, coord=NULL, net$passNetwork, net$passAge, states_pass, plot = "states", main=paste("Pass - Day",day))
      dev.off() 
      #mygplot(seed = 100, coord=NULL, net$passNetwork, net$passAge, states_pass, plot = "states", main=paste("Pass - Day",day))
      
      png(file = paste0("Results/Simulation/", scenario, "_crew_day_", day, ".png"), width = 1000, height = 800)
      mygplot(seed = 200, coord=NULL, net$crewNetwork, net$crewAge, states_crew, plot = "states", main=paste("Crew - Day",day))
      dev.off() 
      #mygplot(seed = 100, coord=NULL, net$crewNetwork, net$crewAge, states_crew, plot = "states", main=paste("Crew - Day",day))
      
    }  
    
    
    
    
    #  -- * -- (10) See if we need to skip rest of simulation  -- * -- 
    
    if (stop_early == TRUE & day == 20) {
      
      cum_detected_sim = summaryhistory_pass[which(summaryhistory_pass[,1] == 20), "6"] + summaryhistory_crew[which(summaryhistory_pass[,1] == 20), "6"]
      error = abs(cum_detected_sim - observed_data_day$`Cumulative Detected Infections`)
      
      if (error > 300) {
        stop_condition <- "error is greater than 100"
        break
      }
    }
    print(paste0("Day ", day, " Scenario ", scenario))
  }
  
  
  #  -- * -- (11) Post-simulation data processing  -- * -- 
  
  #Calculate unique exposed, asymptomatic, pre-symptomatic and symptomatic 
  
  indhistory_pass_summary <- as.data.frame(indhistory_pass) %>%
    filter(!is.na(V1)) %>%
    rename(day = V1, ID = V2, state = V3) %>%
    distinct(ID, state, .keep_all = TRUE) %>%
    group_by(day, state) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    mutate(state = case_when(state== 2 ~ "E_pass",
                             state == 3 ~ "Ia_pass",
                             state == 4 ~ "Ipres_pass",
                             state == 5 ~ "Is_pass")) %>%
    pivot_wider(names_from = state, values_from = count) 
  
  
  #crew
  indhistory_crew_summary <- as.data.frame(indhistory_crew) %>%
    filter(!is.na(V1)) %>%
    rename(day = V1, ID = V2, state = V3) %>%
    distinct(ID, state, .keep_all = TRUE) %>%
    group_by(day, state) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    mutate(state = case_when(state== 2 ~ "E_crew",
                             state == 3 ~ "Ia_crew",
                             state == 4 ~ "Ipres_crew",
                             state == 5 ~ "Is_crew")) %>%
    pivot_wider(names_from = state, values_from = count) 
  
  #  -- * -- (12) Calculate error -- * -- 
  
  #only if we have observed data
  
  
  #if getting data from the observed data set, calculate errors 
  if (!is.null(observed_data)) {
    
    #print(paste0("Day ", day, ": Calculating error df"))
    
    #create error df
    totaldays = length(summaryhistory_pass[,1])
    
    #get cumulative
    error_df = data.frame("day" = summaryhistory_pass[,1],
                          "cum_simulated_detections" =  (summaryhistory_pass[,"6"]+ summaryhistory_crew[,"6"]),
                          "cum_observed_detections" = pull(observed_data[1:totaldays, "Cumulative Detected Infections"]),
                          "daily_observed_detections" =  pull(observed_data[1:totaldays, "Number of positive tests by day"]))
    
    #get daily
    error_df <- error_df %>%
      mutate(daily_simulated_detections = cum_simulated_detections - lag(cum_simulated_detections, n = 1L))
    
    #calculations
    error_df$cum_uw_error = (error_df$cum_simulated_detections - error_df$cum_observed_detections)
    error_df$cum_w_error =  pmax(error_df$cum_simulated_detections - error_df$cum_observed_detections, 0)  + 4*pmax( error_df$cum_observed_detections -  error_df$cum_simulated_detections,0)
    error_df$daily_uw_error = (error_df$daily_simulated_detections - error_df$daily_observed_detections)
    error_df$daily_w_error =  pmax(error_df$daily_simulated_detections - error_df$daily_observed_detections, 0)  + 4*pmax( error_df$daily_observed_detections -  error_df$daily_simulated_detections,0)
    
    #error_summary 
    
    #only if 27 days are collected 
    
    if (stop_condition == "no error") {
      uw_error_agg = sqrt(sum(error_df$daily_uw_error^2, na.rm =TRUE))
      #w_error_agg = sqrt(sum(error_df$daily_w_error^2, na.rm =TRUE)) 
      uw_error_cum = abs(tail(error_df$cum_uw_error, n = 1))
    } else {
      uw_error_agg = NA
      #w_error_agg = NA 
      uw_error_cum = NA
      
    }
    
  } else { #no observed data
    #print(paste0("Day ", day, ": Data is not observed, error df is not being calculated"))
    
    uw_error_agg = NA
    uw_error_cum = NA 
    error_df = NA
  }
  
  #print(paste0("Scenario ", scenario, ", ship type: ", type, " initialization ", initial))
  
  #avg cont prob
  eff_cont_main_pass= mean(contagion_prob_running_main_pass, na.rm = TRUE)
  eff_cont_main_crew = mean(contagion_prob_running_main_crew , na.rm = TRUE)
  eff_cont_act_pass = mean(contagion_prob_running_act_pass, na.rm = TRUE)
  eff_cont_act_kids = mean(contagion_prob_running_act_kids, na.rm = TRUE)
  eff_cont_emp_crew = mean(contagion_prob_running_emp_crew, na.rm = TRUE)
  eff_cont_cp = mean(contagion_prob_running_cp, na.rm = TRUE)
  
  eff_contagion_prob <- c(eff_cont_main_pass, eff_cont_main_crew, 
                          eff_cont_act_pass , eff_cont_act_kids,
                          eff_cont_emp_crew ,  eff_cont_cp)
  
  #calculate ages of exposed/infected individuals
  
  if (detailed_output == TRUE) {
    
    
    
    exposed_pass <- data.frame(indhistory_pass) %>%
      rename(day = X1, ID = X2, state = X3) %>%
      filter(state == 2) %>%
      distinct(ID) %>% pull(ID)
    
    exposed_pass_age <- table(net$passAge[exposed_pass])
    
    asymp_pass <- data.frame(indhistory_pass) %>%
      rename(day = X1, ID = X2, state = X3) %>%
      filter(state == 3) %>%
      distinct(ID) %>% pull(ID)
    
    asymp_pass_age <- table(net$passAge[asymp_pass])
    
    symp_pass <- data.frame(indhistory_pass) %>%
      rename(day = X1, ID = X2, state = X3) %>%
      filter(state == 5) %>%
      distinct(ID) %>% pull(ID)
    
    symp_pass_age <- table(net$passAge[symp_pass])
    
    
  }
  
  if (detailed_output == TRUE) {
    return(list(summaryhistory_crew = summaryhistory_crew,
                summaryhistory_pass = summaryhistory_pass,
                indhistory_crew = indhistory_crew_summary,
                indhistory_pass = indhistory_pass_summary,
                summaryhistory_direction = summaryhistory_direction, 
                summaryhistory_network = summaryhistory_network,
                exposed_pass_age = exposed_pass_age,
                asymp_pass_age = asymp_pass_age,
                symp_pass_age = symp_pass_age,
                error_df = error_df,
                uw_error_agg = uw_error_agg,
                uw_error_cum =  uw_error_cum,
                stop_message = stop_condition,
                mean_eff_cont_prob =   eff_contagion_prob,
                number_tests = number_tests,
                operations_track = operations_track
                
    ))
  } else {
    return(list(summaryhistory_crew = summaryhistory_crew,
                summaryhistory_pass = summaryhistory_pass,
                indhistory_crew = indhistory_crew_summary,
                indhistory_pass = indhistory_pass_summary,
                error_df = error_df,
                uw_error_agg = uw_error_agg,
                uw_error_cum =  uw_error_cum,
                stop_message = stop_condition,
                mean_eff_cont_prob =   eff_contagion_prob,
                number_tests = number_tests,
                operations_track = operations_track
    ))
    
  }
}



### 3. simulate_voyage_onerep_random: run simulation with different seed every time

simulate_voyage_onerep_random <- function(net, cs_input, scenario = 1, seed = NULL, groupseed = NA, covid_input, network_input, plotNetwork = T) {
  
  #create initial seed

  initial <- simulation_initialize(net = net, cs_input = cs_input, scenario = scenario, seed = seed , groupseed = groupseed,  plotNetwork = plotNetwork)
  
  simulation <- simulate_voyage_onerep(net = net, initial = initial, cs_input = cs_input, scenario = scenario,
                                       covid_input = covid_input, network_input = network_input, plotNetwork = F)
  
  return(list(summaryhistory_crew = simulation$summaryhistory_crew,
              summaryhistory_pass = simulation$summaryhistory_pass,
              indhistory_crew = simulation$indhistory_crew,
              indhistory_pass = simulation$indhistory_pass,
              summaryhistory_direction = simulation$summaryhistory_direction, 
              summaryhistory_network = simulation$summaryhistory_network,
              summaryhistory_direction = simulation$summaryhistory_direction))
}




