library(ggplot2)
library(tidyverse)
library(dplyr)
library(hrbrthemes)
library(viridisLite)
library(viridis)
library(ggridges)

# Note: simulation isolated in sims.R
Qpwise <- function(iterations,initcohorttype,initcohortsize,transmissionrate,recoveryrate,arrivalrate,departrate,maxtimeval,maxprox){
  
  # PARAMETERS:
  # iterations - number of overall simulations performed
  # initcohorttype - {0,1} 0 for all initial infected at time zero, 1 for all initial infected at time unif dist(0,1)
  # initcohortsize - total number of initial infected
  # transmissionrate - beta 
  # recoveryrate - gamma
  # arrivalrate - rate at which individuals come within proximity for queue process
  # departrate - rate at which individuals leave proximity for queue process
  # maxtimeeval - max time allowed for queue process
  # maxprox - max number of individuals initially in proximity to uniformally draw from
  
  genbysim <- data.frame(generation = numeric(),
                         simnumb = character())
  
  # Initialize Simulation .............................................................................................
  for (k in 1:iterations){
    # Generate a number of infected individuals and the times at which they were infected  
    # Refer to this as the initial cohort
    if (initcohorttype == 0){
      vector <- rep(0,initcohortsize) # all initial cohort infections occur at time=0
    } else if (initcohorttype == 1){
      vector <- runif(initcohortsize, min=0, max=1) # all initial cohort infections occur between 0 and 1 uniformally
    } else {
      print("ERROR in initcohorttype")
    }
    vector <- sort(vector) # put vector of initial cohort into chronological order 
    
    # Create data frame this goes into
    # contactdata <- data.frame(icmnumber = numeric(),
    #                           icmtime = numeric(),
    #                           secondinfecttime = numeric(),
    #                           genint = numeric())
    contactdata <- data.frame()
    # ...................................................................................................................
    # Execute the contact regime for one generation of secondary infections to be added to data frame
    
    for (i in 1:length(vector)){ # for each init cohort in vector run queue contact simulation
      
      # QUEUE ---------------------------------------------------------------------------------------------------------------------
      # function which takes in arrival and departure rates with an initial population
      # Initiate
      queuedata <- data.frame()
      timein <- 0
      timeout <- 0
      timespent <- 0
      inittime <- vector[i]
      # set init prox for each person
      initprox <- floor(runif(1,min=0,max=maxprox))
      
      ###print(paste0("There are ", initprox, " person(s) initially in proximity to initcohort member ", i, " ." ))
      
      # Initial Proximity **************************************************************************************************************
      if (initprox == 0){ # if no one initially in proximity skip
        ###print(paste0("No one is included in the initial proximity."))
      } else if (initprox > 0) { # if initially in proximity then...
        for (x in 1:initprox){
          arrives = 0 + inittime # initprox all arrive at time = 0
          if (departrate == 0){ # if initprox can't leave set long departure
            departs = maxtimeval + inittime
            ###print(paste0("Person from initprox ", x, " arrived at time ", arrives, " and did NOT depart until ", departs, "."))
          } else if (departrate > 0){ # if initprox can leave then...
            departs = rexp(1, rate = 1/departrate) + inittime # exp dist departures by departrate (since others can leave too -> could eventually impose different rate for init vs pop)
            ###print(paste0("Person from initprox ", x, " arrived at time ", arrives, " and departed at time ", departs, " (total time later)."))
          } else {
            stop("ERROR in departrate, must be >= 0")
          }
          # append results of initprox to queue
          initoutcome = c(x,arrives,arrives,departs, departs)
          queuedata <- rbind(queuedata, initoutcome)
        }
      } else {
        stop("ERROR in initprox, must be >= 0")
      }
      # ********************************************************************************************************************************
      
      # Dynamic Proximity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (arrivalrate == 0){ # if no one enters then skip (accounted for above with initprox)
        ###print(paste0("No one from the population can enter proximity when arrivalrate is ", arrivalrate,"."))
        if (initprox == 0){
          stop("ERROR in arrivalrate = ", arrivalrate, " and initprox = ", initprox, " -> no queue produced because one can interact!")
        }
      } else if (arrivalrate > 0){ # if pop can enter then...
        # for each pop member see when they arrive and depart
        
        index = initprox + 1
        while (timein < maxtimeval){
          arrival = rexp(1, rate = 1/arrivalrate) + inittime # define arrival exp distribution using arrivalrate
          if (timein + arrival > maxtimeval){
            ###print(paste0("ERROR in arrival, exceeds maxtimeval!"))
            break
          }
          if (departrate == 0){ # if pop can't leave (but can enter) then define same late departure for all
            departure = maxtimeval + inittime
            ###print(paste0("Overall person ", index, " arrived ", arrival, " later and did NOT depart ", departure, "."))
            # adjust time counts for departrate = 0
            timein <- timein + arrival 
            timeout <- departure
            if (timeout > maxtimeval){
              timeout <- maxtimeval
            }
            timespent <- departure - arrival
          } else if (departrate > 0){ # if pop can leave then...
            departure = rexp(1, rate = 1/departrate) + inittime # define departure distribution using departrate
            ###print(paste0("Overall person ", index, " arrived ", arrival, " later and departed ", departure, " later."))
            # adjust time counts for departrate > 0
            timein <- timein + arrival 
            timeout <- timein + departure
            if (timeout > maxtimeval){
              timeout <- maxtimeval
            }
            timespent <- departure
          } else {
            stop("ERROR in departrate")
          }
          ###print(paste0("Overall person number ", index, " arrived ", arrival, " later, which is at time ", timein, " and departed ", departure, " later, which is at time ", timeout, "."))
          if (timein > timeout){
            print(paste0("ERROR in last arrival, past maxtimeval!"))
            break
          }
          # append outcomes to df
          outcome = c(index,timein,arrival,timeout,timespent)
          queuedata <- rbind(queuedata, outcome)
          index <- index + 1
        }
      } else {
        stop("ERROR in arrivalrate, must be >= 0")
      }
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      names(queuedata) <- c('newpersonid','arrivaltime','interarrival','departuretime','timespent')
      
      # Continuous Cumulative Proximity @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      # compute number of individuals within proximity for all time
      numindiv <- data.frame()
      # Append 1 to arrivals to represent addition to pop
      addarrivals <- cbind(queuedata$arrivaltime, rep(1, nrow(queuedata)))
      numindiv <- rbind(numindiv,addarrivals)
      # Append -1 to departures to represent removal from pop
      adddeparts <- cbind(queuedata$departuretime, rep(-1, nrow(queuedata)))
      numindiv <- rbind(numindiv,adddeparts)
      names(numindiv) <- c('timing', 'arrivORdept')
      
      # order df to compute arrivals/departures in order they occur
      numindiv <- numindiv[order(as.numeric(numindiv$timing),  decreasing = FALSE),]
      countsus <- 0 
      sus_count <- numeric(nrow(numindiv))
      
      for (y in 1:nrow(numindiv)){
        if (numindiv$arrivORdept[y] == 1){
          countsus <- countsus + 1
          ##print(paste0("The ", j, " th change resulted in ", countsus, " susceptibles after an arrival."))
        } else if (numindiv$arrivORdept[y] == -1){
          countsus <- countsus - 1
          ##print(paste0("The ", j, " th change resulted in ", countsus, " susceptibles after a departure."))
        } else {
          stop("ERROR in count!!")
        }
        sus_count[y] <- countsus
      }
      numindiv$sus_count <- sus_count
      # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
      # END QUEUE  ----------------------------------------------------------------------------------------------------------------------------
      
      # DRAW + APPLY IIP (PWISE) ==============================================================================================================
      
      tau1 = rexp(1,recoveryrate) + inittime # adjusted for initcohort timing exp dist of factor of IIP -> this case timing of infection
      tau2 = tau1 + rexp(1,recoveryrate)
      ###print(paste0("The IIP for init cohort individual ", i, " is determined by tau1 = ", tau1, " and tau2 = ", tau2," which is adjusted by init cohort timing."))
      
      # Iterate through Queue output to determine who overlaps with infectiousness period and when
      interacts = data.frame()
      for (n in 1:nrow(queuedata)){
        if (queuedata[n,2] < tau1){ # strictly less than for type 2 vs. 4
          if (queuedata[n,4] <= tau1){
            ###print(paste0("Queue Person ID: ", queuedata[n,1], " for initcohort member ", i, " did not overlap with the infectiousness period (type 1)."))
            next
          } else if (queuedata[n,4] < tau2){ # strictly less than for type 2 vs. 3
            ###print(paste0("Queue Person ID: ", queuedata[n,1], " for initcohort member ", i, " did overlap with the infectiousness period."))
            overlap = queuedata[n,4] - tau1
            ###print(paste0("Queue Person ID: ", queuedata[n,1], " overlapped for ", overlap, " time units (type 2)."))
            addnew = data.frame(
              newpersonid = queuedata[n, 1],
              overlap = overlap,
              overlap_start = tau1,
              overlap_end = queuedata[n, 4]
            )
            #print(addnew)
            interacts <- rbind(interacts, addnew)
            next
          } else if (queuedata[n,4] >= tau2){ 
            ###print(paste0("Queue Person ID: ", queuedata[n,1], " for initcohort member ", i, " did overlap with the infectiousness period."))
            overlap = tau2 - tau1
            ###print(paste0("Queue Person ID: ", queuedata[n,1], " overlapped for ", overlap, " time units (type 3)."))
            addnew = data.frame(
              newpersonid = queuedata[n, 1],
              overlap = overlap,
              overlap_start = tau1, 
              overlap_end = tau2  
            )
            #print(addnew)
            interacts <- rbind(interacts, addnew)
            next
          } else {
            ###print(paste0("Queue Person ID: ", queuedata[n,1], "did not meet first 3 types."))
          }
        } else if (queuedata[n,2] < tau2){ # this is strictly less than for type 6
          if (queuedata[n,4] <= tau2){
            ###print(paste0("Queue Person ID: ", queuedata[n,1], " for initcohort member ", i, " did overlap with the infectiousness period."))
            overlap = queuedata[n,4] - queuedata[n,2]
            ###print(paste0("Queue Person ID: ", queuedata[n,1], " overlapped for ", overlap, " time units (type 4)."))
            addnew = data.frame(
              newpersonid = queuedata[n, 1],
              overlap = overlap,
              overlap_start = queuedata[n,2],  
              overlap_end = queuedata[n, 4]  
            )
            #print(addnew)
            interacts <- rbind(interacts, addnew)
            next
          } else if (queuedata[n,4] > tau2){
            ###print(paste0("Queue Person ID: ", queuedata[n,1], " for initcohort member ", i, " did overlap with the infectiousness period."))
            overlap = tau2 - queuedata[n,2]
            ###print(paste0("Queue Person ID: ", queuedata[n,1], " overlapped for ", overlap, " time units (type 5)."))
            addnew = data.frame(
              newpersonid = queuedata[n, 1],
              overlap = overlap,
              overlap_start = queuedata[n,2],  
              overlap_end = tau2  
            )
            #print(addnew)
            interacts <- rbind(interacts, addnew)
            next
          }
        } else { # (queuedata[n,2] >= tau2) => (queuedata[n,4] > tau2)
          ###print(paste0("Queue Person ID: ", queuedata[n,1], " for initcohort member ", i, " did not overlap with the infectiousness period (type 6)."))
          next
        }
      }
      ###print(interacts)
      
      # PART 2 - APPLY SURVIVAL ANALYSIS TO QUEUEDATA
      # if no resulting interactions, go to next loop...
      if (nrow(interacts) == 0){
        ###print(paste0("NOTE: No infectious interactions for initcohort member ", i, " so moving on..."))
        next
      } else {
        # Poisson process over length of overlap to determine number of possible hits then draw from a uniform to determine when they are
        for (z in 1:nrow(interacts)){
          # Draw from Poisson 
          hitrate = transmissionrate*(interacts[z,4] - interacts[z,3])
          if (is.na(hitrate)){
            stop("ERROR in hitrate NA!")
          }
          hits = rpois(1,hitrate)
          ###print(paste0("Interaction ", interacts[z,1], " aquired ", hits, " hits with a rate of ", hitrate, "."))
          if (hits == 0){
            next
          } else {
            # Draw from Uniform
            timing = runif(hits,interacts[z,3],interacts[z,4])
            # First hit to occur sequentially is the timing of secondary infection
            sortedtiming = sort(timing)
            secondaryinf = sortedtiming[1]
            ###print(paste0("Interaction ", interacts[z,1]," produced a secondary infection at time ", secondaryinf, " during a ", interacts[z,2], " overlap from time ", interacts[z,3]," to ", interacts[z,4], "."))
            output = data.frame(
              icmnumber = i,
              icmtime = inittime, #vector[i],
              secondinfecttime = secondaryinf,
              genint = secondaryinf - inittime #vector[i]
            )
            contactdata <- rbind(contactdata,output)
          }
        }
      }
      ###print(contactdata)
      if (length(contactdata) == 0){
        ###print(paste0("WARNING: Initcohort ", i, " had no secondary infection contacts produced so moving on..."))
        next
      } else {
      # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      # Data Wrangling __________________________________________________________________________________________________________________________
      names(contactdata) <- c('icmnumber','icmtime','secondinfecttime','genint') # insert names in contactdata
      contactdata <- contactdata[!is.na(contactdata$genint), ] # remove NA to prevent bugs
      contactdata <- contactdata[order(contactdata$genint, decreasing = FALSE),] # order based on timing of generation interval length 
      gidata = contactdata$genint # extract the generation interval data
      namedintervals <- data.frame() # create data frame to input each simulation name to gi data that resets each sim
      if (length(gidata) == 0){
        #print(paste0("NOTE: In simulation number ", k," person ", i," did not produce any secondary infections!"))
        next
      } else {
        namedintervals <- data.frame( # append simulation group
          generation = gidata,
          simnumb = paste0("sim_", k),
          stringsAsFactors = FALSE
        ) 
      }
      }
    }
    
    namedintervals <- namedintervals[!is.na(namedintervals$generation),] # remove any NA GI rows
    # confirm types
    namedintervals$simnumb <- as.character(namedintervals$simnumb)
    namedintervals$generation <- as.numeric(namedintervals$generation)
    
    # Fix ordering issue
    namedintervals <- namedintervals %>%
      group_by(simnumb) %>%
      arrange(generation)
    
    # Add Cumulative Sum by Simulation
    namedintervals$orderofinf <- ave(namedintervals$generation,namedintervals$simnumb, FUN=seq_along) # assign number for order generation intervals occur
    
    # Normalize the Cumulative Sum
    namedintervals <- namedintervals %>%
      group_by(simnumb) %>%
      mutate(max_order = max(orderofinf, na.rm = TRUE)) %>%
      mutate(normcsum = ifelse(max_order > 0, orderofinf / max_order, NA)) %>%
      select(-max_order)
    
    # Now combine computed info for the k simulation with all sims (genbysim)
    genbysim <- rbind(genbysim,namedintervals)
    
    # Specify the Theoretical Distribution Characteristics
    x_vals <- seq(0, max(genbysim$generation), length.out = 1000)
    
    theoretical <- data.frame(
      time = x_vals,
      pdf = recoveryrate * exp(-(recoveryrate) * x_vals),
      cdf = 1 - exp(-(recoveryrate) * x_vals),
      dist = rexp(x_vals, rate = recoveryrate)
    )
    
    #  ________________________________________________________________________________________________________________________________________
    
    # Data References for Plotting
    # queuedata - queue outputs for each person (newpersonid), time at which they arrive (arrivaltime), time selected by the distribution of 
    #             arrivals (interarrival), time at which they depart (departuretime), difference btwn arrival and departure (timespent) aka 
    #             departure selected from distribution
    #             names(queuedata) <- c('newpersonid','arrivaltime','interarrival','departuretime','timespent')
    ##print(queuedata)
    
    # numindiv -  ordered increasing timing of an arrival or departure (timing) from proximity denoted by 1 or -1, respectively (arrivORdept) 
    #             to compute cumulative number of individuals within proximity at each timing (sus_count)
    #             names(numindiv) <- c('timing', 'arrivORdept') + sus_count
    ##print(numindiv)
    
    # contactdata - MAIN OUTPUT* ordered by increasing generation interval, each individual with an identifier (icmnumber) receives an 
    #               infection at time (secondinfecttime) from an individual from the initial cohort infected at time (icmtime) and 
    #               the time between these values (genint)
    #               names(contactdata) <- c('icmnumber','icmtime','secondinfecttime','genint')
    ##print(contactdata)
    
    # namedintervals - intermediate step to compute all outputs found in genbysim for each simulation iteration. All outputs of this after 
    #                  appended to genbysim
    
    # genbysim - contains all the generation intervals values (generation) with an appended number for cumulative purposes (orderofinf) 
    #            which are normalized (normcsum) per simulations performed (simnumb) from the input parameter 'iterations'
    #             names(genbysim) <- c('generation', 'simnumb') + orderofinf + normcsum
    ##print(genbysim)
    
    # Tracking Progress of Large Simulations
    print(paste0("Simulation number: ", k, " finished and produced ", nrow(namedintervals), " secondary infections from an initial cohort of ", initcohortsize, " people."))
  }
  
  # Figures
  mybinsize = 0.01
  
  # Plot distribution of generation interval length density curves of each simulation and the theoretical exp
  p1 <- ggplot(genbysim, aes(x = generation, color = simnumb)) +
    geom_density(linewidth = 1.2) +
    geom_line(data = theoretical, aes(x = time, y = pdf), color = "grey", linetype = "longdash", linewidth = 1.2) +
    geom_density(data = theoretical, aes(x = dist), color = "black", linetype = "longdash", linewidth = 1.2) +
    theme_minimal() +
    labs(
      title = "PWise Density of Generation Intervals Distribution by Simulation",
      x = "Generation Interval",
      y = "Density",
      color = "Simulation"
    )
  
  # Plot histograms of generation interval length for each simulation and the theoretical exp
  p2 <- genbysim %>%
    mutate(simnumb = fct_reorder(simnumb, generation)) %>%
    ggplot(aes(x=generation, color=simnumb, fill=simnumb)) +
    geom_histogram(alpha=0.6, binwidth = mybinsize) +
    geom_histogram(data = theoretical, aes(x = dist), fill = "black", color = NA, alpha = 0.4, binwidth = mybinsize) +
    theme_ipsum() +
    theme() +
    labs(
      title = "PWise Generation Interval Distributions by Simulation",
      x = "Generation Interval Length",
      y = "Number of Secondary Infections"
    ) +
    facet_wrap(~simnumb)
  
  # alt use normalized versions
  p3 <- genbysim %>%
    mutate(simnumb = fct_reorder(simnumb, generation)) %>%
    ggplot(aes(x = generation, color = simnumb, fill = simnumb)) +
    stat_bin(aes(y = after_stat(count / max(count))), alpha = 0.6, binwidth = mybinsize, position = "identity") + # Normalize each histogram within its own facet
    stat_bin(data = theoretical, aes(x = dist, y = after_stat(count / max(count))), fill = "black", color = NA, alpha = 0.4, binwidth = mybinsize, position = "identity") + # Normalize theoretical
    theme_ipsum() +
    labs(
      title = "PWise Normalized Generation Interval Distributions by Simulation",
      x = "Generation Interval",
      y = "Normalized Number of Secondary Infections"
    ) +
    facet_wrap(~simnumb)
  
  p4 <- genbysim %>%
    mutate(simnumb = fct_reorder(simnumb, generation)) %>%
    ggplot(aes(x = generation, y = after_stat(density * width), color = simnumb, fill = simnumb)) +
    geom_histogram(alpha=0.6, binwidth = mybinsize) +
    geom_histogram(data = theoretical, aes(x = dist, y = after_stat(density * width)), fill = "black", color = NA, alpha = 0.4, binwidth = mybinsize) +
    theme_ipsum() +
    labs(
      title = "PWise Relative Generation Interval Distributions by Simulation",
      x = "Generation Interval",
      y = "Relative Secondary Infections"
    ) +
    facet_wrap(~simnumb)
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
}
