########################################################
################## Defining Functions ##################
##################### Do Not Alter #####################
########################################################


#Used to get number of states 
#If want to expand state space increase length of return statement
#i.e. return(c(0:2)) will be 3 states
numOfStates <- function(data){
  return(c(0:1))
}

#Used in forward-backward
#Indicator function, 0 if x and y are not equal, 1 otherwise 
Classification <- function(x_val, y_val){
  if (is.na(x_val)){return(1)}
  if (x_val == y_val){return(1)}
  else {return(0)}
}

#Normalizes matrix so that rows sum to 1
Normalize <- function(data){
  for (i in 1:dim(data)[1]){
    if (sum(data[i,],na.rm = TRUE) != 0){
      data[i,] <- data[i,]/sum(data[i,])
    }
  }
  return(data)
}

#Input: data matrix
#Output: all unique patterns, not in array/matrix form
GetPatterns <- function(states){
  library(plyr)
  options(scipen = 999)
  vals <- rep(NaN,dim(states)[1])
  for (i in 1:dim(states)[1]){
    val <- ""
    for(j in 1:dim(states)[2]){
      if (!is.na(states[i,j])){
        states_val <- states[i,j] 
      } else {
        states_val <- 2 
      }
      val <- paste0(val, states_val)
    }
    
    vals[i] <- val
  }
  return(vals)
}

#Input: Output of GetPatterns
#Output: Matrix of all unique patterns, frequency vector
#ith entry of freq_vec is # of times ith row of state_mat occurs in the original matrix 
Pattern2Data <- function(unique_patterns){
  time_length <- nchar(as.character(unique_patterns$all_patterns[[1]]))
  n <- dim(unique_patterns)[1]
  state_mat <- matrix(NaN,n,time_length)
  freq_vec <- numeric(n)
  for (i in 1:n){
    freq_vec[i] <- unique_patterns$Freq[[i]]
    pattern <- unique_patterns$all_patterns[[i]]
    
    for (time in 1:time_length){
      new_val <- as.integer(substr(pattern, time,time))
      if (new_val == 2){
        new_val <- NaN
      }
      state_mat[i,time] <- new_val
    }
  }
  return(list(state_mat, freq_vec))
}


#Forward algorithm
#Calculates P(X_1, ..., X_time, Y_time)
ForwardIter <- function(data,time,individual,initial_probabilities,transition){
  k <- numOfStates(data)
  alpha_matrix<-vector("list",time)
  alpha_i <- numeric(length(k))
  
  for (i in 1:time){
    alpha_matrix[[i]] <- alpha_i
  }
  
  for (i in 1:length(k)){
    if (!is.na(data[individual,1])){
      alpha_matrix[[1]][i] <- initial_probabilities[i] * Classification(data[individual,1], k[i])
    } else {
      alpha_matrix[[1]][i] <- initial_probabilities[i]
    }
  }
  if (time == 1){
    return(alpha_matrix)
  }
  
  for (i in 2:time){
    for (j in 1:length(k)){
      if (!is.na(data[individual,i])){
        alpha_matrix[[i]][j] <- (alpha_matrix[[i-1]] %*% transition[,j]) * Classification(data[individual,i], k[j])
      } else {
        alpha_matrix[[i]][j] <- (alpha_matrix[[i-1]] %*% transition[,j])
      }
    }
  }
  return(alpha_matrix)
}

#Backward algorithm 
#Calculates P(X_time+1, ..., X_n | Y_time)
BackwardIter <- function(data,time,individual,transition){
  k <- numOfStates(data)
  length <- dim(data)[2]
  beta_i <- numeric(length(k))
  misclass_vector <- numeric(length(k))
  beta_matrix<-vector("list",length)
  
  for (i in 1:length){
    beta_matrix[[i]] <- beta_i
  }
  
  for (i in 1:length(k)){
    beta_matrix[[length]][i] <- 1
  }
  if (time == length){
    return(beta_matrix)
  }
  for (i in (length-1):time){
    if (!is.na(data[individual,i+1])){
      for (l in 1:length(k)){
        misclass_vector[l] <- Classification(data[individual,i+1],k[l])
      }
      for (j in 1:length(k)){
        beta_matrix[[i]][j] <- sum(beta_matrix[[i+1]] * transition[j,] * misclass_vector,na.rm = TRUE)
      }
    } else {
      for (j in 1:length(k)){
        if (sum(beta_matrix[[i+1]]) == length(k)){
          beta_matrix[[i]][j] <- 1
        } else {
          beta_matrix[[i]][j] <- sum(beta_matrix[[i+1]] * transition[j,]) 
        }
      }
    }
  }
  
  return(beta_matrix)
}

#Calculates likelihood of data given all parameters
CalcLikelihoodPattern <- function(data, initial_probabilities,transition, freq_vec, pi_0, pi_1){
  likelihood_sum <- 0
  for (i in 1:dim(data)[1]){
    likelihood <- log(CalcIndLikelihood(data,i,initial_probabilities,transition,pi_0,pi_1)) * freq_vec[i]
    likelihood_sum <- likelihood_sum + likelihood
  }
  return(likelihood_sum)
}

#Calculates likelihood of an individual given all parameters
CalcIndLikelihood <- function(data, individual, initial_probabilities,transition, pi_0, pi_1){
  
  forward <- ForwardIter(data,dim(data)[2],individual,initial_probabilities,transition)
  forward_sum <- sum(forward[[dim(data)[2]]])
  
  if (max(data[individual,], na.rm = T) == 0){
    likelihood <- pi_0 + ((1 - pi_0 - pi_1) * forward_sum)
  } else if (min(data[individual,], na.rm = T) == 1){
    likelihood <- pi_1 + ((1 - pi_1 - pi_0) * forward_sum)
  } else {
    likelihood <- ((1 - pi_0 - pi_1) * forward_sum)
  }
  
  
  
  return(likelihood)
}

#Calculates probability of pi_0, pi_1 (probability of staying in the 0 and 1 state)
CalcStayer <- function(data,init,tran,freq_vec,pi_0, pi_1){
  stayer_vec_0 <- numeric(dim(data)[1])
  stayer_vec_1 <- numeric(dim(data)[1])
  for (i in 1:dim(data)[1]){
    if (max(data[i,], na.rm = T) == 0){
      stayer_vec_0[i] <- (pi_0 / CalcIndLikelihood(data,i,init,tran,pi_0, pi_1))
    } else if (min(data[i,], na.rm = T) == 1){
      stayer_vec_1[i] <- (pi_1 / CalcIndLikelihood(data,i,init,tran,pi_0, pi_1))
    }
    
  }
  
  pi_0 <- (stayer_vec_0 %*% freq_vec)/sum(freq_vec)
  pi_1 <- (stayer_vec_1 %*% freq_vec)/sum(freq_vec)
  return(c(pi_0[1],pi_1[1]))
}

#Calculates initial probability parameter (r_0, r_1)
CalcInitialPattern <- function(data, init, tran, freq_vec,pi_0,pi_1){
  k <- numOfStates(data)
  prob_list <- numeric(length(k))
  
  for (i in 1:dim(data)[1]){
    forward <- ForwardIter(data,1,i,init, tran)
    
    backward <- BackwardIter(data,1,i,tran)
    
    mover <- forward[[1]] * backward[[1]] * (1 - pi_0 - pi_1)  / CalcIndLikelihood(data,i,init,tran,pi_0, pi_1)
    
    
    for (j in k){
      prob_list[j+1] <- prob_list[j+1] + (mover[j+1] * freq_vec[i])
    }
  }
  return(prob_list/sum(prob_list))
} 

#Calculates transition probability parameter (p_00,p_01,p_10,p_11 )
CalcTransitionPattern <- function(data, init, tran,freq_vec, pi_0, pi_1){
  k <- numOfStates(data)
  tran_matrix <- matrix(0L, nrow = length(k), ncol = length(k))
  
  for (ind in 1:dim(data)[1]){
    for (time in 2:dim(data)[2]){
      forward <- ForwardIter(data, time - 1, ind, init, tran)
      backward <- BackwardIter(data, time, ind, tran)
      denom <- CalcIndLikelihood(data,ind,init,tran, pi_0,pi_1)
      for (initial_state in 1:length(k)){
        for(new_state in 1:length(k)){
          
          num <- forward[[time-1]][[initial_state]] * backward[[time]][[new_state]] * tran[initial_state,new_state] 
          num <- num * Classification(data[ind,time],new_state-1) * (1 - pi_0 - pi_1)
          
          if(denom == 0){
            frac <- 0
          } else {
            frac <- (num/denom)
          }
          tran_matrix[initial_state,new_state] <- tran_matrix[initial_state,new_state] + (frac * freq_vec[ind])
        }
      }
    }
  }
  return(Normalize(tran_matrix))
}

#Decodes the NA values for data matrix using parameters
LocalDecode <- function(data,init,tran,pi_0,pi_1){
  local_mat <- matrix(0, dim(data)[1], dim(data)[2])
  for (i in 1:dim(data)[1]){
    local_mat[i,] <- LocalDecodeInd(data,i,init,tran,pi_0,pi_1)
  }
  return(local_mat)
}


#Decodes the NA values for individual using parameters
#Compares whether 0 or 1 is more likely at each time point
#Chooses the higher probablity and imputes the value
LocalDecodeInd <- function(data, ind, init, tran, pi_0, pi_1){
  new_ind_data <- data[ind,] 
  forward_full <- ForwardIter(data,dim(data)[2], ind, init, tran)
  backward_full <- BackwardIter(data,1,ind,tran)
  
  if (min(data[ind,], na.rm = T) == 1){
    pi_0_current <- 0
    pi_1_current <- pi_1
  } else if (max(data[ind,], na.rm = T) == 0){
    pi_0_current <- pi_0
    pi_1_current <- 0
  } else {
    pi_0_current <- 0
    pi_1_current <- 0
  }
  
  for (i in 1:length(data[ind,])){
    denom <- pi_0_current + pi_1_current + ((1 - pi_0_current - pi_1_current) * (forward_full[[i]] %*% backward_full[[i]]))
    zero_prob <- (((1 - pi_0_current - pi_1_current) * (forward_full[[i]][1] * backward_full[[i]][1])) + pi_0_current) /denom
    one_prob <- (((1 - pi_0_current - pi_1_current) * (forward_full[[i]][2] * backward_full[[i]][2])) + pi_1_current) /denom
    
    if (zero_prob > one_prob){
      new_ind_data[i] <- 0
    } else {
      new_ind_data[i] <- 1
    }
    
  }
  
  return(new_ind_data)
}

DataSanitation <- function(initial_data){
  #Finds patterns in data
  #Outputs unique rows and a vector of frequency 
  all_patterns <- as.data.frame(GetPatterns(initial_data))
  unique_patterns <- as.data.frame(table(all_patterns))
  
  pattern_list <- Pattern2Data(unique_patterns)
  state_pattern <- pattern_list[[1]]
  freq_vec <- pattern_list[[2]]
  
  return(list(state_pattern, freq_vec))
}

##################### EM Algorithm #####################
EM <- function(state_pattern, freq_vec, r_0, p_00, p_11, pi_0, pi_1, epsilon){
  init <- numeric(2)
  tran <- matrix(0,2,2)
  
  init[1] <- r_0
  
  tran[1,1] <- p_00
  tran[2,2] <- p_11
  
  #Fills in rest of initial parameters
  tran[1,2] <- 1 - tran[1,1]
  tran[2,1] <- 1 - tran[2,2]
  init[2] <- 1 - init[1]
  
  
  #Running time depends on initial parameters
  #With these initial parameters should take ~1 day
  #Calculates likelihood before M step
  old_likelihood <- CalcLikelihoodPattern(state_pattern,init,tran,freq_vec,pi_0,pi_1)
  
  #maximizes each parameter
  init <- CalcInitialPattern(state_pattern, init, tran, freq_vec, pi_0, pi_1)
  tran <- CalcTransitionPattern(state_pattern,init,tran,freq_vec,pi_0, pi_1)
  stayers <- CalcStayer(state_pattern, init, tran, freq_vec, pi_0, pi_1)
  pi_0 <- stayers[1]
  pi_1 <- stayers[2]
  
  #Calculates likelihood after M step
  new_likelihood <- CalcLikelihoodPattern(state_pattern,init,tran,freq_vec,pi_0, pi_1)
  
  #Prints increase in likelihood
  #Should never be negative
  print(new_likelihood - old_likelihood)
  
  #Same as before but in a while loop
  #Continues until increases between iterations are less than epsilon
  while ((new_likelihood - old_likelihood) > epsilon){
    old_likelihood <- new_likelihood
    
    init <- CalcInitialPattern(state_pattern, init, tran, freq_vec, pi_0, pi_1)
    tran <- CalcTransitionPattern(state_pattern,init,tran,freq_vec,pi_0, pi_1)
    stayers <- CalcStayer(state_pattern, init, tran, freq_vec, pi_0, pi_1)
    pi_0 <- stayers[1]
    pi_1 <- stayers[2]
    new_likelihood <- CalcLikelihoodPattern(state_pattern,init,tran,freq_vec,pi_0, pi_1)
    
    print(new_likelihood - old_likelihood)
  }
  
  return(list(init, tran, pi_0, pi_1))
}

Decode <- function(initial_data, state_pattern, init, tran, pi_0, pi_1){
  #Decodes data using output parameters from EM
  locally_decoded <- LocalDecode(state_pattern,init,tran,pi_0,pi_1)
  #library only used for row.match
  library(prodlim)
  
  #Goes back to original data size (no pattern)
  #Takes ~30 minutes
  decoded_data <- initial_data
  for (ind in 1:dim(decoded_data)[1]){
    true_data <- decoded_data[ind,]
    true_ind <- row.match(true_data,state_pattern)
    decoded_data[ind,] <- locally_decoded[true_ind,]
  }
  return(decoded_data)
}


#Generates data
#n is number of individuals
#t is number of time points
#p is the proportion of stayers in the zero state 
GenerateAll <- function(r_0, p_00, p_11, n,t,p0,p1){
  
  init_initial <- c(r_0, 1 - r_0)
  tran_initial <- matrix(0,2,2)
  tran_initial[1,1] <- p_00
  tran_initial[1,2] <- 1 - p_00
  tran_initial[2,2] <- p_11
  tran_initial[2,1] <- 1 - p_11
  
  #Stayer_vec_true is a vector where a 0/1 at the ith index says that individual i is either a mover/stayer
  stayer_vec_true <- GenerateMS(n,p0,p1)
  #data_mat will be the created data matrix
  data_mat <- matrix(0L, nrow = n, ncol = t)
  
  #Creates only the initial values for the data
  for (i in 1:n){
    #only enters this if statement if mover
    #As data was initialized to all 0's dont need to change stayer data
    if (stayer_vec_true[i] == 0){
      #Have to add a -1 at the end so we end up with a state space of 0/1 rather than 1/2
      data_mat[i,1] <- which(rmultinom(1,1,init_initial) == 1) - 1
    } else if (stayer_vec_true[i] == 2){
      data_mat[i,1] <- 1
    }
  }
  
  #Creates the rest of the data
  for (i in 1:n){
    #starts at 2 since we already created the initial values
    for (j in 2:t){
      #similar as above, only enter if mover
      if(stayer_vec_true[i] == 0){
        #P_current is the row of the transition matrix corresponding to the state at time j-1
        p_current <- tran_initial[data_mat[i,j-1]+1,]
        
        #creates the next state at time j depending only on the state at time j-1
        data_mat[i,j] <- which(rmultinom(1,1,p_current) == 1) - 1
      } else if (stayer_vec_true[i] == 2){
        data_mat[i,j] <- 1
      }
    }
  }
  return(list(data_mat, stayer_vec_true))
}



#Generates a vector of 0/1
#A 1 at the ith index says that the ith individual is a stayer in the zero state 
#n and p are the same as above
#Used only in the function GenerateAll
GenerateMS <- function(n,p0,p1){
  
  
  stayer_vec <- numeric(n)
  for (i in 1:n){
    rand <- runif(1)
    if (rand < p0){
      stayer_vec[i] <- 1
    } else if (rand < (p0+p1)){
      stayer_vec[i] <- 2
    } 
  }
  return(stayer_vec)
}


InitialProbPattern <- function(data,stayer_vec){
  
  k <- numOfStates(data)
  init_prob <- numeric(length(k))
  
  for (i in 1:dim(data)[1]){
    if (!is.na(data[i,1])){
      if (stayer_vec[i] == 0){
        init_prob[data[i,1] + 1] <- (init_prob[data[i,1]+1] + 1)
      }
    }
  }
  
  return(init_prob/sum(init_prob))
}



#Creates a transition matrix from a data matrix
#idea is to for loop thru both dimensions of the matrix
#The transition matrix at (i,j) is increased by one every time the data transitions from i to j
#Normalized sum of rows to 1 at the end
TransitionProbPattern <- function(data, stayer_vec){
  
  k <- numOfStates(data)
  tran <- matrix(0L, nrow = length(k), ncol = length(k))
  
  for (i in 1:dim(data)[1]){
    for (j in 1:(dim(data)[2] - 1)){
      if (stayer_vec[i] == 0){
        if (data[i,j] == 0){
          if (data[i,j+1] == 0){
            tran[1,1] <- tran[1,1] +1
          } else {
            tran[1,2] <- tran[1,2] + 1
          }
        } else {
          if (data[i,j+1] == 0){
            tran[2,1] <- tran[2,1] + 1
          } else {
            tran[2,2] <- tran[2,2] + 1
          }
        }
      }
    }
  }
  
  return(Normalize(tran))
}


#Work on sometimes just missing data entirely
IntroduceMissing <- function(data){
  for (i in 1:dim(data)[1]){
    old_val <- which(rmultinom(1,1,c(.2,.2,.2,.2,.2)) == 1)
    if (old_val > 1){data[i,1:old_val-1] <- NaN}
    if (rbinom(1,1,.1)){data[i,old_val] <- NaN}
    old_val <- old_val + 1
    new_val <- 0
    while (new_val < dim(data)[2] & (old_val) < (dim(data)[2] - 2)){
      new_val <- old_val + 5
      if (new_val > dim(data)[2]){new_val <- dim(data)[2]}
      data[i,old_val:new_val] <- NaN
      if (rbinom(1,1,.1) & new_val + 1 < dim(data)[2]){data[i,new_val+1] <- NaN}
      old_val <- new_val + 2
    }
    data[i,(dim(data)[2] - 1):dim(data)[2]] <- NaN
    
  }
  return(data)
}

GetData <- function(state_data, year_data) {
  dates <- as.matrix(state_data[,6:9])
  states <- as.matrix(state_data[,10:13])
  diagnosed <- as.vector(state_data[,15])
  
  #min should be 83
  min_date <- 83
  #max is 14 but extrapolating to 16
  max_date <- 16
  
  GetDate <- function(date){
    if (date < 20){return(date+2000)} 
    else {return(date+1900)} 
  }
  
  catch_data <- matrix(NaN, dim(dates)[1], (GetDate(max_date) - GetDate(min_date) + 1))
  colnames(catch_data) <- c(GetDate(min_date):GetDate(max_date))
  
  for (i in 1:dim(dates)[1]){
    for (j in 1:dim(dates)[2]){
      if (!is.na(dates[[i,j]])){
        state_bool <- 0 
        year <- GetDate(as.integer(substr(dates[[i,j]],(nchar(dates[[i,j]]) - 1),nchar(dates[[i,j]]))))
        
        if (!is.na(states[[i,j]])){
          if (states[[i,j]] %in% year_data$state){
            row_vals <- year_data[year_data$state == states[[i,j]],]
            min_val <- row_vals[1,2]
            max_val <- row_vals[1,3]
            if (year > min_val & year < max_val){
              state_bool <- 2
              #make the following state_bool 2 if want to categorize by year
            } else {state_bool <- 2}
            #make the following state_bool 2 if want to categorize by state
          } else {state_bool <- 1}
        }
        
        if (state_bool == 2){
          adj_year <- year - GetDate(min_date) + 1
          catch_data[i,adj_year] <- 1
        }
        
        if (state_bool == 1){
          adj_year <- year - GetDate(min_date) + 1
          catch_data[i,adj_year] <- 0
        }
        
      }
    }
    if (!is.na(diagnosed[i])){
      diagnosed_date <- GetDate(as.integer(substr(diagnosed[i],(nchar(diagnosed[i]) - 1),nchar(diagnosed[i]))))
      adj_diagnosed_date <- diagnosed_date - GetDate(min_date) + 1
      catch_data[i,adj_diagnosed_date] <- 1
    }
    
    
    
    
  }
  return(catch_data)
}


########################################################
#################### Function Calls ####################
########################  Alter ########################
########################################################

########################################################
##################### INTRUCTIONS ######################
# If using real data uncomment lines 542-545, 606 and 614 comment lines 550-577, 607, and 615
# If using sumlated data uncomment lines 550-577, 607, and 615 and comment lines 542-545, 606 and 614

########################################################
################### Data Generation ####################
########################################################

###################### Real Data #######################

# #Loads in initial_data matrix, data matrix used in this example
# #If using different matrix dont use this line and change parameter in DataSanitation and Decode
# 
state_data <- read.delim("State_residence_May092019.txt")
year_data <- read.csv("year registry state.csv")
initial_data <- GetData(state_data, year_data)
rownames(initial_data) <- state_data[,1]
full_data <- DataSanitation(initial_data)
state_pattern <- full_data[[1]]
freq_vec <- full_data[[2]]

################### Simulated Data #####################
# #n is the number of individuals and t is the number of time points for the simulated data
# n <- 1000
# t <- 33
# 
# #Initial parameters for simulation
# r_0_sim <- .40
# p_00_sim <- .96
# p_11_sim <- .98
# pi_0_sim <- .22
# pi_1_sim <- .51
# 
# #Generates simulated data
# all_data_sim <- GenerateAll(r_0_sim, p_00_sim, p_11_sim,n,t,pi_0_sim,pi_1_sim)
# sim_data <- all_data_sim[[1]]
# stayer_vec_sim_true <- all_data_sim[[2]]
# 
# #True simulated parameters
# init_sim_true <- InitialProbPattern(sim_data,stayer_vec_sim_true)
# tran_sim_true <- TransitionProbPattern(sim_data,stayer_vec_sim_true)
# pi_0_sim_true <- sum(stayer_vec_sim_true == 1)/n
# pi_1_sim_true <- sum(stayer_vec_sim_true == 2)/n
# 
# #removing observations
# sim_data_obs <- IntroduceMissing(sim_data)
# 
# #Simulated observed data
# full_data_sim <- DataSanitation(sim_data_obs)
# state_pattern_sim <- full_data_sim[[1]]
# freq_vec_sim <- full_data_sim[[2]]

########################################################
#################### EM Algorithm ######################
########################################################


#Inputs: initial_data: data matrix to use, change name if using different matrix
#       r_0: probability of starting in the 0 state
#       p_00: probability of transitioning from 0 to 0
#       p_11: probability of transitioning from 1 to 1
#       pi_0: Probability of being a stayer in the 0 state
#       pi_1: Probability of being a stayer in the 1 state
#       epsilon: If likelihood increases between EM iterations are less than this value EM is considered converged

r_0 <- .4036024
p_00 <- .9696433
p_11 <- .97841330
pi_0 <- .2239507
pi_1 <- .5127309
epsilon <- .000001

#Outputs: state_pattern: unique patterns of initial_data
#         freq_vec: Frequency vector to go along with state_pattern
#                   ith entry corresponds to how many times ith row in state_pattern was observed
#         parameters: initial, transition, pi_0, and pi_1 estimates from EM
#         decoded_data: Imputed initial_data using EM estimates

#Only use one of the following 2 lines, first is for real data, second is for simulated data
###################### Real Data #######################
parameters <- EM(state_pattern,freq_vec, r_0, p_00, p_11, pi_0, pi_1, epsilon)
################### Simulated Data #####################
# parameters <- EM(state_pattern_sim,freq_vec_sim, r_0, p_00, p_11, pi_0, pi_1, epsilon)

# save(parameters, file = "parameters.rda")

########################################################
###################### Decoding ########################
########################################################

#Only use one of the following 2 lines, first is for real data, second is for simulated data
###################### Real Data #######################
decoded_data <- Decode(initial_data, state_pattern, parameters[[1]], parameters[[2]], parameters[[3]], parameters[[4]])
################### Simulated Data #####################
# decoded_data <- Decode(sim_data_obs, state_pattern_sim, parameters[[1]], parameters[[2]], parameters[[3]], parameters[[4]])



##################### Output csv #######################
# #Following line is only for real data
row.names(decoded_data) <- state_data[,1]
write.csv(decoded_data, file = "decoded_data.csv")


