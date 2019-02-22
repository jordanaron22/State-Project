
state_data <- read.csv("State_residence_Dec112018.csv")
year_data <- read.csv("year registry state.csv")
dates <- as.matrix(state_data[,6:9])
states <- as.matrix(state_data[,10:13])

#should be 83
min_date <- 83
max_date <- 14

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
}



GetPatterns <- function(states,stayer_mat){
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
    if (!missing(stayer_mat)){
      val <- paste0(val,stayer_mat[i])
    }
    
    vals[i] <- val
  }
  return(vals)
}



Pattern2Data <- function(unique_patterns, stayer_mat){
  time_length <- nchar(as.character(unique_patterns$all_patterns[[1]]))
  n <- dim(unique_patterns)[1]
  state_mat <- matrix(NaN,n,time_length)
  freq_vec <- numeric(n)
  stayer_mat <- numeric(n)
  for (i in 1:n){
    freq_vec[i] <- unique_patterns$Freq[[i]]
    pattern <- unique_patterns$all_patterns[[i]]
    
    if (!missing(stayer_mat)){
      #stayer_mat[i] <- as.integer(substr(pattern, years_to_impute,years_to_impute))
    }
    
    for (time in 1:time_length){
      new_val <- as.integer(substr(pattern, time,time))
      if (new_val == 2){
        new_val <- NaN
      }
      state_mat[i,time] <- new_val
    }
  }
  if (!missing(stayer_mat)){
    return(list(state_mat, freq_vec,stayer_mat))
  }
  return(list(state_mat, freq_vec))
}





all_patterns <- as.data.frame(GetPatterns(catch_data))
unique_patterns <- as.data.frame(table(all_patterns))

pattern_list <- Pattern2Data(unique_patterns)
state_pattern <- pattern_list[[1]]
freq_vec <- pattern_list[[2]]





numOfStates <- function(data){
  return(c(0:1))
}



Classification <- function(x_val, y_val){
  if (is.na(x_val)){return(1)}
  if (x_val == y_val){return(1)}
  else {return(0)}
}



Normalize <- function(data){
  for (i in 1:dim(data)[1]){
    if (sum(data[i,],na.rm = TRUE) != 0){
      data[i,] <- data[i,]/sum(data[i,])
    }
  }
  return(data)
}



GetPatterns <- function(states,stayer_mat){
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
    if (!missing(stayer_mat)){
      val <- paste0(val,stayer_mat[i])
    }
    
    vals[i] <- val
  }
  return(vals)
}



Pattern2Data <- function(unique_patterns, stayer_mat){
  time_length <- nchar(as.character(unique_patterns$all_patterns[[1]]))
  n <- dim(unique_patterns)[1]
  state_mat <- matrix(NaN,n,time_length)
  freq_vec <- numeric(n)
  stayer_mat <- numeric(n)
  for (i in 1:n){
    freq_vec[i] <- unique_patterns$Freq[[i]]
    pattern <- unique_patterns$all_patterns[[i]]
    
    if (!missing(stayer_mat)){
      #stayer_mat[i] <- as.integer(substr(pattern, years_to_impute,years_to_impute))
    }
    
    for (time in 1:time_length){
      new_val <- as.integer(substr(pattern, time,time))
      if (new_val == 2){
        new_val <- NaN
      }
      state_mat[i,time] <- new_val
    }
  }
  if (!missing(stayer_mat)){
    return(list(state_mat, freq_vec,stayer_mat))
  }
  return(list(state_mat, freq_vec))
}




#JASA paper method 
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




#JASA method 
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



CalcLikelihoodPattern <- function(data, initial_probabilities,transition, freq_vec, pi_0, pi_1){
  likelihood_sum <- 0
  for (i in 1:dim(data)[1]){
    likelihood <- log(CalcIndLikelihood(data,i,initial_probabilities,transition,pi_0,pi_1)) * freq_vec[i]
    likelihood_sum <- likelihood_sum + likelihood
  }
  return(likelihood_sum)
}




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



init <- c(.4,.6)

tran <- matrix(0,2,2)
tran[1,1] <- .53
tran[2,2] <- .39

tran[1,2] <- 1 - tran[1,1]
tran[2,1] <- 1 - tran[2,2]

pi_0 <- .25
pi_1 <- .6



old_likelihood <- CalcLikelihoodPattern(state_pattern,init,tran,freq_vec,pi_0,pi_1)


init <- CalcInitialPattern(state_pattern, init, tran, freq_vec, pi_0, pi_1)
tran <- CalcTransitionPattern(state_pattern,init,tran,freq_vec,pi_0, pi_1)
stayers <- CalcStayer(state_pattern, init, tran, freq_vec, pi_0, pi_1)
pi_0 <- stayers[1]
pi_1 <- stayers[2]
new_likelihood <- CalcLikelihoodPattern(state_pattern,init,tran,freq_vec,pi_0, pi_1)

print(new_likelihood - old_likelihood)

while ((new_likelihood - old_likelihood) > 0.00001){
  old_likelihood <- new_likelihood
  
  init <- CalcInitialPattern(state_pattern, init, tran, freq_vec, pi_0, pi_1)
  tran <- CalcTransitionPattern(state_pattern,init,tran,freq_vec,pi_0, pi_1)
  stayers <- CalcStayer(state_pattern, init, tran, freq_vec, pi_0, pi_1)
  pi_0 <- stayers[1]
  pi_1 <- stayers[2]
  new_likelihood <- CalcLikelihoodPattern(state_pattern,init,tran,freq_vec,pi_0, pi_1)
  
  print(new_likelihood - old_likelihood)
}

state_project_parameters <- list(init, tran, pi_0, pi_1)
save(state_project_parameters, file = "State Project Faster.rda")



while ((new_likelihood - old_likelihood) > 0.00000001){
  old_likelihood <- new_likelihood
  
  init <- CalcInitialPattern(state_pattern, init, tran, freq_vec, pi_0, pi_1)
  tran <- CalcTransitionPattern(state_pattern,init,tran,freq_vec,pi_0, pi_1)
  stayers <- CalcStayer(state_pattern, init, tran, freq_vec, pi_0, pi_1)
  pi_0 <- stayers[1]
  pi_1 <- stayers[2]
  new_likelihood <- CalcLikelihoodPattern(state_pattern,init,tran,freq_vec,pi_0, pi_1)
  
  print(new_likelihood - old_likelihood)
}

state_project_parameters <- list(init, tran, pi_0, pi_1)
save(state_project_parameters, file = "State Project Slower.rda")


