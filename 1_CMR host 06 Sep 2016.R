
#------ Model host system
#--- Host = Multi-state model with capture-mark-recapture analysis
        # Capture-mark-recapture data
#--- Pathogen = Hierarchical model
        # re-infection = infection when you are already infected

#-- Model simulation
        # Modeling the dynamics of a capture-mark-recapture host population

# 4 state system
  # Not entered
  # uninfected
  # infected
  # Dead

# 3 observed states
  # seen uninfected
  # seen infected
  # not seen

#------ Consider revising ------#
  # All individuals enter as uninfected


#--- Survey conditions
n.occasions <- 7

#---- Population parameters
# Super population size
N_U <- 350
N_I <- 300

N <- N_U + N_I

# Number of states
n.states <- 4

# Number of observable states
n.obs <- 3

#--- Define parameter values
# Uninfected survival probability
phi_U <- 0.9

# Infected survival probability
phi_I <- 0.7

# Entry probability     # Must sum to one
gamma_U <- (1/(n.occasions-1))/2
gamma_I <- (1/(n.occasions-1))/2

(gamma_U * (n.occasions-1))*2

# Transition probability 
#(Uninfected to Infected)
beta_UI <- 0.2

#(Infected to Uninfected)
beta_IU <- 0.3


# Detection probability
p_U <- 0.9
p_I <- 0.8

#---- 1. State process

PSI.state <- array(NA, dim = c(n.states, n.states, N, n.occasions-1))
for(i in 1:N){
  for(j in 1:n.occasions-1){
    PSI.state[,,i,j] <- matrix(c(1 - gamma_U - gamma_I, gamma_U,  gamma_I,  0,
                                 
                                 0,     phi_U * (1- beta_UI), phi_U * beta_UI, 1-phi_U,
                                 
                                 0,     phi_I * beta_IU, phi_I * (1- beta_IU), 1-phi_I,
                                 
                                 0,     0,        0,       1), 
                               nrow = 4, ncol = 4, byrow = T)
  }
}

# 2. Observational process
PSI.obs <- array(NA, dim = c(n.states, n.obs, N, n.occasions-1))
for(i in 1:N){
  for(j in 1:n.occasions-1){
    PSI.obs[,,i,j] <- matrix(c(  0,       0,   1,
                                 p_U,     0,   1- p_U,
                                 0,     p_I,   1- p_I,
                                 0,       0,   1), 
                               nrow = 4, ncol = 3, byrow = T)
  }
}


#---- Function to simulate capture-recapture data under the JS model

simul.js <- function(PSI.state, PSI.obs, gamma_U, gamma_I, N, N_U, N_I){
  
  n.occasions <- dim(PSI.state)[4] + 1
  
  B_U <- rmultinom(1, N_U, rep(gamma_U, times = (n.occasions-1)))
  B_I <- rmultinom(1, N_I, rep(gamma_I, times = (n.occasions-1)))
  # Generate no. of entering hosts per occasion
        # N is superpopulation size
  
  B <- B_U + B_I

  CH.sur <- CH.p <- matrix(0, ncol = n.occasions, nrow = N)
  
# Define a vector with the occasion of entering the population
  
  ent.occ_U <- ent.occ_I <- numeric()
  
  for (t in 1:(n.occasions-1)){
    
    ent.occ_U <- c(ent.occ_U, rep(t, B_U[t]))
    ent.occ_I <- c(ent.occ_I, rep(t, B_I[t]))
    
    }

ent.occ <- c(ent.occ_U, ent.occ_I)  

# Simulating survival
  
for (i in 1:length(ent.occ_U)){     # For each individual in the superpopulation
    CH.sur[i, ent.occ_U[i]] <- 2
    CH.p[i, ent.occ[i]] <- 1
}

for(i in (length(ent.occ_U)+1):length(ent.occ)){    
    CH.sur[i, ent.occ[i]] <- 3
    CH.p[i, ent.occ[i]] <- 2
}

for (i in 1:N){           

# Probability of captureing the individual first occasion    

# if(CH.sur[i, ent.occ[i]] == 2){ 
#   CH.p[i, ent.occ[i]] <- 1 * rbinom(1, 1, p_U)
# }else{  
#   CH.p[i, ent.occ[i]] <- 2 * rbinom(1, 1, p_I) 
#   }
    
  if (ent.occ[i] == n.occasions) next   # If the entry occasion = the last occasion, go next
    
    for(t in (ent.occ[i]+1):n.occasions){  # for each occasion from entering 2 last occasion
      
# Determine individual state history
      
      sur <- which(rmultinom(1, 1, PSI.state[CH.sur[i, t-1], , i, t-1]) == 1)
      
      CH.sur[i, t] <- sur
      
# Determine if the individual is captured

      event <- which(rmultinom(1, 1, PSI.obs[CH.sur[i, t], , i, t-1]) == 1)
      
      CH.p[i, t] <- event
      
    } #t
    
  } #i

#---- Reorder
CH.p <-
  CH.p[
    c(which(ent.occ == 1), 
      which(ent.occ == 2),
      which(ent.occ == 3), 
      which(ent.occ == 4)),]

CH.sur <-
  CH.sur[
    c(which(ent.occ == 1), 
      which(ent.occ == 2),
      which(ent.occ == 3), 
      which(ent.occ == 4)),]

# Remove individuals never captured
  
CH <- CH.p * CH.sur

cap.sum <- rowSums(CH)
never <- which(cap.sum == 0)  

if(length(never) > 0) {
CH.p <- CH.p[-never,]  
CH.sur <- CH.sur[-never,]  
}

Nt <- numeric(n.occasions)

for(i in 1:n.occasions){
  
  Nt[i] <- length(which(CH.sur[,i] > 0))

}


#-----

CH <- CH.p
CH[is.na(CH) == TRUE] <- 0
CH[CH == dim(PSI.state)[1]] <- 0
id <- numeric(0)

for(i in 1:dim(CH)[1]){     # For each individual row
  
  z <- min(which(CH[i, ] != 0))   
      # which column is the lowest that does not = 0?
  
  ifelse(z == dim(CH)[2], id <- c(id, i), id <- c(id))
      # keep a list of individuals only captured on last occasion
}
 
CH.sur[CH.sur == 0] <- 1

if(length(id) > 0){
  CH.p <- CH.p[-id, ]
  CH.sur <- CH.sur[-id,]
}

  return(list(CH.p = CH.p, 
              CH.sur = CH.sur, 
              B = B,
              Nt = Nt))
}

# Execute simulation function

sim <- simul.js(PSI.state = PSI.state, PSI.obs = PSI.obs, 
                gamma_U = gamma_U, gamma_I = gamma_I, N = N,
                N_U = N_U, N_I = N_I)

CH <- sim$CH.p

sim$CH.sur
#------------



