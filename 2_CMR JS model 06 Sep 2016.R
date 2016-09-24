# The JS model as a multistate model 
# Specify model in BUGS language

setwd("/Users/Cici/Desktop/")
sink("model_JS.txt")
cat("

model {
    
#--------------------------------------
# Parameters:
  # phi_U: survival probability of uninfected
  # phi_I: survival probability of infected
  # gamma_U: removal entry probability of uninfected
  # gamma_I: removal entry probability of infected
  # p_U: capture probability of uninfected
  # p_I: capture probability of infected
#--------------------------------------
# States (S):
  # 1 not yet entered
  # 2 uninfected
  # 3 infected
  # 4 dead
# Observations (O):
  # 1 seen as uninfected
  # 2 seen as infected
  # 3 not seen
#--------------------------------------
    
# Priors and constraints

phi_U ~ dunif(0, 1)   # Prior for mean survival
phi_I ~ dunif(0, 1)   # Prior for mean survival

for(t in 1:(n.occasions -1)){
  gamma_U[t] ~ dunif(0, 1)   # Prior for entry probabilities
  gamma_I[t] ~ dunif(0, 1)   # Prior for entry probabilities
}

p_U ~ dunif(0, 1)     # Prior for mean capture
p_I ~ dunif(0, 1)     # Prior for mean capture

beta_UI ~ dunif(0, 1)
beta_IU ~ dunif(0, 1)

# State-transition and observation matrices 	

for (i in 1:M){  

  # Define probabilities of state S(t+1) given S(t)

    for (t in 1:(n.occasions-1)){

        ps[1,i,t,1] <- 1 - gamma_U[t] - gamma_I[t]
        ps[1,i,t,2] <- gamma_U[t]
        ps[1,i,t,3] <- gamma_I[t]
        ps[1,i,t,4] <- 0

        ps[2,i,t,1] <- 0
        ps[2,i,t,2] <- phi_U * (1 - beta_UI)
        ps[2,i,t,3] <- phi_U * beta_UI
        ps[2,i,t,4] <- 1 - phi_U

        ps[3,i,t,1] <- 0
        ps[3,i,t,2] <- phi_I * beta_IU
        ps[3,i,t,3] <- phi_I * (1 - beta_IU)
        ps[3,i,t,4] <- 1 - phi_I

        ps[4,i,t,1] <- 0
        ps[4,i,t,2] <- 0
        ps[4,i,t,3] <- 0
        ps[4,i,t,4] <- 1

# Define probabilities of O(t) given S(t)
        po[1,i,t,1] <- 0
        po[1,i,t,2] <- 0
        po[1,i,t,3] <- 1

        po[2,i,t,1] <- p_U
        po[2,i,t,2] <- 0
        po[2,i,t,3] <- 1 - p_U

        po[3,i,t,1] <- 0
        po[3,i,t,2] <- p_I
        po[3,i,t,3] <- 1 - p_I

        po[4,i,t,1] <- 0
        po[4,i,t,2] <- 0
        po[4,i,t,3] <- 1

    } #t

} #i

# Likelihood 
for (i in 1:M){

  # Define latent state at first occasion

    z[i, 1] <- 1   # Make sure that all M individuals are in state 1 at t = 1
  
  for (t in 2:n.occasions){

    # State process: draw S(t) given S(t-1)

      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, ])

    # Observation process: draw O(t) given S(t)

      y[i,t] ~ dcat(po[z[i,t], i, t-1, ])
      y.new[i, t] ~ dcat(po[z[i, t], i, t-1, ])

  } #t

} #i

# Calculate derived population parameters
for (t in 1:(n.occasions-1)){
  qgamma_U[t] <- 1-gamma_U[t]
  qgamma_I[t] <- 1-gamma_I[t]
}

cprob_U[1] <- gamma_U[1]
cprob_I[1] <- gamma_I[1]

for (t in 2:(n.occasions-1)){
  cprob_U[t] <- gamma_U[t] * prod(qgamma_U[1:(t-1)])
  cprob_I[t] <- gamma_I[t] * prod(qgamma_I[1:(t-1)])
} #t

psi_U <- sum(cprob_U[])            # Inclusion probability
psi_I <- sum(cprob_I[])            # Inclusion probability

for (t in 1:(n.occasions-1)){
  b_U[t] <- cprob_U[t] / psi_U      # Entry probability
  b_I[t] <- cprob_I[t] / psi_I      # Entry probability
} #t

for (i in 1:M){
  for (t in 2:n.occasions){
    al_U[i,t-1] <- equals(z[i,t], 2)
    al_I[i,t-1] <- equals(z[i,t], 3)
  } #t
  for (t in 1:(n.occasions-1)){
    d_U[i,t] <- equals(z[i,t] - al_U[i,t], 0)
    d_I[i,t] <- equals(z[i,t] - al_I[i,t], 0)
  } #t   
  alive[i] <- sum(al_U[i,], al_I[i,])
} #i

for (t in 1:(n.occasions-1)){
  N_U[t] <- sum(al_U[,t])        # Actual population size
  N_I[t] <- sum(al_I[,t])        # Actual population size
  B_U[t] <- sum(d_U[,t])         # Number of entries
  B_I[t] <- sum(d_I[,t])         # Number of entries
} #t

for (i in 1:M){
  w[i] <- 1-equals(alive[i],0)
} #i

Nsuper <- sum(w[])            # Superpopulation size

#---------------- Calculate Bayesian posterior predictive check


#---- r is puts a 1 when an individual
#---- was seen in a specific state

for(t in 1:(n.occasions-1)){

  for(i in 1:M){

    for(s in 1:state){

      r[s, i, t] <- ifelse(y[i, t+1] == s, 1, 0)

      r.new[s, i, t] <- ifelse(y.new[i, t+1] == s, 1, 0)

    }
  }
}

#---- R is the total number of times an individual 
#---- was seen in a specific state

for(t in 1:(n.occasions-1)){

  for(s in 1:state){

    # sum across individuals in each state, each time period

      R_state[s, t] <- sum(r[s, , t])

      R_state.new[s, t] <- sum(r.new[s, , t])

  }

}

#---- PO1 adds up the probabilities that an 
#---- an individual was observed in that state = expected

for(t in 1:(n.occasions-1)){

  for(s in 1:state){

    for(i in 1:M){

        muy[i, t, s] <-  ps[z[i, t], i, t, z[i, t+1]] * po[z[i, t], i, t, s]

    } 

    PO_expt[s, t] <- sum(0.01, sum(muy[ , t, s]))

  }

}

#--------- Posterior predictive check

for(t in 1:(n.occasions-1)){

  for(s in 1:state){

      E.act[s, t] <- pow(pow(R_state[s, t], 0.5) - pow(PO_expt[s, t], 0.5), 2)

      E.new[s, t] <- pow(pow(R_state.new[s, t], 0.5) - pow(PO_expt[s, t], 0.5), 2)


  }

}


zzz.fit <- sum(E.act[,])
zzz.fit.new <- sum(E.new[,])

}
",fill = TRUE)
sink()
