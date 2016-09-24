
# Analysis of the JS model as a multistate model
CH <- sim$CH.p
remove1 <- numeric()

for(i in 1:dim(CH)[1]){

  remove1[i] <- length(which(CH[i,] == 0 | CH[i,] == 3)) == 4

}

if(sum(remove1) > 0) {
CH <- CH[-which(remove1 == 1),]
}

# Add dummy occasion
CH.du <- cbind(rep(0, dim(CH)[1]), CH)

# Augment data
nz <- 500
CH.ms <- rbind(CH.du, matrix(0, ncol = dim(CH.du)[2], nrow = nz))

# Recode CH matrix: a 0 is not allowed in WinBUGS!
CH.ms[CH.ms==0] <- 3        # Not seen = 3, seen = 1 or 2

# Bundle data
jags.data <- list(y = CH.ms, 
                  n.occasions = dim(CH.ms)[2], 
                  M = dim(CH.ms)[1],
                  state = 3)

# Initial values
#--------------------- Initial values

ch <- CH.du

js.multistate.init <- function(ch, nz){
  
  ch[ch==0] <- NA
  ch[ch==3] <- NA
  
  state <- ch
  colnames(state) <- 1:ncol(state)
  
  # When the individual is known to be alive between first and last capture fill it in with 2s
  
  for (i in 1:nrow(ch)){
    
    n1 <- min(which(ch[i,] < 3))
    n2 <- max(which(ch[i,] < 3))
    
    fill <- which(is.na(state[i,n1:n2]) == TRUE)
    
    state[i, names(fill)] <- 2
      #sample(c(1, 2), size = length(fill), replace = TRUE)
  }
  
  state <- state + 1
  
  get.first <- function(x) min(which(!is.na(x)))
  get.last <- function(x) max(which(!is.na(x)))   
  f <- apply(state, 1, get.first)
  l <- apply(state, 1, get.last)
  
  for (i in 1:nrow(ch)){
    
# Before initial observation- not observed
    state[i,1:(f[i]-1)] <- 1
    
    if(l[i]!= ncol(ch)){state[i, (l[i]+1):ncol(ch)] <- 4}
    
      if(ch[i, f[i]] == 1){state[i, f[i]] <- 2}
      if(ch[i, f[i]] == 2){state[i, f[i]] <- 3} 
    
  }   
  
  state <- rbind(state, matrix(1, ncol = ncol(ch), nrow = nz))
  state <- cbind(rep(NA, times = nrow(state)), state[,-1])
  return(state)
}

n.occasions <- dim(CH.ms)[2]

zinit <- js.multistate.init(CH.du, nz)

inits <- function(){list(phi_U = runif(1, 0.5, 1), 
                         phi_I = runif(1, 0.5, 1), 
                         
                         p_U = runif(1, 0.5, 1), 
                         p_I = runif(1, 0.5, 1), 
                         
                         beta_UI = runif(1, 0.5, 1),
                         beta_IU = runif(1, 0.5, 1),
                         
                         gamma_U = runif((n.occasions-1), 0, 0.5), 
                         gamma_I = runif((n.occasions-1), 0, 0.5),
                         
                         z = zinit
                         )}    

# Parameters monitored
params <- c("phi_U", 
            "phi_I", 
            "p_U", 
            "p_I",
            "beta_UI",
            "beta_IU",
            "gamma_U",
            "gamma_I",
             "N_U", "N_I", "B_U", "B_I",
            "zzz.fit", "zzz.fit.new")

# MCMC settings
ni <- 10000
nb <- 1000
nt <- 10
nc <- 3

#-- call Library

library("jagsUI")

# Call JAGS from R
setwd("/Users/Cici/Desktop/")

js.ms <- jags(data = jags.data, inits = inits, parameters.to.save = params, 
              model.file = "model_JS.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

#library("R2jags")
#
#system.time(out <- jags(jags.data, inits, params, "model_JS.txt", n.chains = nc,
#                        n.thin = nt, n.iter = ni, n.burnin = nb))



