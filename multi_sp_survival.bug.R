model {
  for(s in 1:nspec){ # Define priors for phi and p parameters
    phi[s] <- exp(lphi[s]) / (1 + exp(lphi[s]))
    lphi[s] ~ dnorm(mu.lphi, tau.lphi)   # Prior for logit(survival): 
    # here we specify a stat. distribution for the logit transformation
    # of the survival probability of species s, AND this statistical
    # distribution has two hyperparameters which we estimate
    
    p[s] <- exp(lp[s]) / (1 + exp(lp[s]))
    lp[s] ~ dnorm(mu.lp, tau.lp)
  }
  
  # The next lines give (hyper-)priors for the hyperparameters that
  #  characterise the entire community of species
  
  mu.lphi <- logit(mean.phi) # Hyperpriors for survival hyperparams
  mean.phi ~ dunif(0,1)       # mean hyperparam. (community avge. survival)
  tau.lphi <- pow(sd.lphi, -2)
  sd.lphi ~ dunif(0, 3)       # sd hyperparam. (community variability surv.)
  
  mu.lp <- logit(mean.p)      # Hyperpriors for recapture hyperparams
  mean.p ~ dunif(0, 1)         # mean hyperparameter (community average recapt.)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 3)         # sd hyperparam. (community heterogeneity recap.)
  
  # Likelihood for both subsets of data at once (i.e., for the single data set)
  for(s in 1:nspec){           # Loop over species
    for(i in 1:nind[s]){       # Loop over individuals (nind varies by species)
      # Define latent state at first capture
      z[i,f[i,s], s] <- 1
      for(t in (f[i,s]+1):n.occ){ # Loop over occasions
        # State process: the latent alive/dead state
        z[i,t,s] ~ dbern(mu1[i,t,s])
        mu1[i,t,s] <- z[i,t-1,s] * phi[s]
        # Observation process: relates true state to observed state, y = ch
        y[i,t,s] ~ dbern(mu2[i,t,s])
        mu2[i,t,s] <- z[i,t,s] * p[s]     # p also indexed by species
      } #t
    } #i
  } #s
  
}
