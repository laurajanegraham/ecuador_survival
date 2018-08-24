model {
  # here i'm going to work out if detectability varies across habitats
  for(s in 1:nspec){ # Define priors for phi and p parameters
    for(a in 1:2) { # loop over marking occasions - we are estimating first banding survival seperately
      for(h in 1:3) { # loop over habitats - (1) introduced, (2) native, (3) scrub
        logit(phi[s,a,h]) <- lphi[s,a,h] + beta.group*group[s] + beta.phi.fam[fam[s]] + epsilon
        lphi[s,a,h] ~ dnorm(mu.lphi[a,h], tau.lphi[a,h])   # Prior for logit(survival): 
      }
    }
    # here we specify a stat. distribution for the logit transformation
    # of the survival probability of species s, AND this statistical
    # distribution has two hyperparameters which we estimate
    for(h in 1:3){ # get a specific detectability per habitat
      p[s,h] <- exp(lp[s,h]) / (1 + exp(lp[s,h]))
      lp[s,h] ~ dnorm(mu.lp[h], tau.lp[h])
    }
  }
  
  # The next lines give (hyper-)priors for the hyperparameters that
  #  characterise the entire community of species
  for(a in 1:2) {
    for(h in 1:3) { 
      mu.lphi[a,h] <- logit(mean.phi[a,h]) # Hyperpriors for survival hyperparams
      mean.phi[a,h] ~ dunif(0,1)       # mean hyperparam. (community avge. survival)
      tau.lphi[a,h] <- pow(sd.lphi[a,h], -2)
      sd.lphi[a,h] ~ dunif(0, 3)       # sd hyperparam. (community variability surv.)
    }
  }
  
  for(h in 1:3){
    mu.lp[h] <- logit(mean.p[h])      # Hyperpriors for recapture hyperparams
    mean.p[h] ~ dunif(0, 1)         # mean hyperparameter (community average recapt.)
    tau.lp[h] <- pow(sd.lp[h], -2)
    sd.lp[h] ~ dunif(0, 3)         # sd hyperparam. (community heterogeneity recap.)
    
  }
  
  # Family effects in the mean
  for (j in 1:nfam){ 
    beta.phi.fam[j] ~ dnorm(0, tau.fam)
  }
  tau.fam <- pow(sd.fam, -2)
  sd.fam ~ dunif(0, 10)
  
  # prior on the group effect
  beta.group ~ dunif(-1,1)
  
  # prior on the regression error
  epsilon ~ dnorm(0, tau.epsilon)
  sd.epsilon ~ dunif(0, 1)          
  tau.epsilon <- pow(sd.epsilon, -2)
  
  # Likelihood for both subsets of data at once (i.e., for the single data set)
  for(s in 1:nspec){           # Loop over species
    for(i in 1:nind[s]){       # Loop over individuals (nind varies by species)
      # Define latent state at first capture
      z[i,f[i,s], s] <- 1
      for(t in (f[i,s]+1):n.occ){ # Loop over occasions
        # State process: the latent alive/dead state
        z[i,t,s] ~ dbern(mu1[i,t,s])
        mu1[i,t,s] <- z[i,t-1,s] * phi[s,tsm[i,t-1,s],hab[i,t-1,s]]
        # Observation process: relates true state to observed state, y = ch
        y[i,t,s] ~ dbern(mu2[i,t,s])
        mu2[i,t,s] <- z[i,t,s] * p[s,hab[i,t,s]]     # p also indexed by species
      } #t
    } #i
  } #s
  
  # mean and sd survival for each species
  for(s in 1:nspec) {
    speciesp[s] <- mean(p[s,])
    speciesphi1[s] <- mean(phi[s,1,])  
    speciesphi1sd[s] <- sd(phi[s,1,])
    speciesphi2[s] <- mean(phi[s,2,])  
    speciesphi2sd[s] <- sd(phi[s,2,])
  }
  
  # differences in survival between habitats
  indiff <- mean.phi[2,1] - mean.phi[2,2]
  isdiff <- mean.phi[2,1] - mean.phi[2,3]
  nsdiff <- mean.phi[2,2] - mean.phi[2,3]
  
  # community mean survival
  allmean.phi <- mean(mean.phi[2,])
  allsd.phi <- mean(sd.lphi[2,])
  
  # # mean survival for specialisms (could not work out how to do in a loop due to diff # species)
  # for(s in 1:gsampsize[1]){
  #   temp1p[s] <- speciesp[group[1,s]]
  #   temp1phi[s] <- speciesphi[group[1,s]]
  # }
  # groupp[1] <- mean(temp1p)
  # groupphi[1] <- mean(temp1phi)
  # 
  # for(s in 1:gsampsize[2]){
  #   temp2p[s] <- speciesp[group[2,s]]
  #   temp2phi[s] <- speciesphi[group[2,s]]
  # }
  # groupp[2] <- mean(temp2p)
  # groupphi[2] <- mean(temp2phi)
  # 
  # for(s in 1:gsampsize[3]){
  #   temp3p[s] <- speciesp[group[3,s]]
  #   temp3phi[s] <- speciesphi[group[3,s]]
  # }
  # groupp[3] <- mean(temp3p)
  # groupphi[3] <- mean(temp3phi)
}
