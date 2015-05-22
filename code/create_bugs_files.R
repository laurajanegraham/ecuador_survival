############ WinBUGS code to estimate species survival #########################

# code for survival with constant parameters -----------------------------------
sink("~/.wine/drive_c/temp/Rtmp/cjs-c-c.bug")
cat("
model{
        for(i in 1:nind){
                for(t in f[i] : (n.occasions - 1)) {
                        phi[i, t] <- mean.phi
                        p[i, t] <- mean.p
                }
        }
        mean.phi ~ dunif(0, 1)
        mean.p ~ dunif(0, 1)

        for(i in 1:nind) {
                z[i, f[i]] <- 1
                for(t in (f[i] + 1):n.occasions) {
                        z[i, t] ~ dbern(mu1[i, t])
                        mu1[i, t] <- phi[i, t - 1] * z[i, t - 1]
                        y[i, t] ~ dbern(mu2[i, t])
                        mu2[i, t] <- p[i, t - 1] * z[i, t]
                }
        }
}
    ", fill = TRUE)
sink()

# code for survival with fixed time effects ------------------------------------
sink("~/.wine/drive_c/temp/Rtmp/cjs-fixedtime.bug")
cat("
model{
        for(i in 1:nind){
                for(t in f[i] : (n.occasions - 1)) {
                        phi[i, t] <- alpha[t]
                        p[i, t] <- beta[t]
                }
        }
        for (t in 1:n.occasions - 1) {
                alpha[t] ~ dunif(0, 1)
                beta[t] ~ dunif(0, 1)
        }

        for(i in 1:nind) {
                z[i, f[i]] <- 1
                for(t in (f[i] + 1):n.occasions) {
                        z[i, t] ~ dbern(mu1[i, t])
                        mu1[i, t] <- phi[i, t - 1] * z[i, t - 1]
                        y[i, t] ~ dbern(mu2[i, t])
                        mu2[i, t] <- p[i, t - 1] * z[i, t]
                }
        }
}
    ", fill = TRUE)
sink()

# code for survival with fixed group effects -----------------------------------
sink("~/.wine/drive_c/temp/Rtmp/cjs-group.bug")
cat("
    model {
    
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    phi[i,t] <- phi.g[group[i]]
    p[i, t] <- mean.p
    } #t
    } #i
    for (u in 1:g){
    phi.g[u] ~ dunif(0, 1)              # Priors for group-specific survival
    }
    mean.p ~ dunif(0, 1)
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- phi[i,t-1] * z[i,t-1]
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# code for survival with random time effects (can then build in covariates) ----
sink("~/.wine/drive_c/temp/Rtmp/cjs-time-raneff.bug")
cat("
model {

    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu + epsilon[t]
    p[i,t] <- mean.p
    } #t
    } #i
    for (t in 1:(n.occasions-1)){
    epsilon[t] ~ dnorm(0, tau)
    }
    
    #mu ~ dnorm(0, 0.001)                    # Prior for logit of mean survival
    #mean.phi <- 1 / (1+exp(-mu))            # Logit transformation
    mean.phi ~ dunif(0, 1)                   # Prior for mean survival
    mu <- log(mean.phi / (1-mean.phi))       # Logit transformation
    sigma ~ dunif(0, 10)                     # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)                  # Temporal variance
    mean.p ~ dunif(0, 1)                     # Prior for mean recapture
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- phi[i,t-1] * z[i,t-1]
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# code for multinomial formulation ---------------------------------------------
sink("~/.wine/drive_c/temp/Rtmp/cjs-mnl.bug")
cat("
model {
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
    phi[t] ~ dunif(0, 1)         # Priors for survival
    p[t] ~ dunif(0, 1)           # Priors for recapture
    }
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
    marr[t,1:n.occasions] ~ dmulti(pr[t, ], r[t])
    }
    # Calculate the number of birds released each year
    for (t in 1:(n.occasions-1)){
    r[t] <- sum(marr[t, ])
    }
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]                # Probability of non-recapture
    pr[t,t] <- phi[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
    } #t
    
    # Assess model fit using Freeman-Tukey statistic
    # Compute fit statistics for observed data
    for (t in 1:(n.occasions-1)){
    for (j in 1:n.occasions){
    expmarr[t,j] <- r[t]*pr[t,j]
    E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
    } #j
    } #t
    # Generate replicate data and compute fit stats from them
    for (t in 1:(n.occasions-1)){
    marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], r[t])
    for (j in 1:n.occasions){
    E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
    } #j
    } #t
    fit <- sum(E.org[,])
    fit.new <- sum(E.new[,])
    }
    ",fill = TRUE)
sink()

# code for multinomial formulation with random time effects --------------------
sink("~/.wine/drive_c/temp/Rtmp/cjs-mnl-ran-time.bug")
cat("
model {
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
        logit(phi[t]) <- mu + epsilon[t]
        epsilon[t] ~ dnorm(0, tau)
        p[t] <- mean.p
        phi.est[t] <- 1 / (1 + exp())
    }
    mean.phi ~ dunif(0, 1)             # Prior for mean survival
    mu <- log(mean.phi / (1-mean.phi)) # Logit transformation
    sigma ~ dunif(0, 5)                # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)
    # Temporal variance on real scale
    sigma2.real <- sigma2 * pow(mean.phi, 2) * pow((1-mean.phi), 2) 
    mean.p ~ dunif(0, 1)                   # Prior for mean recapture
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
        marr[t,1:n.occasions] ~ dmulti(pr[t,], r[t])
    }
    # Calculate the number of birds released each year
    for (t in 1:(n.occasions-1)){
        r[t] <- sum(marr[t,])
    }
    # Define the cell probabilities of the m-array:
    # Main diagonal
    for (t in 1:(n.occasions-1)){
        q[t] <- 1-p[t]
        pr[t,t] <- phi[t]*p[t]	
        # Above main diagonal
        for (j in (t+1):(n.occasions-1)){
            pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
        } #j	
    # Below main diagonal
        for (j in 1:(t-1)){
            pr[t,j]<-0
        } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
        pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
    } # t

    # Assess model fit using Freeman-Tukey statistic
    # Compute fit statistics for observed data
    for (t in 1:(n.occasions-1)){
        for (j in 1:n.occasions){
            expmarr[t,j] <- r[t]*pr[t,j]
            E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
        }
    }
    # Generate replicate data and compute fit stats from them
    for (t in 1:(n.occasions-1)){
        marr.new[t,1:n.occasions] ~ dmulti(pr[t,], r[t])
        for (j in 1:n.occasions){
            E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
        }
    }
    fit <- sum(E.org[,])
    fit.new <- sum(E.new[,])
}
    ",fill = TRUE)
sink()

# code for multinomial formulation with fixed group effects --------------------
sink("~/.wine/drive_c/temp/Rtmp/cjs-mnl-habitat.bug")
cat("
    model {
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
    phi.intro[t] <- mean.phiintro
    phi.native[t] <- mean.phinative
    phi.scrub[t] <- mean.phiscrub
    p[t] <- mean.p
    }
    mean.phiintro ~ dunif(0, 1)          # Prior for mean survival in introduced habitat
    mean.phinative ~ dunif(0, 1)           # Prior for mean survival in native habitat
    mean.phiscrub ~ dunif(0, 1)           # Prior for mean survival in scrub
    mean.p ~ dunif(0, 1)               # Prior for mean recapture
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
    marr.i[t,1:n.occasions] ~ dmulti(pr.i[t,], r.i[t])
    marr.n[t,1:n.occasions] ~ dmulti(pr.n[t,], r.n[t])
    marr.s[t,1:n.occasions] ~ dmulti(pr.s[t,], r.s[t])
    }
    # Calculate the number of birds released each year
    for (t in 1:(n.occasions-1)){
    r.i[t] <- sum(marr.i[t,])
    r.n[t] <- sum(marr.n[t,])
    r.s[t] <- sum(marr.s[t,])
    }
    # Define the cell probabilities of the m-arrays
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]            # Probability of non-recapture
    pr.i[t,t] <- phi.intro[t]*p[t]
    pr.n[t,t] <- phi.native[t]*p[t]
    pr.s[t,t] <- phi.scrub[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr.i[t,j] <- prod(phi.intro[t:j])*prod(q[t:(j-1)])*p[j]
    pr.n[t,j] <- prod(phi.native[t:j])*prod(q[t:(j-1)])*p[j]
    pr.s[t,j] <- prod(phi.scrub[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.i[t,j] <- 0
    pr.n[t,j] <- 0
    pr.s[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr.i[t,n.occasions] <- 1-sum(pr.i[t,1:(n.occasions-1)])
    pr.n[t,n.occasions] <- 1-sum(pr.n[t,1:(n.occasions-1)])
    pr.s[t,n.occasions] <- 1-sum(pr.s[t,1:(n.occasions-1)])
    } #t

    # Assess model fit using Freeman-Tukey statistic
    # Compute fit statistics for observed data
    for (t in 1:(n.occasions-1)){
        for (j in 1:n.occasions){
            expmarr.i[t,j] <- r.i[t]*pr.i[t,j]
            expmarr.n[t,j] <- r.n[t]*pr.n[t,j]
            expmarr.s[t,j] <- r.s[t]*pr.s[t,j]
            E.org.i[t,j] <- pow((pow(marr.i[t,j], 0.5)-pow(expmarr.i[t,j], 0.5)), 2)
            E.org.n[t,j] <- pow((pow(marr.n[t,j], 0.5)-pow(expmarr.n[t,j], 0.5)), 2)
            E.org.s[t,j] <- pow((pow(marr.s[t,j], 0.5)-pow(expmarr.s[t,j], 0.5)), 2)
        }
    }
    # Generate replicate data and compute fit stats from them
    for (t in 1:(n.occasions-1)){
        marr.new.i[t,1:n.occasions] ~ dmulti(pr.i[t,], r.i[t])
        marr.new.n[t,1:n.occasions] ~ dmulti(pr.n[t,], r.n[t])
        marr.new.s[t,1:n.occasions] ~ dmulti(pr.s[t,], r.s[t])
        for (j in 1:n.occasions){
            E.new.i[t,j] <- pow((pow(marr.new.i[t,j], 0.5)-pow(expmarr.i[t,j], 0.5)), 2)
            E.new.n[t,j] <- pow((pow(marr.new.n[t,j], 0.5)-pow(expmarr.n[t,j], 0.5)), 2)
            E.new.s[t,j] <- pow((pow(marr.new.s[t,j], 0.5)-pow(expmarr.s[t,j], 0.5)), 2)
        }
    }
    fit <- sum(E.org.i[,]) + sum(E.org.n[,]) + sum(E.org.s[,])
    fit.new <- sum(E.new.i[,]) + sum(E.new.n[,]) + sum(E.new.s[,])
    }
    ",fill = TRUE)
sink()
