################################################################################
# WinBUGS code to estimate species survival ------------------------------------
################################################################################
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
    p[i,t] <- p.g[group[i]]
    } #t
    } #i
    for (u in 1:g){
    phi.g[u] ~ dunif(0, 1)              # Priors for group-specific survival
    p.g[u] ~ dunif(0, 1)                # Priors for group-specific recapture
    }
    
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