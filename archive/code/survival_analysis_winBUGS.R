require(R2WinBUGS)
require(mcmcplots)
source("code/fnCleanBandingDat.R")
source("code/fnEncounterHistory.R")
source("code/fnKnownStateInits.R")

# Load and format banding data -------------------------------------------------
banding.dat.clean <- CleanBandingDat()

# get a count of how many individual captures exist for each species and then 
# use those species with more than 100 for the initial analyses (there are 10 of
# these species)
species.list <- group_by(banding.dat.clean, Specie.Name) %>%
        summarise(count = length(Specie.Name)) %>%
        arrange(-count)

# this is the log file that will record how long each model run took. 
logfile.name <- paste0("logs/winbugs_log_", Sys.time(), ".txt")
file.create(logfile.name)

# initialise lists for storing the results
cjs.group.habitat <- list()
cjs.ran.time <- list()

for (species in species.list$Specie.Name) {
        
        # get records for species, get necessary columns
        sp_dat <- filter(banding.dat.clean, Specie.Name == species) %>%
                select(Band.Number, session_new, habitat)
        
        sp_eh <- EncounterHistory(sp_dat, 
                                  "session_new", 
                                  "Band.Number",
                                  "habitat")
        
        sp_bugs <- sp_eh$eh.full
        
        # Create group variable as a number
        group <- as.numeric(as.factor(sp_bugs$habitat))
        CH <- as.matrix(select(sp_bugs, -Band.Number, -habitat))
        
        # Create vector with occasion of marking
        get.first <- function(x) min(which(x!=0))
        f <- apply(CH, 1, get.first)
        
        # Analysis of survival by habitat type -----------------------------------------
        
        # Bundle data
        bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g = length(unique(group)), group = group)
        
        # Initial values
        inits <- function(){list(z = cjs.init.z(CH, f), phi.g = runif(length(unique(group)), 0, 1), mean.p = runif(1, 0, 1))}  
        
        # Parameters monitored
        parameters <- c("phi.g", "mean.p")
        
        # MCMC settings
        ni <- 10000
        nt <- 3
        nb <- 2000
        nc <- 3
        
        # Call WinBUGS from R
        write(paste("Analysis by habitat for", species, sep=" "), 
                logfile.name, append = TRUE)
        write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), 
              logfile.name, append = TRUE)
        strt<-Sys.time()
        cjs.group.habitat[[species]] <- bugs(bugs.data, inits, parameters, "cjs-group.bug", n.chains = nc, 
                          n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, 
                          working.directory='~/.wine/drive_c/temp/Rtmp/', clearWD=TRUE)
        write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), 
                    sep = " "), logfile.name, append = TRUE)
        
        # Analysis of survival across time (time as random factor) ---------------------
        # bundle data
        bugs.data <- list(y= CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], 
                          z = known.state.cjs(CH))
        
        # initial values
        inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), 
                                 sigma = runif(1, 0, 10), mean.p = runif(1, 0, 1))}
        
        # parameters monitored
        parameters <- c("mean.phi", "mean.p", "sigma2")
        
        # MCMC settings
        ni <- 60000
        nt <- 6
        nb <- 5000
        nc <- 3
        
        write(paste("Analysis with time as random effect for", species, sep=" "), 
              logfile.name, append = TRUE)
        write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), 
              logfile.name, append = TRUE)
        strt <- Sys.time()
        cjs.ran.time[[species]] <- bugs(bugs.data, inits, parameters, "cjs-time-raneff.bug", n.chains = nc, 
                        n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, 
                        working.directory='~/.wine/drive_c/temp/Rtmp/', clearWD=TRUE)
        write(paste("Model run took", round(Sys.time()-strt, 2), "minutes", 
                    sep = " "), logfile.name, append = TRUE)
}

save(cjs.group.habitat, file="results/cjs.group.habitat.rda")
save(cjs.ran.time, file="results/cjs.ran.time.rda")