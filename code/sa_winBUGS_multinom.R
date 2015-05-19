require(R2WinBUGS)
require(mcmcplots)
source("code/fnCleanBandingDat.R")
source("code/fnEncounterHistory.R")
source("code/fnKnownStateInits.R")

# Load and format banding data -------------------------------------------------
banding.dat.clean <- CleanBandingDat()

# this is the log file that will record how long each model run took. 
logfile.name <- paste0("logs/winbugs_log_", Sys.time(), ".txt")
file.create(logfile.name)

# get a count of how many individual captures exist for each species and then 
# use those species with more than 100 for the initial analyses (there are 10 of
# these species)
species.list <- group_by(banding.dat.clean, Specie.Name) %>%
        summarise(count = length(Specie.Name)) %>%
        arrange(-count)

species.list <- species.list[1:5,1]
for (species in species.list$Specie.Name){

        sp_dat <- filter(banding.dat.clean, Specie.Name == species) %>%
                select(Band.Number, session_new, habitat)
        
        sp_eh <- EncounterHistory(sp_dat, "session_new", "Band.Number", "habitat")
        
        marr <- sp_eh$m.array
        marr.gp <- sp_eh$m.array.gp
        
        
        # Code for the random time effects analysis ------------------------------------
        bugs.data <- list(marr = marr, n.occasions = ncol(marr))
        
        inits <- function(){list(mean.phi = runif(1, 0, 1), sigma = runif(1, 0, 5), 
                                 mean.p = runif(1, 0, 1))}  
        
        # Parameters monitored
        parameters <- c("phi", "mean.p", "mean.phi", "sigma2", "sigma2.real", "fit", "fit.new")
        
        ni <- 10000
        nt <- 3
        nb <- 2000
        nc <- 3
        
        write(paste("Analysis by habitat for", species, sep=" "), 
              logfile.name, append = TRUE)
        write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), 
              logfile.name, append = TRUE)
        strt <- Sys.time()       
        cjs <- bugs(bugs.data, inits, parameters, "cjs-mnl-ran-time.bug", n.chains = nc, n.thin = nt, 
                    n.iter = ni, n.burnin = nb, debug = TRUE, 
                    working.directory='~/.wine/drive_c/temp/Rtmp/', clearWD=TRUE)
        
        write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), 
                    sep = " "), logfile.name, append = TRUE)
}

