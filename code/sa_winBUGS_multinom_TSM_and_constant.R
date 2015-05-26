require(R2WinBUGS)
require(mcmcplots)
require(dplyr)
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
    filter(count > 50)

cjs.mnl.tsm <- list()

for (species in species.list$Specie.Name){
    sp_dat <- filter(banding.dat.clean, Specie.Name == species) %>%
        select(Band.Number, session_new, habitat)

    sp_eh <- EncounterHistory(sp_dat, "session_new", "Band.Number", "habitat")
    
    ni <- 10000
    nt <- 3
    nb <- 2000
    nc <- 3
    
    # Code for the TSM analysis ------------------------------------------------
    marr.TSM1 <- sp_eh$m.array.TSM1
    marr.TSM2 <- sp_eh$m.array.TSM2
    
    bugs.data <- list(marr.TSM1 = marr.TSM1, marr.TSM2 = marr.TSM2, n.occasions = ncol(marr.TSM1))
    
    inits <- function(){list(mean.TSM1 = runif(1, 0, 1), 
                             mean.TSM2 = runif(1, 0, 1),
                             mean.p = runif(1, 0, 1))}  
    
    # Parameters monitored
    parameters <- c("mean.p", "mean.phitsm1", "mean.phitsm2",
                    "fit", "fit.new")
    
    
    write(paste("TSM model for ", species, sep=" "), 
          logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), 
          logfile.name, append = TRUE)
    strt <- Sys.time()       
    cjs.mnl.tsm[[species]] <- bugs(bugs.data, inits, parameters, "cjs-mnl-TSM.bug", n.chains = nc, n.thin = nt, 
                        n.iter = ni, n.burnin = nb, debug = FALSE, 
                        working.directory='~/.wine/drive_c/temp/Rtmp/', clearWD=TRUE)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), 
                sep = " "), logfile.name, append = TRUE)
    
    # Code for the constant parameter analysis ------------------------------------------------
    marr = sp_eh$m.array
    
    bugs.data <- list(marr = marr, n.occasions = ncol(marr))
    
    inits <- function(){list(mean.phi = runif(1, 0, 1), 
                             mean.p = runif(1, 0, 1))}  
    
    # Parameters monitored
    parameters <- c("mean.p", "mean.phi", "fit", "fit.new")
    
    
    write(paste("Constant model for", species, sep=" "), 
          logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), 
          logfile.name, append = TRUE)
    strt <- Sys.time()       
    cjs.mnl.tsm[[species]] <- bugs(bugs.data, inits, parameters, "cjs-mnl-fixed.bug", n.chains = nc, n.thin = nt, 
                                   n.iter = ni, n.burnin = nb, debug = FALSE, 
                                   working.directory='~/.wine/drive_c/temp/Rtmp/', clearWD=TRUE)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), 
                sep = " "), logfile.name, append = TRUE)
    
}
