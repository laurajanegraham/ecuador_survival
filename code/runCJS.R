# runCJS.R----------------------------------------------------------------------
# Laura Graham December 2015

# The purpose of this script is to run a set of different models in the Cormack-Jolly-Seber 
# formation to estimate bird survival. The models are as follows: 

# 1. Null model (survival rate constant over time and individuals)
# 2. Time since marked (takes into account species which never return - transients)
# 3. Habitat model (survival rate varies across habitat type)
# 4. Time model (survival rate varies across time - time is a random factor so we get estimate of variance)
#       a) no covariate
#       b) evi as a covariate
#       c) ndvi as a covariate
#       d) temperature as a covariate

#load packages
require(R2jags)
require(R2WinBUGS)
require(mcmcplots)
require(dplyr)
require(Hmisc)
require(tidyr)
require(lubridate)
library(abind)
library(snowfall)

# load functions
source("code/fnCleanBandingDat.R")
source("code/fnCleanRSData.R")
source("code/fnEncounterHistory.R")
source("code/fnKnownStateInits.R")
source("code/JAGSParallel.R")

# Load and format banding and remote sensing data ------------------------------
banding.dat.clean <- CleanBandingDat()
sessions <- unique(banding.dat.clean$session_new)
#RS.dat <- CleanRSData(12)

# this is the log file that will record how long each model run took. 
logfile.name <- paste0("logs/winbugs_log_", Sys.Date(), ".txt")
file.create(logfile.name)

# get a count of how many individual captures exist for each species and then 
# use those species with more than 100 for the initial analyses (there are 10 of
# these species)
species.list <- group_by(banding.dat.clean, Specie.Name) %>%
    summarise(count = length(Specie.Name)) %>%
    filter(count > 50)

# this will create a list of species where the habitat analysis wasn't run
# (because it's not present in all 3 habitats)

# same MCMC options for each model - may need to think about changing this
ni <- 10000
nt <- 3
nb <- 2000
nc <- 3

for (species in species.list$Specie.Name){
    # if there are already results for the species (i.e. for a different model), load so they are all together
    if(file.exists(paste0("results/", species, "_CJS_model_output.rda"))) {
        load(paste0("results/", species, "_CJS_model_output.rda"))
    } else {
        # otherwise create a list to put results into
        modelout <- list()
    }
    # get species data ---------------------------------------------------------
    sp_dat <- filter(banding.dat.clean, Specie.Name == species) %>%
        select(Band.Number, session_new, habitat)
    
    # create encounter histories -----------------------------------------------
    sp_eh <- EncounterHistory(sp_dat, "session_new", "Band.Number", "habitat", sessions=sessions)
    
    # get the different marrays required by each of the models -----------------
    marr <- sp_eh$m.array
    marr.gp <- sp_eh$m.array.gp
    #marr.TSM1 <- sp_eh$m.array.TSM1
    #marr.TSM2 <- sp_eh$m.array.TSM2
    
    # run the models and output to logfile
    # Constant/null ----------------------------------------------------------------------------------------
    write(paste("Null model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    constant.data <- list(marr = marr, n.occasions = ncol(marr), r = rowSums(marr))
    constant.inits <- function(){list(mean.phi = runif(1, 0, 1), 
                                      mean.p = runif(1, 0, 1))}  
    constant.parameters <- c("mean.p", "mean.phi", "fit", "fit.new")
    modelout[["null"]] <- JAGSParallel(nc, data=constant.data, inits=constant.inits, params=constant.parameters, 
                                       model.file="cjs-mnl-fixed.bug", n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
#     # TSM --------------------------------------------------------------------------------------------------
#     write(paste("TSM model for ", species, sep=" "), logfile.name, append = TRUE)
#     write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
#     strt <- Sys.time()       
#     tsm.data <- list(marr.TSM1 = marr.TSM1, marr.TSM2 = marr.TSM2, n.occasions = ncol(marr.TSM1))
#     tsm.inits <- function(){list(mean.TSM1 = runif(1, 0, 1), 
#                                  mean.TSM2 = runif(1, 0, 1), 
#                                  mean.p = runif(1, 0, 1))}  
#     tsm.parameters <- c("mean.p", "mean.phitsm1", "mean.phitsm2","fit", "fit.new")
#     modelout[["TSM"]] <- JAGSParallel(nc, data=tsm.data, inits=tsm.inits, params=tsm.parameters, 
#                                       model.file="cjs-mnl-TSM.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
#     write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
#     
    # Habitat --------------------------------------------------------------------------------------------------
    # Habitat is more complicated because some species are only in two habitats
    write(paste("Habitat model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    if(length(marr.gp)==3) {
        habitat.data <- list(marr.i = marr.gp[["Introduced"]], marr.n = marr.gp[["Native"]],
                             marr.s = marr.gp[["Scrub"]], n.occasions = ncol(marr),
                             r.i = rowSums(marr.gp[["Introduced"]], r.n = rowSums(marr.gp[["Native"]]), 
                                           r.s = rowSums(marr.gp[["Scrub"]])))
        
        habitat.inits <- function(){list(mean.phiintro = runif(1, 0, 1), 
                                         mean.phinative = runif(1, 0, 1),
                                         mean.phiscrub = runif(1, 0, 1),
                                         mean.p = runif(1, 0, 1))} 
        habitat.parameters <- c("mean.p", "mean.phinative", "mean.phiscrub", "mean.phiintro", "fit", "fit.new")
        modelout[["habitat"]] <- JAGSParallel(nc, data=habitat.data, inits=habitat.inits, params=habitat.parameters, 
                                              model.file="cjs-mnl-habitat.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
        write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), 
                    sep = " "), logfile.name, append = TRUE)
    }
    if(length(marr.gp)==2) {
        # these cases will need a little more manipulation at the analysis stage to label habitat type
        habitat.data <- list(marr.1 = marr.gp[[1]], marr.2 = marr.gp[[2]], n.occasions = ncol(marr),
                             r.1 = rowSums(marr.gp[[1]]), r.2 = rowSums(marr.gp[[2]]))
        
        habitat.inits <- function(){list(mean.phi1 = runif(1, 0, 1), 
                                         mean.phi2 = runif(1, 0, 1),
                                         mean.p = runif(1, 0, 1))} 
        habitat.parameters <- c("mean.p", "mean.phi1", "mean.phi2", "fit", "fit.new")
        modelout[["habitat"]] <- JAGSParallel(nc, data=habitat.data, inits=habitat.inits, params=habitat.parameters, 
                                              model.file="cjs-mnl-habitat-2groups.bug", n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb)
        write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), 
                    sep = " "), logfile.name, append = TRUE)
    }
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # Time -----------------------------------------------------------------------------------------------------
    write(paste("Time model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    time.data <- list(marr = marr, n.occasions = ncol(marr), r=rowSums(marr))
    time.inits <- function(){list(mean.phi = runif(1, 0, 1), 
                                  sigma = runif(1, 0, 5),
                                  mean.p = runif(1, 0, 1))}  
    time.parameters <- c("phi", "mean.p", "mean.phi", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["time"]] <- JAGSParallel(nc, data=time.data, inits=time.inits, params=time.parameters, 
                                       model.file="cjs-mnl-ran-time.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
#     # EVI ---------------------------------------------------------------------------------------------------
#     write(paste("EVI model for ", species, sep=" "), logfile.name, append = TRUE)
#     write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
#     strt <- Sys.time()       
#     evi.data <- list(marr = marr, n.occasions = ncol(marr), x = RS.dat$evi)
#     timecov.inits <- function(){list(mean.phi = runif(1, 0, 1), 
#                                      sigma = runif(1, 0, 5), 
#                                      mean.p = runif(1, 0, 1), 
#                                      beta = runif(1, -5, 5))}  
#     timecov.parameters <- c("phi", "mean.p", "mean.phi", "beta", "sigma2", "sigma2.real", "fit", "fit.new")
#     modelout[["evi"]] <- JAGSParallel(nc, data=evi.data, inits=timecov.inits, params=timecov.parameters, 
#                                       model.file="cjs-mnl-ran-time-cov.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
#     write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
#     
#     # NDVI ---------------------------------------------------------------------------------------------------
#     write(paste("NDVI model for ", species, sep=" "), logfile.name, append = TRUE)
#     write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
#     strt <- Sys.time()       
#     ndvi.data <- list(marr = marr, n.occasions = ncol(marr), x = RS.dat$ndvi)
#     modelout[["ndvi"]] <- JAGSParallel(nc, data=ndvi.data, inits=timecov.inits, params=timecov.parameters, 
#                                        model.file="cjs-mnl-ran-time-cov.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
#     write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
#     
#     # temp --------------------------------------------------------------------------------------------------
#     write(paste("Temperature model for ", species, sep=" "), logfile.name, append = TRUE)
#     write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
#     strt <- Sys.time()       
#     temp.data <- list(marr = marr, n.occasions = ncol(marr), x = RS.dat$temp)
#     modelout[["temp"]] <- JAGSParallel(nc, data=temp.data, inits=timecov.inits, params=timecov.parameters, 
#                                        model.file="cjs-mnl-ran-time-cov.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
#     write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
     
    save(modelout, file=paste0("results/", species, "_CJS_model_output.rda"))
}