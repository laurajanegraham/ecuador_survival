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
#       c) temperature as a covariate
#       d) evi and temperature as covariates
#       e) evi, temperature and the intersection as covariates

#load packages
library("R2jags")
library("R2WinBUGS")
library("mcmcplots")
library("dplyr")
library("Hmisc")
library("tidyr")
library("lubridate")
library("abind")
library("snowfall")
library("RODBC")

# load functions
source("code/fnCleanBandingDat.R")
source("code/fnCleanRSData.R")
source("code/fnEncounterHistory.R")
source("code/fnKnownStateInits.R")
source("code/JAGSParallel.R")

# change this depending on whether data should be grouped or for individual species
grouped = FALSE

# Load and format banding and remote sensing data ------------------------------
banding.dat.clean <- CleanBandingDat(inc.juvenile=FALSE)
sessions <- unique(banding.dat.clean$session_new)
#bird_groups <- read.csv("data/bird_groups_AN.csv")
#banding.dat.clean <- merge(banding.dat.clean, bird_groups, by.x="Specie.Name", by.y="SP.")
RS.dat <- CleanRSData(12)
# this is the log file that will record how long each model run took. 
logfile.name <- paste0("logs/winbugs_log_", Sys.Date(), ".txt")
file.create(logfile.name)

# Here we create a count of recaptures (how many individuals captured once, twice ... n times)
if(grouped) { 
  banding.dat.clean$tax_level <- banding.dat.clean$Group
} else {
  banding.dat.clean$tax_level <- banding.dat.clean$Specie.Name
}

dat_summary <- group_by(banding.dat.clean, tax_level, Band.Number) %>%
    summarise(recap=n()-1)

dat_recaps <- mutate(dat_summary, 
                     N=1) %>%
    spread(recap, N) %>%
    select(-Band.Number) %>%
    group_by(tax_level) %>%
    summarise_each(funs(sum(., na.rm=TRUE))) %>%
    mutate(nind=rowSums(.[-1]))

dat_nind <- group_by(dat_summary, tax_level) %>%
    summarise(nind=n())

dat_recaps$recaps <- rowSums(dat_recaps[,-c(1, 2, ncol(dat_recaps))])

# we will model the survival for any species for which 2 or more individuals 
# have recaptures. This is likely to be not enough, but we can assess by 
# checking convergence and precision of estimates. The information in this table
# can be put into the paper
usable_dat <- filter(dat_recaps, recaps > 1)
species.list <- usable_dat$tax_level

for (species in species.list){
    # create a list to put results into
    modelout <- list()
    
    # get species data ---------------------------------------------------------
    sp_dat <- filter(banding.dat.clean, tax_level == species) %>%
        select(Band.Number, session_new, habitat)
    
    # create encounter histories -----------------------------------------------
    sp_eh <- EncounterHistory(sp_dat, "session_new", "Band.Number", "habitat", sessions=sessions)
    
    # get the different marrays required by each of the models -----------------
    marr <- sp_eh$m.array
    marr.gp <- sp_eh$m.array.gp
    marr.TSM1 <- sp_eh$m.array.TSM1
    marr.TSM2 <- sp_eh$m.array.TSM2
    marr.TSM1.gp <- sp_eh$m.array.TSM1.gp
    marr.TSM2.gp <- sp_eh$m.array.TSM2.gp
    
    # run the models and output to logfile
    # Constant/null ----------------------------------------------------------------------------------------
    # MCMC options
    ni <- 10000
    nt <- 3
    nb <- 2000
    nc <- 3
    
    write(paste("Null model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    constant.data <- list(marr = marr, n.occasions = ncol(marr), r = rowSums(marr))
    constant.inits <- function(){list(mean.phi2 = runif(1, 0, 1), 
                                      mean.p = runif(1, 0, 1))}  
    constant.parameters <- c("mean.p", "mean.phi2", "fit", "fit.new")
    modelout[["null"]] <- JAGSParallel(nc, data=constant.data, inits=constant.inits, params=constant.parameters, 
                                       model.file="cjs-mnl-fixed.bug", n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
     
    # TSM --------------------------------------------------------------------------------------------------
    write(paste("Null TSM model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()
    tsm.data <- list(marr.TSM1 = marr.TSM1, marr.TSM2 = marr.TSM2, n.occasions = ncol(marr.TSM1),
                     r.TSM1 = rowSums(marr.TSM1), r.TSM2 = rowSums(marr.TSM2))
    tsm.inits <- function(){list(mean.phi1 = runif(1, 0, 1),
                                 mean.phi2 = runif(1, 0, 1),
                                 mean.p = runif(1, 0, 1))}
    tsm.parameters <- c("mean.p", "mean.phi1", "mean.phi2","fit", "fit.new")
    modelout[["TSM"]] <- JAGSParallel(nc, data=tsm.data, inits=tsm.inits, params=tsm.parameters,
                                      model.file="cjs-mnl-TSM.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)

    #Habitat --------------------------------------------------------------------------------------------------
    #Habitat is more complicated because some species are only in two habitats

    write(paste("Habitat model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    if(length(marr.gp)==3) {
        habitat.data <- list(marr.i = marr.gp[["Introduced"]], marr.n = marr.gp[["Native"]],
                             marr.s = marr.gp[["Scrub"]], n.occasions = ncol(marr),
                             r.i = rowSums(marr.gp[["Introduced"]]), r.n = rowSums(marr.gp[["Native"]]), 
                                           r.s = rowSums(marr.gp[["Scrub"]]))
        
        habitat.inits <- function(){list(mean.phi2.intro = runif(1, 0, 1), 
                                         mean.phi2.native = runif(1, 0, 1),
                                         mean.phi2.scrub = runif(1, 0, 1),
                                         mean.p = runif(1, 0, 1))} 
        habitat.parameters <- c("mean.p", "mean.phi2.native", "mean.phi2.scrub", "mean.phi2.intro", "fit", "fit.new")
        modelout[["habitat"]] <- JAGSParallel(nc, data=habitat.data, inits=habitat.inits, params=habitat.parameters, 
                                              model.file="cjs-mnl-habitat.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    }
    if(length(marr.gp)==2) {
        # these cases will need a little more manipulation at the analysis stage to label habitat type
        habitat.data <- list(marr.1 = marr.gp[[1]], marr.2 = marr.gp[[2]], n.occasions = ncol(marr),
                             r.1 = rowSums(marr.gp[[1]]), r.2 = rowSums(marr.gp[[2]]))
        
        habitat.inits <- function(){list(mean.phi2.hab1 = runif(1, 0, 1), 
                                         mean.phi2.hab2 = runif(1, 0, 1),
                                         mean.p = runif(1, 0, 1))} 
        habitat.parameters <- c("mean.p", "mean.phi2.hab1", "mean.phi2.hab2", "fit", "fit.new")
        modelout[["habitat"]] <- JAGSParallel(nc, data=habitat.data, inits=habitat.inits, params=habitat.parameters, 
                                              model.file="cjs-mnl-habitat-2groups.bug", n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb)
        }
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # Habitat TSM model ------------------------------------------------------------------------------------
    write(paste("Habitat TSM model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    if(length(marr.gp)==3) {
        habitat.data <- list(marr.i1 = marr.TSM1.gp[["Introduced"]], 
                             marr.n1 = marr.TSM1.gp[["Native"]],
                             marr.s1 = marr.TSM1.gp[["Scrub"]], 
                             marr.i2 = marr.TSM2.gp[["Introduced"]], 
                             marr.n2 = marr.TSM2.gp[["Native"]],
                             marr.s2 = marr.TSM2.gp[["Scrub"]], 
                             n.occasions = ncol(marr),
                             r.i1 = rowSums(marr.TSM1.gp[["Introduced"]]), 
                             r.n1 = rowSums(marr.TSM1.gp[["Native"]]), 
                             r.s1 = rowSums(marr.TSM1.gp[["Scrub"]]),
                             r.i2 = rowSums(marr.TSM2.gp[["Introduced"]]), 
                             r.n2 = rowSums(marr.TSM2.gp[["Native"]]), 
                             r.s2 = rowSums(marr.TSM2.gp[["Scrub"]]))
        
        habitat.inits <- function(){list(mean.phi1.intro = runif(1, 0, 1), 
                                         mean.phi1.native = runif(1, 0, 1),
                                         mean.phi1.scrub = runif(1, 0, 1),
                                         mean.phi2.intro = runif(1, 0, 1), 
                                         mean.phi2.native = runif(1, 0, 1),
                                         mean.phi2.scrub = runif(1, 0, 1),
                                         mean.p = runif(1, 0, 1))} 
        habitat.parameters <- c("mean.p", "mean.phi1.native", "mean.phi1.scrub", "mean.phi1.intro", "mean.phi2.native", "mean.phi2.scrub", "mean.phi2.intro2", "fit", "fit.new")
        modelout[["habitat.tsm"]] <- JAGSParallel(nc, data=habitat.data, inits=habitat.inits, params=habitat.parameters, 
                                              model.file="cjs-mnl-TSM-habitat.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    }
    if(length(marr.gp)==2) {
        # these cases will need a little more manipulation at the analysis stage to label habitat type
        habitat.data <- list(marr.i1 = marr.TSM1.gp[[1]], 
                             marr.n1 = marr.TSM1.gp[[2]],
                             marr.i2 = marr.TSM2.gp[[1]], 
                             marr.n2 = marr.TSM2.gp[[2]],
                             n.occasions = ncol(marr),
                             r.i1 = rowSums(marr.TSM1.gp[[1]]), 
                             r.n1 = rowSums(marr.TSM1.gp[[2]]), 
                             r.i2 = rowSums(marr.TSM2.gp[[1]]), 
                             r.n2 = rowSums(marr.TSM2.gp[[2]]))
        
        habitat.inits <- function(){list(mean.phi1.hab1 = runif(1, 0, 1), 
                                         mean.phi1.hab2 = runif(1, 0, 1),
                                         mean.phi2.hab1 = runif(1, 0, 1), 
                                         mean.phi2.hab2 = runif(1, 0, 1),
                                         mean.p = runif(1, 0, 1))} 
        habitat.parameters <- c("mean.p", "mean.phi1.hab1", "mean.phi1.hab2", "mean.phi2.hab1", "mean.phi2.hab2", "fit", "fit.new")
        modelout[["habitat.tsm"]] <- JAGSParallel(nc, data=habitat.data, inits=habitat.inits, params=habitat.parameters, 
                                              model.file="cjs-mnl-TSM-habitat-2groups.bug", n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    }
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # MCMC options for time model - needs more iterations and burnin
    ni <- 20000
    nt <- 3
    nb <- 4000
    nc <- 3
    
    # Time -----------------------------------------------------------------------------------------------------
    write(paste("Time model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    time.data <- list(marr = marr, n.occasions = ncol(marr), r=rowSums(marr))
    time.inits <- function(){list(mean.phi2 = runif(1, 0, 1), 
                                  sigma = runif(1, 0, 5),
                                  mean.p = runif(1, 0, 1))}  
    time.parameters <- c("mean.p", "mean.phi2", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["time"]] <- JAGSParallel(nc, data=time.data, inits=time.inits, params=time.parameters, 
                                       model.file="cjs-mnl-ran-time.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # Time TSM ----------------------------------------------------------------------------------------------
    write(paste("Time TSM model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    time.data <- list(marr.TSM1 = marr.TSM1, marr.TSM2 = marr.TSM2, n.occasions = ncol(marr), r.TSM1=rowSums(marr.TSM1), r.TSM2=rowSums(marr.TSM2))
    time.tsm.inits <- function(){list(mean.phi1 = runif(1, 0, 1), 
                                  mean.phi2 = runif(1, 0, 1), 
                                  sigma = runif(1, 0, 5),
                                  mean.p = runif(1, 0, 1))}  
    time.tsm.parameters <- c("mean.p", "mean.phi1", "mean.phi2", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["time.tsm"]] <- JAGSParallel(nc, data=time.data, inits=time.tsm.inits, params=time.tsm.parameters, 
                                       model.file="cjs-mnl-TSM-ran-time.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # EVI ---------------------------------------------------------------------------------------------------
    write(paste("EVI model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    evi.data <- list(marr = marr, n.occasions = ncol(marr), x = RS.dat$evi, r=rowSums(marr))
    evi.parameters <- c("mean.p", "mean.phi2", "beta_evi", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["evi"]] <- JAGSParallel(nc, data=evi.data, inits=time.inits, params=evi.parameters, 
                                      model.file="cjs-mnl-ran-time-evi.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)

    # EVI TSM ---------------------------------------------------------------------------------------------------
    write(paste("EVI TSM model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    evi.data <- list(marr.TSM1 = marr.TSM1, marr.TSM2 = marr.TSM2, n.occasions = ncol(marr), x = RS.dat$evi, r.TSM1=rowSums(marr.TSM1), r.TSM2=rowSums(marr.TSM2))
    evi.tsm.parameters <- c("mean.p", "mean.phi1", "mean.phi2", "beta_evi", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["evi.tsm"]] <- JAGSParallel(nc, data=evi.data, inits=time.tsm.inits, params=evi.tsm.parameters, 
                                      model.file="cjs-mnl-TSM-ran-time-evi.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # temp --------------------------------------------------------------------------------------------------
    write(paste("Temperature model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    temp.data <- list(marr = marr, n.occasions = ncol(marr), x = RS.dat$temp, r=rowSums(marr))
    temp.parameters <- c("mean.p", "mean.phi2", "beta_temp", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["temp"]] <- JAGSParallel(nc, data=temp.data, inits=time.inits, params=temp.parameters, 
                                       model.file="cjs-mnl-ran-time-temp.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # temp TSM--------------------------------------------------------------------------------------------------
    write(paste("Temperature TSM model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    temp.data <- list(marr.TSM1 = marr.TSM1, marr.TSM2 = marr.TSM2, n.occasions = ncol(marr), x = RS.dat$temp, r.TSM1=rowSums(marr.TSM1), r.TSM2=rowSums(marr.TSM2))
    temp.tsm.parameters <- c("mean.p", "mean.phi1", "mean.phi2", "beta_temp", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["temp.tsm"]] <- JAGSParallel(nc, data=temp.data, inits=time.inits, params=temp.tsm.parameters, 
                                           model.file="cjs-mnl-TSM-ran-time-temp.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # EVI + temp --------------------------------------------------------------------------------------------
    write(paste("EVI + Temperature model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    evitemp.data <- list(marr = marr, n.occasions = ncol(marr), x1 = RS.dat$evi, x2 = RS.dat$temp, r=rowSums(marr))
    evitemp.parameters <- c("mean.p", "mean.phi2", "beta_evi", "beta_temp", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["evi_temp"]] <- JAGSParallel(nc, data=evitemp.data, inits=time.inits, params=evitemp.parameters, 
                                       model.file="D:/ecuador_survival/cjs-mnl-ran-time-2cov.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # EVI + temp TSM--------------------------------------------------------------------------------------------
    write(paste("EVI + Temperature TSM model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    evitemp.data <- list(marr.TSM1 = marr.TSM1, marr.TSM2 = marr.TSM2, n.occasions = ncol(marr), x1 = RS.dat$evi, x2 = RS.dat$temp, r.TSM1=rowSums(marr.TSM1), r.TSM2=rowSums(marr.TSM2))
    evitemp.parameters <- c("mean.p", "mean.phi1", "mean.phi2", "beta_evi", "beta_temp", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["evi_temp.tsm"]] <- JAGSParallel(nc, data=evitemp.data, inits=time.tsm.inits, params=evitemp.parameters, 
                                               model.file="D:/ecuador_survival/cjs-mnl-TSM-ran-time-2cov.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # EVI * temp --------------------------------------------------------------------------------------------
    write(paste("EVI + Temperature model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    evitempint.data <- list(marr = marr, n.occasions = ncol(marr), x1 = RS.dat$evi, x2 = RS.dat$temp, r=rowSums(marr))
    evitempint.parameters <- c("mean.p", "mean.phi2", "beta_evi", "beta_temp", "beta_int", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["evi_temp_int"]] <- JAGSParallel(nc, data=evitempint.data, inits=time.inits, params=evitempint.parameters, 
                                       model.file="cjs-mnl-ran-time-2covint.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    # EVI * temp TSM--------------------------------------------------------------------------------------------
    write(paste("EVI + Temperature TSM model for ", species, sep=" "), logfile.name, append = TRUE)
    write(paste0("ni = ", ni, ", nt = ", nt, " , nb = ", nb, ", nc = ", nc), logfile.name, append = TRUE)
    strt <- Sys.time()       
    evitempint.data <- list(marr.TSM1 = marr.TSM1, marr.TSM2 = marr.TSM2, n.occasions = ncol(marr), x1 = RS.dat$evi, x2 = RS.dat$temp, r.TSM1=rowSums(marr.TSM1), r.TSM2=rowSums(marr.TSM2))
    evitempint.parameters <- c("mean.p", "mean.phi1", "mean.phi2", "beta_evi", "beta_temp", "beta_int", "sigma2", "sigma2.real", "fit", "fit.new")
    modelout[["evi_temp_int.tsm"]] <- JAGSParallel(nc, data=evitempint.data, inits=time.inits, params=evitempint.parameters, 
                                                   model.file="cjs-mnl-TSM-ran-time-2covint.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    write(paste("Model run took", round(Sys.time()-strt, 2),  units(Sys.time()-strt), sep = " "), logfile.name, append = TRUE)
    
    save(modelout, file=paste0("results/", species, "_CJS_nojuv_model_output.rda"))
    
}

# runs the code to check model convergence (output is plots of model convergence and the max rhat value)
source("code/check_convergence.R")
