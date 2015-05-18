require(R2WinBUGS)
require(mcmcplots)
source("code/fnCleanBandingDat.R")
source("code/fnEncounterHistory.R")
source("code/fnKnownStateInits.R")

# Load and format banding data -------------------------------------------------
banding.dat.clean <- CleanBandingDat()

sp_dat <- filter(banding.dat.clean, Specie.Name == species) %>%
        select(Band.Number, session_new, habitat)

sp_eh <- EncounterHistory(sp_dat, "session_new", "Band.Number", "habitat")

marr <- sp_eh$m.array
marr.gp <- sp_eh$m.array.gp

bugs.data <- list(marr = marr, n.occasions = ncol(marr))

inits <- function(){list(mean.phi = runif(1, 0, 1), sigma = runif(1, 0, 5), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("phi", "mean.p", "mean.phi", "sigma2", "sigma2.real", "fit", "fit.new")

ni <- 5000
nt <- 3
nb <- 1000
nc <- 3

cjs <- bugs(bugs.data, inits, parameters, "cjs-mnl-ran-time.bug", n.chains = nc, n.thin = nt, 
            n.iter = ni, n.burnin = nb, debug = TRUE, 
            working.directory='~/.wine/drive_c/temp/Rtmp/', clearWD=TRUE)
