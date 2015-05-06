require(R2WinBUGS)
source("code/fnCleanBandingDat.R")
source("code/fnEncounterHistory.R")

# Load and format banding data -------------------------------------------------
banding.dat.clean <- CleanBandingDat()

# initially going to be testing on the C. iris records
species <- "COELIGENA IRIS"

# Create group variable
group <- c(rep(1, dim(CH.F)[1]), rep(2, dim(CH.M)[1]))

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language


# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g = length(unique(group)), group = group)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), phi.g = runif(length(unique(group)), 0, 1), p.g = runif(length(unique(group)), 0, 1))}  

# Parameters monitored
parameters <- c("phi.g", "p.g")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call WinBUGS from R (BRT 2 min)
cjs.group <- bugs(bugs.data, inits, parameters, "cjs-group.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, working.directory='~/.wine/drive_c/temp/Rtmp/', clearWD=TRUE)

# Summarize posteriors
print(cjs.group, digits = 3)
