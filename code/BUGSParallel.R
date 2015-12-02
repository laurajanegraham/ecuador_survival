# BUGSParallel.R ---------------------------------------------------------------
# Laura Graham, adapted from Heather Lynch

# Allows for MCMC chains to be run in parallel


BUGSParallelized = function(Index, bugs.seed.vec, inits, model.file, data, working.directory=NULL, params, n.thin, n.iter, n.burnin, DIC, bugs.directory){
    RNGset = c("Mersenne-Twister","Marsaglia-Multicarry","Super-Duper","Knuth-TAOCP-2002","Knuth-TAOCP","Wichmann-Hill","L'Ecuyer-CMRG")
    set.seed(1,kind=RNGset[Index])
    random.seed = runif(1,1,1e6)
    mod = bugs(inits=inits, n.chains=1, model.file=model.file, working.directory=working.directory, data=data, parameters.to.save=params, n.thin=n.thin, n.iter=n.iter, n.burnin=n.burnin, DIC=DIC, clearWD = TRUE, bugs.directory = bugs.directory)
    return(mod)
}

BUGSParallel = function(n.cores,data,inits,params,model.file,debug,n.chains,n.iter,n.burnin,n.thin,DIC=FALSE,working.directory=NULL, bugs.directory=bugs.directory){
    # Start snowfall
    sfInit(parallel=TRUE, cpus=n.cores)
    sfLibrary(R2WinBUGS)
    sfLibrary(snowfall)
    sfExportAll()
    sfClusterCall( runif, 4 )
    bugs.seed.vec = ceiling(runif(n.chains,1,1e6))
    # Run BUGS
    BUGSList <- sfLapply(1:n.chains, BUGSParallelized, bugs.seed.vec =bugs.seed.vec, data=data, inits=inits, params=params,
                         model.file= model.file, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC=DIC, working.directory=working.directory, bugs.directory=bugs.directory)
# End snowfall
    sfStop()
    result <- NULL
    model <- NULL
    for (ch in 1:n.chains) {
        result <- abind(result, BUGSList[[ch]]$sims.array, along = 2)
        model[[ch]] <- BUGSList[[ch]]$model
    }
    result <- as.bugs.array(result, model.file = model.file,program = "bugs", DIC = DIC, n.iter = n.iter, n.burnin = n.burnin,n.thin = n.thin)
    out <- list(model = model, BUGSoutput = result, parameters.to.save = params, model.file = model.file, n.iter = n.iter, DIC = DIC)
    class(out) <- c("rbugs.parallel", "rbugs")
    return(out)
}
