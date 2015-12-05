JAGSParallelized = function(Index, jags.seed.vec, inits, model.file, data, params, n.thin, n.iter, n.burnin, DIC){
    RNGset = c("Mersenne-Twister","Marsaglia-Multicarry","Super-Duper","Knuth-TAOCP-2002","Knuth-TAOCP","Wichmann-Hill","L'Ecuyer-CMRG")
    set.seed(1,kind=RNGset[Index])
    random.seed = runif(1,1,1e6)
    Jags = jags(inits=inits, n.chains=1, model.file=model.file, data=data, parameters.to.save=params, n.thin=n.thin, n.iter=n.iter, n.burnin=n.burnin, DIC=DIC)
    return(Jags)
}

#jags.seed=random.seed,

JAGSParallel = function(n.cores,data,inits,params,model.file,n.chains,n.iter,n.burnin,n.thin,DIC=FALSE){
    # Start snowfall
    sfInit(parallel=TRUE, cpus=n.cores)
    sfLibrary(R2jags)
    sfLibrary(snowfall)
    sfExportAll()
    sfClusterCall( runif, 4 )
    jags.seed.vec = ceiling(runif(n.chains,1,1e6))
    # Run JAGS
    JAGSList <- sfLapply(1:n.chains, JAGSParallelized, jags.seed.vec=jags.seed.vec, data=data, inits=inits, params=params,
                         model.file= model.file, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC=DIC)
    # End snowfall
    sfStop()
    result <- NULL
    model <- NULL
    for (ch in 1:n.chains) {
        result <- abind(result, JAGSList[[ch]]$BUGSoutput$sims.array, along = 2)
        model[[ch]] <- JAGSList[[ch]]$model
    }
    result <- as.bugs.array(result, model.file = model.file,program = "jags", DIC = DIC, n.iter = n.iter, n.burnin = n.burnin,n.thin = n.thin)
    out <- list(model = model, JAGSoutput = result, parameters.to.save = params, model.file = model.file, n.iter = n.iter, DIC = DIC)
    class(out) <- c("rjags.parallel", "rjags")
    return(out)
}