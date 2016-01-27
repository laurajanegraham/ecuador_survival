# load packages
require(ggplot2)
require(R2WinBUGS)
require(dplyr)
require(tidyr)
require(taxize)
require(mcmcplots)
require(stringr)
require(Hmisc)

source("code/fnCleanBandingDat.R")

banding.dat.clean <- CleanBandingDat()
sphab <- unique(banding.dat.clean[,c('Specie.Name', 'habitat')])

jags2plot <- function(x) {
    load(x)
    species <- gsub("_", " ", str_match(x,pattern="results/(\\w+.\\w+)_CJS")[,2])
    phi.null <- data.frame(modelout[[1]]$JAGSoutput$summary[c('mean.p', 'mean.phi'),c('mean', '2.5%','97.5%')])
    phi.null$model <- "Null"
    phi.hab <- data.frame(modelout[[2]]$JAGSoutput$summary[-(1:2),c('mean', '2.5%','97.5%')])
    # need to account for species that haven't been banded in all habitats - done by adding in a row of NAs for the missing habitat
    if(nrow(phi.hab)==3) {
        habitats <- filter(sphab, Specie.Name==species) %>% select(habitat) %>% unlist(.)
        if("Introduced" %in% habitats && "Native" %in% habitats) {
            phi.hab <- rbind(phi.hab, c(NA,NA,NA))
        }
        if("Introduced" %in% habitats && "Scrub" %in% habitats) {
            phi.hab <- rbind(phi.hab[1:2,], c(NA, NA, NA), phi.hab[3,])
        }
        if("Native" %in% habitats && "Scrub" %in% habitats) {
            phi.hab <- rbind(phi.hab[1,], c(NA, NA, NA), phi.hab[2:3,])
        }
    }
    phi.hab$model <- "Habitat"
    phi.time <- data.frame(modelout[[3]]$JAGSoutput$summary[c('mean.p', 'mean.phi', 'sigma2.real'),c('mean', '2.5%','97.5%')])
    phi.time$model <- "Time"
    param.est <- rbind(phi.null,phi.hab,phi.time)
    colnames(param.est) <- c("mean", "lci", "uci", "model")
    
    param.est$param <- factor(c('pnull', 'phinull', 'phab', 'phiintro', 'phinative', 'phishrub', 'ptime', 'phitime', 'sigma2'),
                              levels=c('pnull', 'phinull', 'phab', 'phiintro', 'phinative', 'phishrub', 'ptime', 'phitime', 'sigma2'))
    
    param.est$Model <- factor(param.est$model, levels=c('Null', 'Habitat', 'Time'))
    param.est$species <- species
    return(param.est)
}

getFit <- function(x) {
    load(x)
    species <- gsub("_", " ", str_match(x,pattern="results/(\\w+.\\w+)_CJS")[,2])
    fitnull <- data.frame(fit=modelout[[1]]$JAGSoutput$sims.list$fit, fit.new=modelout[[1]]$JAGSoutput$sims.list$fit.new) %>%
        summarise(p.val = round(mean(fit.new > fit), 2)) %>%
        mutate(model="Null")
    fithab <- data.frame(fit=modelout[[2]]$JAGSoutput$sims.list$fit, fit.new=modelout[[2]]$JAGSoutput$sims.list$fit.new) %>%
        summarise(p.val = round(mean(fit.new > fit), 2)) %>%
        mutate(model="Habitat")
    fittime <- data.frame(fit=modelout[[3]]$JAGSoutput$sims.list$fit, fit.new=modelout[[3]]$JAGSoutput$sims.list$fit.new) %>%
        summarise(p.val = round(mean(fit.new > fit), 2)) %>%
        mutate(model="Time")
    fit <- rbind(fitnull, fithab, fittime) %>%
        mutate(species=species)
}

# list results files
getRes <- function(fname) {
    files <- list.files("results", pattern=fname, full.names = TRUE)
    
    dat <- lapply(files, jags2plot)
    dat <- do.call("rbind", dat)
    
    fit <- lapply(files, getFit)
    fit <- do.call("rbind", fit)
    
    plots <- ggplot(dat, aes(x=param, y=mean, colour=Model)) + 
        geom_point() + 
        geom_errorbar(aes(ymin=lci, ymax=uci), width=0.1) + 
        facet_wrap(~species, ncol=3) +
        xlab(expression("Parameter")) + ylab(expression(paste("Mean"%+-%"95% Credible Interval"))) + 
        scale_x_discrete(breaks=c('pnull', 'phinull', 'phab', 'phiintro', 'phinative', 'phishrub', 'ptime', 'phitime', 'sigma2'),
                         labels=c(expression('p'['null']), expression(phi['null']), 
                                  expression('p'['habitat']), expression(phi['introduced']), expression(phi['native']), expression(phi['shrub']),
                                  expression('p'['time']), expression(phi['time'], sigma^2)))
    
    return(list(dat=dat, fit=fit, plots=plots))
}

fullmod <- getRes("CJS_full")
adultmod <- getRes("CJS_nojuv")

save_plot(adultmod$plots, filename = "adult_survival.png", base_height = 12, base_width = 20)
save_plot(fullmod$plots, filename = "all_survival.png", base_height = 12, base_width = 20)

fullfit <- spread(fullmod$fit, model, p.val)
adultfit <- spread(adultmod$fit, model, p.val)

comp <- merge(fullfit, adultfit, by="species", all=TRUE)
comp$bestmod <- apply(comp[2:7], 1, function(x) {
    x <- 0.5-x
    names(which.min(abs(x)))
})


write.csv(comp, file="results/full_comparison.csv")


