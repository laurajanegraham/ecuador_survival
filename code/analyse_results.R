# load packages
library(ggplot2)
library(R2WinBUGS)
library(dplyr)
library(tidyr)
library(taxize)
library(mcmcplots)
library(stringr)
library(Hmisc)
library(cowplot)

source("code/fnCleanBandingDat.R")

banding.dat.clean <- CleanBandingDat()
sphab <- unique(banding.dat.clean[,c('Specie.Name', 'habitat')])

jags2plot <- function(x) {
    load(x)
    species <- gsub("_", " ", str_match(x,pattern="results/(\\w+.\\w+)_CJS")[,2])
    phi.null <- data.frame(modelout$null$JAGSoutput$summary[c('mean.p', 'mean.phi'),c('mean', '2.5%','97.5%')])
    phi.hab <- data.frame(modelout$habitat$JAGSoutput$summary[-(1:2),c('mean', '2.5%','97.5%')])
    # need to account for species that haven't been banded in all habitats - done by adding in a row of NAs for the missing habitat
    if(nrow(phi.hab)<=3) {
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
        if(nrow(phi.hab)==0) {
          phi.hab <- data.frame(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA))
          names(phi.hab) <- names(phi.null)
        }
    }
    phi.time <- data.frame(modelout$time$JAGSoutput$summary[c('mean.p', 'mean.phi', 'sigma2'),c('mean', '2.5%','97.5%')])
    phi.evi <- data.frame(modelout$evi$JAGSoutput$summary[c('mean.p', 'mean.phi', 'sigma2', 'beta'),c('mean', '2.5%','97.5%')])
    phi.temp <- data.frame(modelout$temp$JAGSoutput$summary[c('mean.p', 'mean.phi', 'sigma2', 'beta'),c('mean', '2.5%','97.5%')])
    phi.evi_temp <- data.frame(modelout$evi_temp$JAGSout$summary[c('mean.p', 'mean.phi', 'sigma2', 'beta1', 'beta2'),c('mean', '2.5%','97.5%')])
    phi.evi_temp_int <- data.frame(modelout$evi_temp_int$JAGSoutput$summary[c('mean.p', 'mean.phi', 'sigma2', 'beta1', 'beta2', 'beta3'),c('mean', '2.5%','97.5%')])
    phi.null$model <- "Null"
    phi.hab$model <- "Habitat"
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
    fitnull <- data.frame(fit=modelout$null$JAGSoutput$sims.list$fit, fit.new=modelout$null$JAGSoutput$sims.list$fit.new) %>%
        summarise(p.val = round(mean(fit.new > fit), 2)) %>%
        mutate(model="Null")
    fithab <- data.frame(fit=modelout$habitat$JAGSoutput$sims.list$fit, fit.new=modelout$habitat$JAGSoutput$sims.list$fit.new) %>%
        summarise(p.val = round(mean(fit.new > fit), 2)) %>%
        mutate(model="Habitat")
    fittime <- data.frame(fit=modelout$time$JAGSoutput$sims.list$fit, fit.new=modelout$time$JAGSoutput$sims.list$fit.new) %>%
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

mod <- getRes("CJS_nojuv")

save_plot(mod$plots, filename = "figures/adult_survival.png", base_height = 12, base_width = 20)

fit <- spread(mod$fit, model, p.val)

fit$bestmod <- apply(fit[2:4], 1, function(x) {
    x <- 0.5-x
    names(which.min(abs(x)))
})


write.csv(fit, file="results/model_comparison.csv")

write.csv(mod$dat[which(complete.cases(mod$dat)),], file="results/nojuv_raw_results.csv")
write.csv(usable_dat, file="results/nojuv_recapture_summaries.csv")

# Models with time varying covariates
