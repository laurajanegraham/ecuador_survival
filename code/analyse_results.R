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

getRes <- function(x) {
    load(x)
    modelout <- modelout[which(names(modelout) != "habitat.tsm")]
    res <- lapply(names(modelout), function(x) {
        dat <- data.frame(modelout[[x]]$JAGSoutput$summary)
        dat$model <- x
        dat$param <- rownames(dat)
        return(dat)
    })
    
    res <- do.call("rbind", res)
    res <- filter(res, !param %in% c("fit", "fit.new", paste0("phi[",1:29,"]")))
    species <- gsub("_", " ", str_match(x,pattern="results/(\\w+.\\w+)_CJS")[,2])
    res$species <- species
    # 
    # phi.null <- data.frame(modelout$null$JAGSoutput$summary[c('mean.p', 'mean.phi'),c('mean', '2.5%','97.5%')])
    # phi.hab <- data.frame(modelout$habitat$JAGSoutput$summary[-(1:2),c('mean', '2.5%','97.5%')])
    # # need to account for species that haven't been banded in all habitats - done by adding in a row of NAs for the missing habitat
    # if(nrow(phi.hab)<=3) {
    #     habitats <- filter(sphab, Specie.Name==species) %>% select(habitat) %>% unlist(.)
    #     if("Introduced" %in% habitats && "Native" %in% habitats) {
    #         phi.hab <- rbind(phi.hab, c(NA,NA,NA))
    #     }
    #     if("Introduced" %in% habitats && "Scrub" %in% habitats) {
    #         phi.hab <- rbind(phi.hab[1:2,], c(NA, NA, NA), phi.hab[3,])
    #     }
    #     if("Native" %in% habitats && "Scrub" %in% habitats) {
    #         phi.hab <- rbind(phi.hab[1,], c(NA, NA, NA), phi.hab[2:3,])
    #     }
    #     if(nrow(phi.hab)==0) {
    #       phi.hab <- data.frame(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA))
    #       names(phi.hab) <- names(phi.null)
    #     }
    # }
    # phi.time <- data.frame(modelout$time$JAGSoutput$summary[c('mean.p', 'mean.phi', 'sigma2'),c('mean', '2.5%','97.5%')])
    # phi.evi <- data.frame(modelout$evi$JAGSoutput$summary[c('mean.p', 'mean.phi', 'sigma2', 'beta'),c('mean', '2.5%','97.5%')])
    # phi.temp <- data.frame(modelout$temp$JAGSoutput$summary[c('mean.p', 'mean.phi', 'sigma2', 'beta'),c('mean', '2.5%','97.5%')])
    # phi.evi_temp <- data.frame(modelout$evi_temp$JAGSout$summary[c('mean.p', 'mean.phi', 'sigma2', 'beta1', 'beta2'),c('mean', '2.5%','97.5%')])
    # phi.evi_temp_int <- data.frame(modelout$evi_temp_int$JAGSoutput$summary[c('mean.p', 'mean.phi', 'sigma2', 'beta1', 'beta2', 'beta3'),c('mean', '2.5%','97.5%')])
    # phi.null$model <- "Null"
    # phi.hab$model <- "Habitat"
    # phi.time$model <- "Time"
    # param.est <- rbind(phi.null,phi.hab,phi.time)
    # colnames(param.est) <- c("mean", "lci", "uci", "model")
    # 
    # param.est$param <- factor(c('pnull', 'phinull', 'phab', 'phiintro', 'phinative', 'phishrub', 'ptime', 'phitime', 'sigma2'),
    #                           levels=c('pnull', 'phinull', 'phab', 'phiintro', 'phinative', 'phishrub', 'ptime', 'phitime', 'sigma2'))
    # 
    # param.est$Model <- factor(param.est$model, levels=c('Null', 'Habitat', 'Time'))
    # param.est$species <- species
    # return(param.est)
    return(res)
}

getFit <- function(x) {
    load(x)
    # this is just while habitat.tsm is broken - make sure to remove when fixed
    modelout <- modelout[which(names(modelout) != "habitat.tsm")]
    species <- gsub("_", " ", str_match(x,pattern="results/(\\w+.\\w+)_CJS")[,2])
    fit <- lapply(modelout, function(x) {
        fit <- data.frame(fit=x$JAGSoutput$sims.list$fit, fit.new=x$JAGSoutput$sims.list$fit.new) %>%
            summarise(p.val = round(mean(fit.new > fit), 2))
        })
    fit <- do.call("rbind", fit)
    fit <- data.frame(t(fit))
    rownames(fit) <- species
    # not all species could fit the habitat model (only present in one habitat) - will need to add another line for habitat.tsm once that is working.
    if(!"habitat" %in% colnames(fit)) fit$habitat <- NA
    return(fit)
}


# list results files
# getRes <- function(fname) {
#     files <- list.files("results", pattern=fname, full.names = TRUE)
#     dat <- lapply(files, jags2plot)
#     dat <- do.call("rbind", dat)
#     
#     fit <- lapply(files, getFit)
#     fit <- do.call("rbind", fit)
#     
#     plots <- ggplot(dat, aes(x=param, y=mean, colour=Model)) + 
#         geom_point() + 
#         geom_errorbar(aes(ymin=lci, ymax=uci), width=0.1) + 
#         facet_wrap(~species, ncol=3) +
#         xlab(expression("Parameter")) + ylab(expression(paste("Mean"%+-%"95% Credible Interval"))) + 
#         scale_x_discrete(breaks=c('pnull', 'phinull', 'phab', 'phiintro', 'phinative', 'phishrub', 'ptime', 'phitime', 'sigma2'),
#                          labels=c(expression('p'['null']), expression(phi['null']), 
#                                   expression('p'['habitat']), expression(phi['introduced']), expression(phi['native']), expression(phi['shrub']),
#                                   expression('p'['time']), expression(phi['time'], sigma^2)))
#     
#     return(list(dat=dat, fit=fit, plots=plots))
# }

files <- list.files("results", pattern="CJS_nojuv", full.names = TRUE)

#save_plot(mod$plots, filename = "figures/adult_survival.png", base_height = 12, base_width = 20)
res <- lapply(files, getRes)
res <- do.call("rbind", res)

fit <- lapply(files, getFit)
fit <- do.call("rbind", fit)
fit$bestmod <- apply(fit, 1, function(x) {
    x <- 0.5-x
    names(which.min(abs(x)))
})

fit$species <- rownames(fit)

bestmodstats <- merge(fit, res, by.x=c("bestmod", "species"), by.y=c("model", "species")) %>%
    mutate(res, val=paste0(sprintf("%.2f", round(mean,2)), " (", sprintf("%.2f", round(lci,2)), "-", sprintf("%.2f", round(uci,2)), ")")) %>%
    select(species, param, val) %>%
    spread(param, val) 
write.csv(fit, file="results/model_comparison.csv")

write.csv(mod$dat[which(complete.cases(mod$dat)),], file="results/nojuv_raw_results.csv")
write.csv(usable_dat, file="results/nojuv_recapture_summaries.csv")

# Models with time varying covariates



lapply(fit, function(x) length(x))
