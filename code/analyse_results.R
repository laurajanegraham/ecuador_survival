# load packages
library(ggplot2)
library(R2WinBUGS)
library(dplyr)
library(tidyr)
library(mcmcplots)
library(stringr)
library(Hmisc)
library(cowplot)

source("code/fnCleanBandingDat.R")

banding.dat.clean <- CleanBandingDat()
sphab <- unique(banding.dat.clean[,c('Specie.Name', 'habitat')])

getRes <- function(x) {
    load(x)
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
    return(res)
}

getFit <- function(x) {
    load(x)
    species <- gsub("_", " ", str_match(x,pattern="results/(\\w+.\\w+)_CJS")[,2])
    fit <- lapply(modelout, function(x) {
        fit <- data.frame(fit=x$JAGSoutput$sims.list$fit, fit.new=x$JAGSoutput$sims.list$fit.new) %>%
            summarise(p.val = round(mean(fit.new > fit), 2))
        })
    fit <- do.call("rbind", fit)
    fit <- data.frame(t(fit))
    rownames(fit) <- species
    # not all species could fit the habitat model (only present in one habitat) 
    if(!"habitat" %in% colnames(fit)) fit$habitat <- NA
    if(!"habitat.tsm" %in% colnames(fit)) fit$habitat.tsm <- NA
    return(fit)
}

files <- list.files("results", pattern="CJS_nojuv", full.names = TRUE)

res <- lapply(files, getRes)
res <- do.call("rbind", res)
write.csv(res[which(complete.cases(res)),], file="results/nojuv_raw_results.csv")

fit <- lapply(files, getFit)
fit <- do.call("rbind", fit)
fit$bestmod <- apply(fit, 1, function(x) {
    x <- 0.5-x
    names(which.min(abs(x)))
})

fit$species <- rownames(fit)

write.csv(fit, file="results/model_comparison.csv")

write.csv(usable_dat, file="results/nojuv_recapture_summaries.csv")

# Models with time varying covariates
sigma <- filter(res.el, param=="sigma2") %>%
    separate(model, c("model", "tsm"), sep="[.]", fill="right")

if(bestmod == "tsm") sigma <- filter(sigma, !is.na(tsm))
if(bestmod != "tsm") sigma <- filter(sigma, is.na(tsm))

sigma.null <- filter(sigma, model=="time")
sigma.other <- filter(sigma, model!="time")

sigma.other <- group_by(sigma.other, species) %>%
    mutate(rsquare = (mean - sigma.null$mean) / sigma.null$mean)
