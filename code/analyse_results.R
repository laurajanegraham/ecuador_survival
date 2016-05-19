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
sphab <- unique(banding.dat.clean[,c('Specie.Name', 'habitat')]) %>%
    arrange(Specie.Name, habitat) %>%
    group_by(Specie.Name) %>%
    mutate(newhab=habitat,
           habitat=paste0("hab",row_number()),
           species=Specie.Name,
           count=n()) %>%
    mutate(habitat = ifelse(count==3 & newhab=="Introduced", "intro", habitat),
           habitat = ifelse(count==3 & newhab=="Native", "native", habitat),
           habitat = ifelse(count==3 & newhab=="Scrub", "scrub", habitat))

sphab_all <- filter(sphab)

getRes <- function(x) {
    load(x)
    res <- lapply(names(modelout), function(x) {
        dat <- data.frame(modelout[[x]]$JAGSoutput$summary)
        dat$model <- x
        dat$param <- rownames(dat)
        return(dat)
    })
    
    res <- do.call("rbind", res)
    res <- filter(res, !param %in% c("fit", "fit.new", "sigma2.real", paste0("phi[",1:29,"]")))
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

res.full <- lapply(files, getRes)
res.full <- do.call("rbind", res.full)
write.csv(res.full[which(complete.cases(res.full)),], file="results/nojuv_raw_results.csv")

fit <- lapply(files, getFit)
fit <- do.call("rbind", fit)

# we only want to know between null/habitat/time
fit$bestmod <- apply(fit[1:6], 1, function(x) {
    x <- 0.5-x
    names(which.min(abs(x)))
})

fit$species <- rownames(fit)

# get only results for best fitting model for table 1
res <- merge(res.full, fit, by.x=c("species", "model"), by.y=c("species", "bestmod"))

# rename the null TSM model to make column separation easier
res$model <- ifelse(res$model=="TSM", "null.tsm", as.character(res$model))

# create results table so that it requires minimal adjustment in excel
res_mod <- mutate(res, val=paste0(sprintf("%.2f", round(mean,2)), " (", sprintf("%.2f", round(X2.5.,2)), "-", sprintf("%.2f", round(X97.5.,2)), ")"),
                  param=gsub("mean.", "", param)) %>%
    select(species, param, val, model) %>%
    separate(param, c("param", "habitat"), sep="[.]", fill="right") %>%
    spread(param, val) %>%
    arrange(model, species, habitat) %>%
    merge(sphab, all.x=TRUE) %>%
    mutate(habitat = newhab) %>%
    arrange(model, species, habitat) %>%
    select(model, species, p, habitat, phi1, phi2, sigma2)

write.csv(res_mod, file="results/model_results.csv")

write.csv(fit, file="results/model_comparison.csv")

write.csv(usable_dat, file="results/nojuv_recapture_summaries.csv")

# Models with time varying covariates
time_species <- filter(fit, bestmod %in% c("time", "time.tsm")) %>%
    select(species, bestmod) %>%
    separate(bestmod, c("time", "tsm"), sep="[.]", fill="right")

res.full <- separate(res.full, model, c("model", "tsm"), sep="[.]", fill="right")

sigmanull <- filter(res.full, model=="time", param=="sigma2.real") %>%
    mutate(sigma2null=mean,
           sigma2min=X2.5.,
           sigma2max=X97.5.) %>%
    select(species, tsm, sigma2null, sigma2min, sigma2max) 

    sigmaother <- filter(res.full, param=="sigma2.real", !model %in% c("null", "habitat", "time")) %>%
        merge(time_species, by=c("species", "tsm")) %>%
        merge(sigmanull) %>%
        mutate(rsq=(sigma2null-mean)/sigma2null,
               rsqmin=(sigma2min-X2.5.)/sigma2min,
               rsqmax=(sigma2max-X97.5.)/sigma2max)


rs.mods <- filter(res.full, !model %in% c("null", "habitat", "time")) %>%
    merge(time_species, by=c("species", "tsm")) %>%
    mutate(val=paste0(sprintf("%.2f", round(mean,2)), " (", sprintf("%.2f", round(X2.5.,2)), "-", sprintf("%.2f", round(X97.5.,2)), ")"),
           param=gsub("mean.", "", param)) %>%
    select(species, param, val, model, tsm) %>%
    spread(param, val) %>%
    merge(sigmanull) %>%
    mutate(rsq=(sigma2null-as.numeric(substr(sigma2, 1, 4)))/sigma2null,
           rsqmin=(sigma2min-as.numeric(substr(sigma2, 1, 4)))/sigma2null,)
    




if(bestmod == "tsm") sigma <- filter(sigma, !is.na(tsm))
if(bestmod != "tsm") sigma <- filter(sigma, is.na(tsm))

sigma.null <- filter(sigma, model=="time")
sigma.other <- filter(sigma, model!="time")

sigma.other <- group_by(sigma.other, species) %>%
    mutate(rsquare = (mean - sigma.null$mean) / sigma.null$mean)
