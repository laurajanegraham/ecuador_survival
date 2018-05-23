# load packages
library(ggplot2)
library(R2WinBUGS)
library(plyr)
library(dplyr)
library(tidyr)
library(mcmcplots)
library(stringr)
library(Hmisc)
library(cowplot)


source("code/fnCleanBandingDat.R")

# we have done two levels of analysis: one at the species level, and one at the forest specialisation level
#analysis_level <- "specialisms"
analysis_level <- "species"

data_input <- c("data/banding_sheet.csv", "data/2016_data.csv")
banding.dat.clean <- ldply(data_input, function(f) CleanBandingDat(f, inc.juvenile = FALSE))
guilds <- read.csv("data/forest_specialization.csv")
guilds$Specie.Name <- toupper(guilds$Species)


if(analysis_level == "specialisms") { 
  banding.dat.clean <- merge(banding.dat.clean, guilds)
  banding.dat.clean$tax_level <- banding.dat.clean$Forest.specialization
} else if(analysis_level == "species") {
  banding.dat.clean <- merge(banding.dat.clean, guilds, all.x = TRUE)
  banding.dat.clean$tax_level <- banding.dat.clean$Specie.Name
} else {
  stop("Make sure an analysis level is included at the top of the script!")
}

sphab <- unique(banding.dat.clean[,c('tax_level', 'habitat')]) %>%
    arrange(tax_level, habitat) %>%
    group_by(tax_level) %>%
    mutate(newhab=habitat,
           habitat=paste0("hab",1:n()),
           tax_identity=tax_level,
           count=n()) %>%
    mutate(habitat = ifelse(count==3 & newhab=="Introduced", "intro", habitat),
           habitat = ifelse(count==3 & newhab=="Native", "native", habitat),
           habitat = ifelse(count==3 & newhab=="Scrub", "scrub", habitat))

getRes <- function(x) {
    load(x)
    res <- lapply(names(modelout), function(x) {
        dat <- data.frame(modelout[[x]]$JAGSoutput$summary)
        dat$model <- x
        dat$param <- rownames(dat)
        return(dat)
    })
    
    res <- do.call("rbind", res)
    res <- filter(res, !param %in% c("fit", "fit.new", "sigma2.real", c(paste0("phi[",1:29,"]"), paste0("phi.TSM1[",1:29,"]"), paste0("phi.TSM2[",1:29,"]"))))
    tax_identity <- gsub("_", " ", str_match(x,pattern="/(\\w+.\\w+)_CJS")[,2])
    res$tax_identity <- tax_identity
    return(res)
}

getFit <- function(x) {
  load(x)
  tax_identity <- gsub("_", " ", str_match(x,pattern="/(\\w+.\\w+)_CJS")[,2])
  fit <- lapply(modelout, function(x) {
    fit <- data.frame(fit=x$JAGSoutput$sims.list$fit, fit.new=x$JAGSoutput$sims.list$fit.new) %>%
      summarise(p.val = round(mean(fit.new > fit), 2))
  })
  fit <- do.call("rbind", fit)
  fit <- data.frame(t(fit))
  rownames(fit) <- tax_identity
  # not all species could fit the habitat model (only present in one habitat) 
  if(!"habitat" %in% colnames(fit)) fit$habitat <- NA
  if(!"habitat.tsm" %in% colnames(fit)) fit$habitat.tsm <- NA
  return(fit)
}


files <- list.files(paste0("results/", analysis_level), pattern="CJS_nojuv", full.names = TRUE)  

res.full <- lapply(files, getRes)
res.full <- do.call("rbind", res.full) %>% 
  mutate(param=gsub("mean.", "", param)) %>%
  separate(param, c("param", "habitat"), sep="[.]", fill="right") %>%
  left_join(sphab) %>%
  select(-tax_level, -habitat, -count) %>%
  write_csv(paste0("results/", analysis_level, "/nojuv_raw_results.csv"))

fit <- lapply(files, getFit)
fit <- do.call("rbind", fit)

# we only want to know between null/habitat
fit$bestmod <- apply(fit[,c("null", "null.tsm", "habitat", "habitat.tsm")], 1, function(x) {
    x <- 0.5-x
    names(which.min(abs(x)))
})

fit$tax_identity <- rownames(fit)

# get only results for best fitting model for table 1
res <- merge(res.full, fit, by.x=c("tax_identity", "model"), by.y=c("tax_identity", "bestmod"))

# rename the null TSM model to make column separation easier
res <- separate(res, model, c("model", "tsm"), sep="[.]", fill="right")

# create results table so that it requires minimal adjustment in excel
res_mod <- mutate(res, val=paste0(sprintf("%.2f", round(mean,2)), " (", sprintf("%.2f", round(X2.5.,2)), "-", sprintf("%.2f", round(X97.5.,2)), ")"),
                  param=gsub("mean.", "", param)) %>%
    select(tax_identity, param, val, model) %>%
    separate(param, c("param", "habitat"), sep="[.]", fill="right") %>%
    spread(param, val) %>%
    arrange(model, tax_identity, habitat) %>%
    merge(sphab, all.x=TRUE) %>%
    mutate(habitat = newhab) %>%
    arrange(model, tax_identity, habitat)

# create the recaps summaries
dat_summary <- group_by(banding.dat.clean, tax_level, Band.Number) %>%
  summarise(recap=n()-1)

dat_recaps <- mutate(dat_summary, 
                     N=1) %>%
  spread(recap, N) %>%
  select(-Band.Number) %>%
  group_by(tax_level) %>%
  summarise_all(funs(sum(., na.rm=TRUE))) %>%
  mutate(nind=rowSums(.[-1]))

dat_nind <- group_by(dat_summary, tax_level) %>%
  summarise(nind=n())

dat_recaps$recaps <- rowSums(dat_recaps[,-c(1, 2, ncol(dat_recaps))])

# output results
write_csv(res_mod, paste0("results/", analysis_level, "/model_results.csv"))

write_csv(fit, paste0("results/", analysis_level, "/model_comparison.csv"))

write_csv(dat_recaps, paste0("results/", analysis_level, "/nojuv_recapture_summaries.csv"))

#### 
# BELOW HERE IS RUBBISH
####

# Models with time varying covariates - this is rubbish because nearly all slope CIs contain 0. 
res.full <- read.csv("results/nojuv_raw_results.csv")
fit <- read.csv("results/model_comparison.csv")

time_species <- filter(fit, bestmod %in% c("time", "time.tsm")) %>%
    select(species, bestmod) %>%
    separate(bestmod, c("time", "tsm"), sep="[.]", fill="right")

sigmanull <- filter(res.full, model=="time", param=="sigma2") %>%
    mutate(sigma2null=mean,
           sigma2min=X2.5.,
           sigma2max=X97.5.) %>%
    select(species, tsm, sigma2null, sigma2min, sigma2max) 

    sigmaother <- filter(res.full, param=="sigma2", !model %in% c("null", "habitat", "time")) %>%
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
    merge(sigmaother)
    




if(bestmod == "tsm") sigma <- filter(sigma, !is.na(tsm))
if(bestmod != "tsm") sigma <- filter(sigma, is.na(tsm))

sigma.null <- filter(sigma, model=="time")
sigma.other <- filter(sigma, model!="time")

sigma.other <- group_by(sigma.other, species) %>%
    mutate(rsquare = (mean - sigma.null$mean) / sigma.null$mean)

# Plot of through time estimates of survival 
time.species <- filter(fit, bestmod %in% c("time", "time.tsm")) %>%
  select(species, bestmod)

get.time.res <- function(time.species, x) { 
  species <- time.species[x, "species"]
  bestmod <- time.species[x, "bestmod"]
  file <- list.files("results", pattern = species, full.names = TRUE)
  load(file)
  mod <- modelout[[bestmod]]
  dat <- data.frame(mod$JAGSoutput$summary)
  dat$param <- rownames(dat)
  dat <- filter(dat, param %in% c(paste0("phi[",1:29,"]"), paste0("phi.TSM2[",1:29,"]")))
  dat$session <- sort(unique(banding.dat.clean$session_new))[-1]
  dat$species <- species
  return(dat)
}

time.res <- lapply(rownames(time.species), function(x) get.time.res(time.species, x))
time.res <- do.call("rbind", time.res)

ggplot(time.res, aes(x = session, y = mean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), width=.1, position = position_dodge(0.1)) + 
  geom_line(position = position_dodge(0.1)) +
  facet_wrap(~species, ncol = 2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

