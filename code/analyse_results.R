# load packages
require(ggplot2)
require(R2WinBUGS)
require(dplyr)
require(taxize)
require(mcmcplots)
require(stringr)
source("code/fnCleanBandingDat.R")

banding.dat.clean <- CleanBandingDat()
sphab <- unique(banding.dat.clean[,c('Specie.Name', 'habitat')])

# list results files
files <- list.files("results", pattern="CJS_adult_model_output.rda", full.names = TRUE)

jags2plot <- function(x) {
    load(x)
    species <- str_match(x,pattern="results/(\\w+.\\w+)_CJS")[,2]
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
    species <- str_match(x,pattern="results/(\\w+.\\w+)_CJS")[,2]
    fitnull <- data.frame(fit=modelout[[1]]$JAGSoutput$sims.list$fit, fit.new=modelout[[1]]$JAGSoutput$sims.list$fit.new) %>%
        summarise(p.val = round(mean(fit.new > fit), 2)) %>%
        mutate(model="Null")
    fithab <- data.frame(fit=modelout[[2]]$JAGSoutput$sims.list$fit, fit.new=modelout[[1]]$JAGSoutput$sims.list$fit.new) %>%
        summarise(p.val = round(mean(fit.new > fit), 2)) %>%
        mutate(model="Habitat")
    fittime <- data.frame(fit=modelout[[3]]$JAGSoutput$sims.list$fit, fit.new=modelout[[1]]$JAGSoutput$sims.list$fit.new) %>%
        summarise(p.val = round(mean(fit.new > fit), 2)) %>%
        mutate(model="Time")
    fit <- rbind(fitnull, fithab, fittime) %>%
        mutate(species=species)
}

plot.dat <- lapply(files, jags2plot)
plot.dat <- do.call("rbind", plot.dat)

fit <- lapply(files, getFit)
fit <- do.call("rbind", fit)

out <- ggplot(plot.dat, aes(x=param, y=mean, colour=Model)) + 
    geom_point() + 
    geom_errorbar(aes(ymin=lci, ymax=uci), width=0.1) + 
    facet_wrap(~species, ncol=3) +
    xlab(expression("Parameter")) + ylab(expression(paste("Mean"%+-%"95% Credible Interval"))) + 
    scale_x_discrete(breaks=c('pnull', 'phinull', 'phab', 'phiintro', 'phinative', 'phishrub', 'ptime', 'phitime', 'sigma2'),
                     labels=c(expression('p'['null']), expression(phi['null']), 
                              expression('p'['habitat']), expression(phi['introduced']), expression(phi['native']), expression(phi['shrub']),
                              expression('p'['time']), expression(phi['time'], sigma^2)))
    
save_plot(out, filename = "adult_survival.png", base_height = 12, base_width = 20)

# Results for model with random time effects -----------------------------------
survival.res <- lapply(modelout, function(x) {
    out <- data.frame(x$BUGSoutput$summary[,c('mean', '2.5%', '97.5%')])
    rownames(out) <- gsub("[.]", "_", rownames(out))
    return(out)
})

survival.res <- do.call("rbind", survival.res)
mod_stat <- t(data.frame(strsplit(rownames(survival.res), "[.]")))
survival.res <- cbind(survival.res, mod_stat)
colnames(survival.res) <- c('mean', 'lci', 'uci', "model", "stat")
statlookup <- data.frame(stat = c("mean_p", "mean_phi", "sigma2", "sigma2_real", "beta"),
                         Stat = c("p", "Phi", "Variance (logit)", "Variance (real)", "Beta"))
survival.res <- merge(survival.res, statlookup)



# get the fit stats out of the BUGS objects and into a dataframe with species names
fit.stats <- lapply(modelout, function(x){
    data.frame(fit=x$BUGSoutput$sims.list$fit, 
               fit.new=x$BUGSoutput$sims.list$fit.new)
})

fit.stats <- do.call("rbind", fit.stats)
model <- rownames(fit.stats)
model <- gsub("[[:digit:]]", "", model)
model <- gsub("[.]", "", model)
fit.stats$model <- model

# plot fit against fit.new with the 1:1 line
ggplot(fit.stats, aes(x = fit, y = fit.new)) + 
    geom_point(shape = 1) + 
    geom_abline(intercept = 0, slope = 1) + 
    scale_x_continuous(name = "Discrepancy simulated data") + 
    scale_y_continuous(name = "Discrepancy observed data") + 
    facet_wrap(~model) + 
    theme_classic()

ggsave("results/time_model_fit.png", width = 12, height = 8)

# get the Bayesian p-values
p.vals <- group_by(fit.stats, model) %>%
    summarise(p.val = round(mean(fit.new > fit), 2))

# get and plot the posterior distribution for sigma (to check it's defined)
sigma <- lapply(cjs.mnl.time.ran, function(x) {
    data.frame(sigma = sqrt(x$sims.list$sigma2))
})

sigma <- do.call("rbind", sigma)
sp <- rownames(sigma)
sp <- gsub("[[:digit:]]", "", sp)
sp <- gsub("[.]", "", sp)
sigma$species <- sp

ggplot(sigma, aes(x = sigma)) + 
    geom_density() + 
    facet_wrap(~species) + 
    geom_hline(yintercept = 0.1) + 
    theme_classic()

ggsave("results/sigma_identifiable.png", width = 12, height = 8)

# output plot of phi estimates and 95% CRI between each sampling occasion
phi.est <- lapply(cjs.mnl.time.ran, function(x) {
    data.frame(sampling.occasion = 1:ncol(x$sims.list$phi),
               lower = apply(x$sims.list$phi, 2, function(x) quantile(x, 0.025)),
               upper = apply(x$sims.list$phi, 2, function(x) quantile(x, 0.975)),
               mean.t = colMeans(x$sims.list$phi), 
               mean.phi = x$mean$mean.phi)
    }
    )
    
phi.est <- do.call("rbind", phi.est)
sp <- rownames(phi.est)
sp <- gsub("[[:digit:]]", "", sp)
sp <- gsub("[.]", "", sp)
phi.est$species <- sp

ggplot(phi.est, aes(x=sampling.occasion, y=mean.t)) + 
    facet_wrap(~species) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width = 0, col = "gray45") +
    geom_point(size = 3) + geom_line() + 
    geom_hline(aes(yintercept = mean.phi), col = "red") + 
    labs(x = "Sampling session", y = "Survival probability (with 95% CRI)") + 
    theme_classic()

ggsave("results/time_model_vals.png", width = 12, height = 8)

# Results for model with fixed group effects -----------------------------------
habitat.survival.res <- lapply(cjs.mnl.habitat, function(x) {
    out <- data.frame(t(rbind(data.frame(x$mean), data.frame(x$sd))))
    rownames(out) <- gsub("[.]", "_", rownames(out))
    return(out)
})

habitat.survival.res <- do.call("rbind", habitat.survival.res)
sp_stat <- t(data.frame(strsplit(rownames(habitat.survival.res), "[.]")))
habitat.survival.res <- cbind(habitat.survival.res, sp_stat)
colnames(habitat.survival.res) <- c("Mean", "SD", "Species", "stat")
statlookup <- data.frame(stat = c("mean_p", "mean_phinative", "mean_phiintro", "mean_phiscrub"),
                         Stat = c("p", "Phi (native)", "Phi (introduced)", "Phi (native shrubs)"))
habitat.survival.res <- merge(habitat.survival.res, statlookup)
habitat.survival.res <- arrange(habitat.survival.res, Species, Stat) %>%
    select(Species, Stat, Mean, SD) %>%
    mutate(Mean = round(Mean, 2), 
           SD = round(SD, 2))

# get the fit stats out of the BUGS objects and into a dataframe with species names
fit.stats.hab <- lapply(cjs.mnl.habitat, function(x){
    data.frame(fit=x$sims.list$fit, 
               fit.new=x$sims.list$fit.new)
})

fit.stats.hab <- do.call("rbind", fit.stats.hab)
sp <- rownames(fit.stats.hab)
sp <- gsub("[[:digit:]]", "", sp)
sp <- gsub("[.]", "", sp)
fit.stats.hab$species <- sp

# plot fit against fit.new with the 1:1 line
ggplot(fit.stats.hab, aes(x = fit, y = fit.new)) + 
    geom_point(shape = 1) + 
    geom_abline(intercept = 0, slope = 1) + 
    scale_x_continuous(name = "Discrepancy simulated data") + 
    scale_y_continuous(name = "Discrepancy observed data") + 
    facet_wrap(~species) + 
    theme_classic()

ggsave("results/habitat_model_fit.png", width = 12, height = 8)

# get the Bayesian p-values
p.vals.hab <- group_by(fit.stats.hab, species) %>%
    summarise(p.val.hab = round(mean(fit.new > fit), 2)) %>%
    select(p.val.hab)

# output plot of phi estimates and 95% CRI between each sampling occasion
phi.est.hab <- lapply(cjs.mnl.habitat, function(x) {
    data.frame(habitat = c("Native", "Introduced", "Native shrubs"),
               mean = c(mean(x$sims.list$mean.phinative), 
                        mean(x$sims.list$mean.phiintro),
                        mean(x$sims.list$mean.phiscrub)),
               lower = c(quantile(x$sims.list$mean.phinative, 0.025), 
                         quantile(x$sims.list$mean.phiintro, 0.025),
                         quantile(x$sims.list$mean.phiscrub, 0.025)),
               upper = c(quantile(x$sims.list$mean.phinative, 0.975), 
                         quantile(x$sims.list$mean.phiintro, 0.975),
                         quantile(x$sims.list$mean.phiscrub, 0.975))
               )
}
)

phi.est.hab <- do.call("rbind", phi.est.hab)
sp <- rownames(phi.est.hab)
sp <- gsub("[[:digit:]]", "", sp)
sp <- gsub("[.]", "", sp)
phi.est.hab$species <- sp

ggplot(phi.est.hab, aes(x=habitat, y=mean)) + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.2, col = "gray45") +
    geom_point(size = 3) + 
    facet_wrap(~species) +
    labs(x = "Habitat", y = "Survival probability (with 95% CRI)") + 
    theme_classic()

ggsave("results/habitat_model_vals.png", width = 12, height = 8)

# Results for constant model ---------------------------------------------------
constant.survival.res <- lapply(cjs.mnl.constant, function(x) {
    out <- data.frame(t(rbind(data.frame(x$mean), data.frame(x$sd))))
    rownames(out) <- gsub("[.]", "_", rownames(out))
    return(out)
})

constant.survival.res <- do.call("rbind", constant.survival.res)
sp_stat <- t(data.frame(strsplit(rownames(constant.survival.res), "[.]")))
constant.survival.res <- cbind(constant.survival.res, sp_stat)
colnames(constant.survival.res) <- c("Mean", "SD", "Species", "stat")
statlookup <- data.frame(stat = c("mean_p", "mean_phi"),
                         Stat = c("p", "Phi"))
constant.survival.res <- merge(constant.survival.res, statlookup)
constant.survival.res <- arrange(constant.survival.res, Species, Stat) %>%
    select(Species, Stat, Mean, SD) %>%
    mutate(Mean = round(Mean, 2), 
           SD = round(SD, 2))

# get the fit stats out of the BUGS objects and into a dataframe with species names
fit.stats.constant <- lapply(cjs.mnl.constant, function(x){
    data.frame(fit=x$sims.list$fit, 
               fit.new=x$sims.list$fit.new)
})

fit.stats.constant <- do.call("rbind", fit.stats.constant)
sp <- rownames(fit.stats.constant)
sp <- gsub("[[:digit:]]", "", sp)
sp <- gsub("[.]", "", sp)
fit.stats.constant$species <- sp

# plot fit against fit.new with the 1:1 line
ggplot(fit.stats.constant, aes(x = fit, y = fit.new)) + 
    geom_point(shape = 1) + 
    geom_abline(intercept = 0, slope = 1) + 
    scale_x_continuous(name = "Discrepancy simulated data") + 
    scale_y_continuous(name = "Discrepancy observed data") + 
    facet_wrap(~species) + 
    theme_classic()

ggsave("results/constant_model_fit.png", width = 12, height = 8)

# get the Bayesian p-values
p.vals.constant <- group_by(fit.stats.constant, species) %>%
    summarise(p.val.constant = round(mean(fit.new > fit), 2))  %>%
    select(p.val.constant)

# output plot of phi estimates and 95% CRI between each sampling occasion
phi.est.constant <- lapply(cjs.mnl.constant, function(x) {
    data.frame(mean = x$mean$mean.phi,
               lower = quantile(x$sims.list$mean.phi, 0.025), 
               upper = quantile(x$sims.list$mean.phi, 0.975))
}
)

phi.est.constant <- do.call("rbind", phi.est.constant)
sp <- rownames(phi.est.constant)
sp <- gsub("[[:digit:]]", "", sp)
sp <- gsub("[.]", "", sp)
phi.est.constant$species <- sp

ggplot(phi.est.constant, aes(x=species, y=mean)) + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.2, col = "gray45") +
    geom_point(size = 3) + 
    labs(x = "Species", y = "Survival probability (with 95% CRI)") + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("results/constant_model_vals.png", width = 12, height = 8)

# Results for TSM model --------------------------------------------------------
tsm.survival.res <- lapply(cjs.mnl.tsm, function(x) {
    out <- data.frame(t(rbind(data.frame(x$mean), data.frame(x$sd))))
    rownames(out) <- gsub("[.]", "_", rownames(out))
    return(out)
})

tsm.survival.res <- do.call("rbind", tsm.survival.res)
sp_stat <- t(data.frame(strsplit(rownames(tsm.survival.res), "[.]")))
tsm.survival.res <- cbind(tsm.survival.res, sp_stat)
colnames(tsm.survival.res) <- c("Mean", "SD", "Species", "stat")
statlookup <- data.frame(stat = c("mean_p", "mean_phitsm1", "mean_phitsm2"),
                         Stat = c("p", "Phi (1)", "Phi (2+)"))
tsm.survival.res <- merge(tsm.survival.res, statlookup)
tsm.survival.res <- arrange(tsm.survival.res, Species, Stat) %>%
    select(Species, Stat, Mean, SD) %>%
    mutate(Mean = round(Mean, 2), 
           SD = round(SD, 2))

# get the fit stats out of the BUGS objects and into a dataframe with species names
fit.stats.tsm <- lapply(cjs.mnl.tsm, function(x){
    data.frame(fit=x$sims.list$fit, 
               fit.new=x$sims.list$fit.new)
})

fit.stats.tsm <- do.call("rbind", fit.stats.tsm)
sp <- rownames(fit.stats.tsm)
sp <- gsub("[[:digit:]]", "", sp)
sp <- gsub("[.]", "", sp)
fit.stats.tsm$species <- sp

# plot fit against fit.new with the 1:1 line
ggplot(fit.stats.tsm, aes(x = fit, y = fit.new)) + 
    geom_point(shape = 1) + 
    geom_abline(intercept = 0, slope = 1) + 
    scale_x_continuous(name = "Discrepancy simulated data") + 
    scale_y_continuous(name = "Discrepancy observed data") + 
    facet_wrap(~species) + 
    theme_classic()

ggsave("results/tsm_model_fit.png", width = 12, height = 8)

# get the Bayesian p-values
p.vals.tsm <- group_by(fit.stats.tsm, species) %>%
    summarise(p.val.tsm = round(mean(fit.new > fit), 2)) %>%
    select(p.val.tsm)

# output plot of phi estimates and 95% CRI between each sampling occasion
phi.est.tsm <- lapply(cjs.mnl.tsm, function(x) {
    data.frame(cap.period = c("Phi.1", "Phi.2"),
               mean = c(x$mean$mean.phitsm1, x$mean$mean.phitsm2),
               lower = c(quantile(x$sims.list$mean.phitsm1, 0.025),
                         quantile(x$sims.list$mean.phitsm2, 0.025)),
               upper = c(quantile(x$sims.list$mean.phitsm2, 0.975)))
}
)

phi.est.tsm <- do.call("rbind", phi.est.tsm)
sp <- rownames(phi.est.tsm)
sp <- gsub("[[:digit:]]", "", sp)
sp <- gsub("[.]", "", sp)
phi.est.tsm$species <- sp

ggplot(phi.est.tsm, aes(x=cap.period, y=mean)) + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.2, col = "gray45") +
    geom_point(size = 3) + 
    facet_wrap(~species, nrow=4) +
    labs(x = "Species", y = "Survival probability (with 95% CRI)") + 
    theme_classic() 

ggsave("results/tsm_model_vals.png", width = 12, height = 8)

p.vals <- cbind(p.vals.time, p.vals.hab, p.vals.constant, p.vals.tsm)
species.family <- tax_name(query = p.vals$species,get = "family", db = "itis")
p.vals <- cbind(p.vals, species.family)
p.vals$family[2:3] <- "Parulidae" # these species were classified differently
p.vals <- arrange(p.vals, family, species) %>%
    select(1, 6, 4, 5, 3, 2) 
test <- abs(p.vals[,3:6] - 0.5)
p.vals$model <- apply(test, 1, which.min)
p.vals$model <- factor(p.vals$model, levels = c(1, 2, 3, 4), labels = c("Constant", "TSM", "Habitat", "Time"))

best.mod.count <- group_by(p.vals, model) %>%
    summarise(count = length(model))

write.csv(habitat.survival.res, "results/habitat_mod.csv")
write.csv(survival.res, "results/time_mod.csv")
write.csv(p.vals, "results/model_comparison.csv")
write.csv(tsm.survival.res, "results/tsm_mod.csv")
write.csv(constant.survival.res, "results/constant_mod.csv")
