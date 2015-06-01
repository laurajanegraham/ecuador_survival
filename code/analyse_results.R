# load packages
require(ggplot2)
require(R2WinBUGS)
require(dplyr)
require(taxize)

# load results for each model
load("results/cjs.mnl.habitat.rda")
load("results/cjs.mnl.time.ran.rda")
load("results/cjs.mnl.tsm.rda")
load("results/cjs.mnl.constant.rda")

# Results for model with randon time effects -----------------------------------
survival.res <- lapply(cjs.mnl.time.ran, function(x) {
    out <- data.frame(t(rbind(data.frame(x$mean[-1]), data.frame(x$sd[-1]))))
    rownames(out) <- gsub("[.]", "_", rownames(out))
    return(out)
})

survival.res <- do.call("rbind", survival.res)
sp_stat <- t(data.frame(strsplit(rownames(survival.res), "[.]")))
survival.res <- cbind(survival.res, sp_stat)
colnames(survival.res) <- c("Mean", "SD", "Species", "stat")
statlookup <- data.frame(stat = c("mean_p", "mean_phi", "sigma2", "sigma2_real"),
                         Stat = c("p", "Phi", "Variance (logit)", "Variance (real)"))
survival.res <- merge(survival.res, statlookup)
survival.res <- arrange(survival.res, Species, Stat) %>%
    select(Species, Stat, Mean, SD) %>%
    mutate(Mean = round(Mean, 2), 
           SD = round(SD, 2))


# get the fit stats out of the BUGS objects and into a dataframe with species names
fit.stats <- lapply(cjs.mnl.time.ran, function(x){
    data.frame(fit=x$sims.list$fit, 
               fit.new=x$sims.list$fit.new)
})

fit.stats <- do.call("rbind", fit.stats)
sp <- rownames(fit.stats)
sp <- gsub("[[:digit:]]", "", sp)
sp <- gsub("[.]", "", sp)
fit.stats$species <- sp

# plot fit against fit.new with the 1:1 line
ggplot(fit.stats, aes(x = fit, y = fit.new)) + 
    geom_point(shape = 1) + 
    geom_abline(intercept = 0, slope = 1) + 
    scale_x_continuous(name = "Discrepancy simulated data") + 
    scale_y_continuous(name = "Discrepancy observed data") + 
    facet_wrap(~species) + 
    theme_classic()

ggsave("results/time_model_fit.png", width = 12, height = 8)

# get the Bayesian p-values
p.vals.time <- group_by(fit.stats, species) %>%
    summarise(p.val.time = round(mean(fit.new > fit), 2))

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
               mean = colMeans(x$sims.list$phi))
    }
    )
    
phi.est <- do.call("rbind", phi.est)
sp <- rownames(phi.est)
sp <- gsub("[[:digit:]]", "", sp)
sp <- gsub("[.]", "", sp)
phi.est$species <- sp

ggplot(phi.est, aes(x=sampling.occasion, y=mean)) + 
    facet_wrap(~species) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width = 0, col = "gray45") +
    geom_point(size = 3) + geom_line() + 
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
