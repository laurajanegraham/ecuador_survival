# load packages
require(ggplot2)
require(R2WinBUGS)
require(dplyr)

# load results for each model
load("results/cjs.mnl.habitat.rda")
load("results/cjs.mnl.time.ran.rda")

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
p.vals <- group_by(fit.stats, species) %>%
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
                         Stat = c("p", "Phi (native)", "Phi (introduced)", "Phi (scrub)"))
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
    summarise(p.val.hab = round(mean(fit.new > fit), 2))

p.vals <- merge(p.vals, p.vals.hab)

p.vals$best.mod <- ifelse(abs(0.5 - p.vals$p.val) > abs(0.5 - p.vals$p.val.hab), "habitat", "time") 

write.csv(habitat.survival.res, "results/habitat_mod.csv")
write.csv(survival.res, "results/time_mod.csv")
write.csv(p.vals, "results/model_comparison.csv")
