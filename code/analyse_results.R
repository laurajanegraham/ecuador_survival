# load packages
require(ggplot2)
require(R2WinBUGS)
require(dplyr)

# load results for each model
load("results/cjs.mnl.habitat.rda")
load("results/cjs.mnl.time.ran.rda")

# overall survival probabilities for the random time effects
# with sd
survival.res <- sapply(cjs.mnl.time.ran, function(x) {
        c(mean.phi = round(x$mean$mean.phi, 2), 
          sd.phi = round(x$sd$mean.phi, 2), 
          mean.p = round(x$mean$mean.p, 2), 
          sd.p = round(x$sd$mean.p, 2), 
          mean.sigma2 = round(x$mean$sigma2.real, 2), 
          sd.sigma2 = round(x$sd$sigma2.real, 2))
})

survival.res <- t(survival.res)

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
    facet_wrap(~species)

p.vals <- group_by(fit.stats, species) %>%
    summarise(p.val = mean(fit.new > fit))
mean(bugs.res$sims.list$fit.new > bugs.res$sims.list$fit)

# for same species look at the posterior distributions for each individual year
phi.post <- bugs.res$sims.list$phi
phi.post <- data.frame(phi.post)
phi.post$id <- rownames(phi.post)
phi.post <- gather(phi.post, "time", "posterior", -id)

require(ggplot2)
ggplot(phi.post, aes(x = posterior)) + geom_density() + facet_wrap(~time) + geom_hline(yintercept = 1)

sigma <- sqrt(bugs.res$sims.list$sigma2)
hist(sigma, freq = FALSE, main = "Prior = Uniform(0, 10)", xlab=expression(sigma))
abline(h = 0.1, lty = 2, lwd = 2)
