load("results/cjs.group.habitat.rda")
load("results/cjs.mnl.time.ran.rda")


cjs.group.habitat[[1]]$sims.list$phi.g
cjs.group.habitat[[1]]$mean$phi.g
group.dic <- data.frame(sapply(cjs.group.habitat, function(x) round(x$DIC, 2)))
group.dic$Species <- rownames(group.dic)
time.dic <- data.frame(sapply(cjs.mnl.time.ran, function(x) round(x$DIC, 2)))
time.dic$Species <- rownames(time.dic)
dic.values <- merge(time.dic, group.dic)
names(dic.values) <- c("Species", "time", "habitat")
save(dic.values, file = "results/dicvalues.rda")

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

# means only
survival.res <- sapply(cjs.mnl.time.ran, function(x) {
        c(mean.phi = round(x$mean$mean.phi, 2), 
          mean.p = round(x$mean$mean.p, 2), 
          mean.sigma2 = round(x$mean$sigma2.real, 2)) 
})

survival.res <- t(survival.res)
save(survival.res, file="results/survival.res")

# for one species (DIGLOSSOPIS CYANEA - best fitting model by DIC) do some model fitting
species <- "METALLURA TYRIANTHINA"   

bugs.res <- cjs.mnl.time.ran[[species]]

plot(bugs.res$sims.list$fit, bugs.res$sims.list$fit.new, xlab = "Discrepancy actual data", 
     ylab = "Discrepancy simulated data", las = 1)
abline(0, 1, lwd = 2)

mean(bugs.res$sims.list$fit.new > bugs.res$sims.list$fit)

# for same species look at the posterior distributions for each individual year
phi.post <- bugs.res$sims.list$phi
phi.post <- data.frame(phi.post)
phi.post$id <- rownames(phi.post)
phi.post <- gather(phi.post, "time", "posterior", -id)

require(ggplot2)
ggplot(phi.post, aes(x = posterior)) + geom_density() + facet_wrap(~time) + geom_hline(yintercept = 1)

