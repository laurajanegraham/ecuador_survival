library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# function to extract the key parameters from the Bayesian models. 
extract_par <- function(mod) {
    res <- apply(mod$JAGSoutput$sims.array, 2, function(x) {
        mutate(data.frame(x), draw=1:nrow(x)) %>%
        gather(par, estimate, -draw) %>%
            filter(par %in% c("mean.p", "mean.phi", "mean.phi1", "mean.phi2", "mean.phiintro", "mean.phinative", "mean.phiscrub", "sigma2.real"))
    }
    )
    names(res) <- 1:length(res)
    for(x in 1:length(res)) {
        res[[x]]$chain <- x
    }
    res <- do.call("rbind", res)
    return(res)
}

files <- list.files("results", pattern="CJS_full", full.names = TRUE)

pars <- lapply(files, function(x) {
    load(x)
    pars.x <- lapply(modelout, extract_par)
    pars.x <- do.call("rbind", pars.x)
    pars.x$model <- unlist(lapply(strsplit(rownames(pars.x), split="[.]"), function(x) x[1]))
    return(pars.x)
})

species <- gsub("_", " ", str_match(x,pattern="results/(\\w+.\\w+)_CJS")[,2])
names(pars) <- species
pars.all <- do.call("rbind", pars)
pars.all$species <- unlist(lapply(strsplit(rownames(pars.all), split="[.]"), function(x) x[1]))

dat.null <- filter(pars.all, model=="null")
ggplot(dat.null,aes(x=draw,y=estimate,col=as.factor(chain))) + 
    geom_line() + 
    facet_grid(species~par, scales = "free") + 
    theme_classic() + 
    labs(col="Chain") + 
    ggtitle("Null Model")

ggsave("results/convergence_null_mod.pdf", height=15)

dat.hab <- filter(pars.all, model=="habitat")
ggplot(dat.hab,aes(x=draw,y=estimate,col=as.factor(chain))) + 
    geom_line() + 
    facet_grid(species~par, scales = "free") + 
    theme_classic() + 
    labs(col="Chain") + 
    ggtitle("Habitat model")

ggsave("results/convergence_habitat_mod.pdf", height=15)

dat.time <- filter(pars.all, model=="time")
ggplot(dat.time,aes(x=draw,y=estimate,col=as.factor(chain))) + 
    geom_line() + 
    facet_grid(species~par, scales = "free") + 
    theme_classic() + 
    labs(col="Chain") + 
    ggtitle("Null Model")

ggsave("results/convergence_time_mod.pdf", height=15)



