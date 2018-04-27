#load packages
library(plyr)
library(tidyverse)
library(data.table)
library(R2jags)


# Load and format banding data ----
data_input <- c("data/banding_sheet.csv", "data/2016_data.csv")

species_info <- read_csv("data/evaluated_species.csv")

dat <- ldply(data_input, fread) %>% 
    # remove juveniles and records with no species or band number
    filter(Age != "Y", 
         !(`Band Number` %in% c("99999", "9999", "", "SIN ANILLO")), 
         !(`Specie Name` %in% c("", "U", "99999", "9"))) %>% 
  # add a habitat and year column
  mutate(habitat = ifelse(Location %in% c("LLAV", "SANA"), "Scrub", ifelse(Location == "MASE", "Native", "Introduced")),
         year = str_sub(Session, start = -2) %>% as.factor %>% as.numeric) %>% 
  # select required columns
  select(id = "Band Number", species = "Specie Name", habitat, year) %>% 
  # join to species info to get the specialisms and remove species we are not analysing
  inner_join(species_info) %>% 
  # get unique yearly records (combining across the three sampling sessions)
  group_by(year, id, species, specialization) %>% 
  summarise(habitat = names(which.max(table(habitat)))) %>% 
  mutate(recap = 1) %>% ungroup()

# create sampsize vector (list of nind for each species) ----
sampsize <- select(dat, species, id) %>% unique %>% group_by(species) %>% summarise(count = n()) %>% pull(count)

# create the y matrix (nind * nocc * nspecies recapture history) ----
y_2d <- select(dat, species, id, year, recap) %>% 
  spread(year, recap, fill = 0)

y_list <- split(y_2d, f = y_2d$species) %>% 
  lapply(., function(x) {
    df = as.matrix(select(x, -species, -id))
    na_rows = matrix(NA, nrow = max(sampsize) - nrow(df), ncol = ncol(df))
    out = rbind(df, na_rows)
  }
  )

y <- simplify2array(y_list)

# create the TSM matrix (nind * nocc * nspecies banding occasion (1 - first, 2 - all others)) ----
tsm_2d <- select(dat, species, id, year) %>% 
  mutate(year = str_pad(year, 2, pad = "0")) %>% 
  # first group-arrange-mutate gets labels the first (1) and all following (2) banding occasions
  group_by(species, id) %>% 
  arrange(species, id, year) %>% 
  mutate(recap = row_number()) %>% 
  mutate(recap = ifelse(recap == 1, 1, 2)) %>% 
  # spread and gather gets the years with no recaps included as 0
  spread(year, recap, fill = 0) %>% 
  gather(key = year, value = recap, -species, -id) %>% 
  # second group-arrange-mutate fills in the blanks
  group_by(species, id) %>% 
  arrange(species, id, year) %>% 
  mutate(recap = ifelse(recap == 0, NA, recap)) %>% 
  fill(recap) %>% 
  spread(year, recap) %>% ungroup

tsm_list <- split(tsm_2d, f = tsm_2d$species) %>% 
  lapply(., function(x) {
    df = as.matrix(select(x, -species, -id))
    na_rows = matrix(NA, nrow = max(sampsize) - nrow(df), ncol = ncol(df))
    out = rbind(df, na_rows)
  }
  )

tsm <- simplify2array(tsm_list)

# create the habitat matrix (nind * nocc * nspecies habitat where recaptured) ----
hab_2d <- mutate(dat, fhabitat = habitat %>% as.factor %>% as.numeric) %>% 
  select(species, id, year, fhabitat) %>% 
  mutate(year = str_pad(year, 2, pad = "0")) %>% 
  # spread and gather gets the years with no recaps included as 0
  spread(year, fhabitat, fill = 0) %>% 
  gather(key = year, value = fhabitat, -species, -id) %>% 
  # second group-arrange-mutate fills in the blanks
  group_by(species, id) %>% 
  arrange(species, id, year) %>% 
  mutate(fhabitat = ifelse(fhabitat == 0, NA, fhabitat)) %>% 
  fill(fhabitat) %>% 
  spread(year, fhabitat) %>% ungroup
  

hab_list <- split(hab_2d, f = hab_2d$species) %>% 
  lapply(., function(x) {
    df = as.matrix(select(x, -species, -id))
    na_rows = matrix(NA, nrow = max(sampsize) - nrow(df), ncol = ncol(df))
    out = rbind(df, na_rows)
  }
  )

hab <- simplify2array(hab_list)

# create f matrix (nind * nspecies occasion when banded) ----
f <- group_by(dat, species, id) %>% summarise(year = min(year)) %>% 
  group_by(species) %>% mutate(id = row_number()) %>% 
  spread(species, year) %>% select(-id) %>% as.matrix()

# create vector of specialism group (list of specialism for each species) ----
g <- arrange(species_info, species) %>% pull(specialization)

# bundle data for jags ----
null_data <- list(y = y, f = f, nind = sampsize, n.occ = dim(y)[2], nspec = dim(y)[3])
tsm_data <- list(y = y, f = f, tsm = tsm, nind = sampsize, n.occ = dim(y)[2], nspec = dim(y)[3])
hab_data <- list(y = y, f = f, tsm = tsm, hab = hab, nind = sampsize, n.occ = dim(y)[2], nspec = dim(y)[3])
full_data <- list(y = y, f = f, tsm = tsm, hab = hab, g = g, nind = sampsize, n.occ = dim(y)[2], nspec = dim(y)[3])

# set initial values ---- 
# initial values for z will be set as 1 for each occasion following the occasion
# on which an individual was first banded
zinit <- array(dim = dim(y))
for(s in 1:dim(y)[3]){
  for (i in 1:sampsize[s]) {
    for (t in (f[i,s] + 1):dim(y)[2]) {
      if(f[i,s]+1 <= dim(y)[2]) zinit[i, t, s] <- 1
    }
  }
}

inits <- function(){list(z = zinit)}

# run model in jags ----
# estimate per species, no further delineation of parameters estimated
null_out <- jags(data = null_data, model.file = "multi_sp_survival.bug.R", inits = inits, n.chains = 3, parameters.to.save = c("speciesphi", "speciesp", "mean.phi", "mean.p", "sd.lphi", "sd.lp"))

# estimate per species for survival between first and second banding, and all other
# bandings (TSM model, accounts for transience)
tsm_out <- jags(data = tsm_data, model.file = "multi_sp_survival_tsm.bug.R", inits = inits, n.chains = 3, parameters.to.save = c("speciesphi", "speciesp", "mean.phi", "mean.p", "sd.lphi", "sd.lp"))

# as tsm but with an estimate per habitat too - allows us to look at the effect of habitat
# on survival
hab_out <- jags(data = hab_data, model.file = "multi_sp_survival_hab.bug.R", inits = inits, n.chains = 3, parameters.to.save = c("speciesphi", "speciesp", "mean.phi", "mean.p", "sd.lphi", "sd.lp", "indiff", "isdiff", "nsdiff"))

# as hab, but with p estimated per habitat also - allows us to look at (and control for?)
# the difference in detectability between habitats
habp_out <- jags(data = hab_data, model.file = "multi_sp_survival_habp.bug.R", inits = inits, n.chains = 3, parameters.to.save = c("phi", "p", "mean.phi", "mean.p", "sd.lphi", "sd.lp", "indiff", "isdiff", "nsdiff"))

# as habp, but the specialisms are included which means we can look at how specialism
# survival varies across habitats
full_out <- jags(data = full_data, model.file = "multi_sp_survival_full.bug.R", inits = inits, n.chains = 3, parameters.to.save = c("speciesphi", "speciesp", "mean.phi", "mean.p", "sd.lphi", "sd.lp", "indiff", "isdiff", "nsdiff"))

# Get output into useful tables (summary and posterior simulations) ----
# posterior simulation array
full_sims <- habp_out$BUGSoutput$sims.array
# will need the names because lost in next step
params <- names(habp_out[1,1,])
# convert to 2d array
dim(full_sims) <- c(prod(dim(full_sims)[1:2]), dim(full_sims)[3])
# re-add the parameter names
colnames(full_sims) <- params

full_sims_narrow <- as.tibble(full_sims) %>% 
  gather(key = param, value = value) %>% 
  mutate(param = str_replace(param, "\\[", "_"),
         param = str_replace(param, "\\]", "")) %>% 
  separate(param, into = c("param_name", "param_index"), sep = "_") 

# posterior summaries
full_summary <- as_tibble(full_out$BUGSoutput$summary, rownames = "param") %>% 
  mutate(param = str_replace(param, "\\[", "_"),
         param = str_replace(param, "\\]", "")) %>% 
  separate(param, into = c("param_name", "param_index"), sep = "_")

# Survival of specialisms across habitat types ----
# get the posterior simulations for mean.phi for the second banding occasion onwards in
# order to plot comparison between habitat survival distributions for each specialism
hab_comp <- filter(full_sims_narrow, param_name == "mean.phi") %>% 
  separate(param_index, into = c("a", "h", "g")) %>% 
  filter(a == 2) %>% 
  mutate(Habitat = factor(h, levels=c(2, 1, 3), labels=c("Native", "Introduced", "Scrub")),
         Specialism = factor(g, labels=c("Open habitat generalist", "Scrubby generalist", 
                                         "Forest Specialist"))) %>% 
  select(Habitat, Specialism, value)

# plot the results
theme_set(theme_classic() + theme(strip.background = element_blank()))

ggplot(hab_comp, aes(x=value, fill=Habitat)) + 
  geom_density(alpha=0.5, position="identity") + 
  scale_fill_viridis_d() + 
  facet_grid(. ~ Specialism) + 
  xlab(expression("Mean "*phi))

# Variability of detection across habitat types ----
habp_comp <- filter(full_sims_narrow, param_name == "mean.p") %>% 
  mutate(Habitat = factor(param_index, levels=c(2, 1, 3), labels=c("Native", "Introduced", "Scrub"))) %>% 
  select(Habitat, value)

ggplot(habp_comp, aes(x=value, fill=Habitat)) + 
  geom_density(alpha=0.5, position="identity") + 
  scale_fill_viridis_d() + 
  xlab("Mean p")

# Individual species results (mean)
survival_comp <- filter(full_summary, param_name == "speciesphi") %>% 
  bind_cols(arrange(species_info, species)) %>% 
  select(species, specialization, mean_val = mean, LCI = `2.5%`, UCI = `97.5%`)

ggplot(survival_comp, aes(x = species, y = mean_val)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI)) + 
  scale_y_continuous(limits = c(0, 1))
  
