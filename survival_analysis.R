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

# TSM data
tsm_2d <- select(dat, species, id, year) %>% 
  group_by(species, id) %>% 
  arrange(species, id, year) %>% 
  mutate(recap = row_number()) %>% 
  mutate(recap = ifelse(recap == 1, 1, 2)) %>% 
  spread(year, recap, fill = 0)


# create the X matrix (nind * nocc * nspecies habitat where recaptured) ----
x_2d <- mutate(dat, fhabitat = habitat %>% as.factor %>% as.numeric) %>% 
  select(species, id, year, fhabitat) %>% 
  spread(year, fhabitat, fill = NA)

x_list <- split(x_2d, f = x_2d$species) %>% 
  lapply(., function(x) {
    df = as.matrix(select(x, -species, -id))
    na_rows = matrix(NA, nrow = max(sampsize) - nrow(df), ncol = ncol(df))
    out = rbind(df, na_rows)
  }
  )

X <- simplify2array(x_list)

# create f matrix (nind * nspecies occasion when banded) ----
f <- group_by(dat, species, id) %>% summarise(year = min(year)) %>% 
  group_by(species) %>% mutate(id = row_number()) %>% 
  spread(species, year) %>% select(-id) %>% as.matrix()


# bundle data for jags ----
bugs.data <- list(y = y, f = f, nind = sampsize, n.occ = dim(y)[2], nspec = dim(y)[3])

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

inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), 
                         sd.lphi = runif(1, 0, 3), sd.lp = runif(1, 0, 3), z = zinit)}

# run model in jags ----
out <- jags(data = bugs.data, model.file = "multi_sp_survival.bug.R", inits = inits, n.chains = 3, parameters.to.save = c("phi", "p", "mean.phi", "mean.p", "sd.lphi", "sd.lp")

