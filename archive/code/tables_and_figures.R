library(tidyverse)
library(stringr)
library(cowplot)
# get species names in the correct case for display and matching
sp_case_convert <- function(x) {
  x <- tolower(x)
  str_sub(x, 1, 1) <- toupper(str_sub(x, 1, 1))
  return(x)
}

# We have removed some of the species from the displayed results 
# inner join this file to each of the results tables to filter
sp_eval <- read_csv("data/evaluated_species.csv")

# Table 1 - null model results for every species analysed
dat <- read_csv(paste0("results/species/nojuv_raw_results.csv")) %>%
  filter(model %in% c("null", "null.tsm")) %>%
  inner_join(sp_eval) %>% 
  mutate(val=paste0(sprintf("%.2f", round(mean,2)), " (", sprintf("%.2f", round(X2.5.,2)), "-", sprintf("%.2f", round(X97.5.,2)), ")"),
         param=gsub("mean.", "", param)) %>%
  select(tax_identity, param, val, model) %>%
  spread(param, val) %>%
  inner_join(read_csv("results/species/model_comparison.csv")) %>%
  select(-habitat, -habitat.tsm, -bestmod) %>%
  arrange(tax_identity)

dat$bestmod <- apply(dat[,c("null", "null.tsm")], 1, function(x) {
  x <- 0.5-x
  names(which.min(abs(x)))
})

table1 <- filter(dat, model == bestmod) %>% rename(Species = tax_identity) %>%
  mutate(Species = sp_case_convert(Species)) %>%
  left_join(read_csv("data/forest_specialization.csv")) %>% 
  select(Species, specialization, p, phi1, phi2) %>%
  arrange(specialization, Species) %>% write_csv("doc/tables/table1_nullmod_allspecies.csv")

# Table 2 - model comparison table
table2 <- read_csv("results/species/model_comparison.csv") %>%
  inner_join(sp_eval) %>% 
  rename(Species = tax_identity) %>%
  mutate(Species = sp_case_convert(Species)) %>%
  left_join(read_csv("data/forest_specialization.csv")) %>% 
  separate(bestmod, c("Model", "TSM"), fill="right") %>%
  mutate(TSM = ifelse(is.na(TSM), "N", "Y")) %>%
  select(Species, specialization, null, null.tsm, habitat, habitat.tsm, Model, TSM) %>%
  arrange(specialization, Species) %>% write_csv("doc/tables/table2_model_comparison.csv")

# Figure 1 - habitat results
figure1_tab <- read_csv("results/species/nojuv_raw_results.csv") %>%
  inner_join(read_csv("results/species/model_comparison.csv")) %>%
  inner_join(sp_eval) %>% 
  filter(model == bestmod & model %in% c("habitat", "habitat.tsm")) %>%
  rename(Species = tax_identity) %>%
  mutate(Species = sp_case_convert(Species)) %>%
  left_join(read_csv("data/forest_specialization.csv"), all.x = TRUE) %>%
  filter(param == "phi2") %>%
  select(Species, specialization, newhab, mean, X2.5., X97.5.) %>%
  rename(Habitat = newhab)

# get the species names in the right order
sp_names <- distinct(figure1_tab, Species, specialization) %>%
  arrange(specialization, Species)

sp_names$specialization <- factor(sp_names$specialization, labels=c("Open habitat generalists", "Shrubby generalists", "Forest specialists"))
figure1_tab$Species <- factor(figure1_tab$Species, levels=sp_names$Species)

# need to know where lines between specializations are
sp_last <- group_by(sp_names, specialization) %>%
  filter(row_number()==n())

sp_last_position <- which(sp_names$Species %in% sp_last$Species)

sp_first <- group_by(sp_names, specialization) %>%
  filter(row_number()==1)

sp_first_position <- which(sp_names$Species %in% sp_first$Species)

sp_text_position <- (sp_last_position - sp_first_position)/2 + sp_first_position

specializations <- levels(sp_names$specialization)


pd <- position_dodge(1)
figure1 <- ggplot(figure1_tab, aes(x = Species, y = mean, group=factor(Habitat))) + 
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), colour="black", width=.1, position=pd) +
  geom_point(position=pd, aes(shape=factor(Habitat)), size=3) + 
  geom_vline(xintercept = 1:nrow(sp_names) + 0.5, linetype = "dotted") +
  geom_vline(xintercept = sp_last_position[1:length(sp_last_position)-1] + 0.5, linetype = "F1") + 
  annotate("text", x = sp_text_position, y = 1, label = specializations) + 
  labs(x = "Species", y = "Apparent Survival", shape = "Habitat") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

save_plot("doc/figures/figure1_habitatresults.png", figure1, base_width = 19.92, base_height = 7.34)

# supp table 1 - recapture history
# This is created in 'analyse_results.R": "results/species/nojuv_recapture_summaries.csv"
tableS1 <- read_csv("results/species/nojuv_recapture_summaries.csv") %>% 
  mutate(Species = sp_case_convert(tax_level)) %>%
  left_join(read_csv("data/forest_specialization.csv")) %>%
  select(Species, specialization, `0`, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `12`, nind, recaps) %>%
  arrange(specialization, Species) %>% write_csv("doc/tables/tableS1_recapture_history.csv")
  
  
# supp table 2 - all results for the habitat model (incl. the ones which are not well fitting)
dat <- read_csv("results/species/nojuv_raw_results.csv") %>%
  filter(model %in% c("habitat", "habitat.tsm")) %>%
  inner_join(sp_eval) %>% 
  mutate(val=paste0(sprintf("%.2f", round(mean,2)), " (", sprintf("%.2f", round(X2.5.,2)), "-", sprintf("%.2f", round(X97.5.,2)), ")"),
         param=gsub("mean.", "", param)) %>%
  select(tax_identity, newhab, param, val, model) %>%
  spread(param, val) %>%
  inner_join(read_csv("results/species/model_comparison.csv")) %>%
  select(-null, -null.tsm, -bestmod) %>%
  arrange(tax_identity)

dat$bestmod <- apply(dat[,c("habitat", "habitat.tsm")], 1, function(x) {
  x <- 0.5-x
  names(which.min(abs(x)))
})

tableS2 <- filter(dat, model == bestmod) %>% rename(Species = tax_identity) %>%
  mutate(Species = sp_case_convert(Species)) %>%
  left_join(read_csv("data/forest_specialization.csv"), all.x = TRUE) %>% 
  rename(Habitat = newhab) %>%
  select(Species, Habitat, specialization, p, phi1, phi2) %>%
  arrange(specialization, Species) %>%write_csv("doc/tables/tableS2_habitatmod_allspecies.csv")
