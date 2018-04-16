library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

res <- read.csv("results/species/nojuv_raw_results.csv")
res <- filter(res, species %in% c("ERIOCNEMIS LUCIANI", "METALLURA TYRIANTHINA", "DIGLOSSA HUMERALIS"), model %in% c("null", "habitat"))
res <- droplevels(res)
res$Probabilities <- c("Recapture (null)", "Survival (all habitats)", "Recapture (habitat)", "Survival (introduced)", "Survival (native)", "Survival (scrub)", "Recapture (null)", "Survival (all habitats)", "Recapture (habitat)", "Survival (introduced)", "Survival (native)", "Survival (scrub)", "Recapture (null)", "Survival (all habitats)", "Recapture (habitat)", "Survival (introduced)", "Survival (native)", "Survival (scrub)")
res$Probabilities <- factor(res$Probabilities, levels=c("Recapture (null)", "Survival (all habitats)", "Recapture (habitat)", "Survival (introduced)", "Survival (native)", "Survival (scrub)"))
res_el <- filter(res, species == "ERIOCNEMIS LUCIANI")
res_mt <- filter(res, species == "METALLURA TYRIANTHINA")
res_dh <- filter(res, species == "DIGLOSSA HUMERALIS")

el <- ggplot(res_el, aes(x = Probabilities, y = mean)) + 
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=.1) +
  geom_point() + 
  ylab(expression(paste("Mean " %+-% " 95% CI"))) + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

save_plot("../../SCALEFORES/CMR Spatial Modelling Talk/eluciani_null.png", el, base_width=6)

mt <- ggplot(res_mt, aes(x = Probabilities, y = mean)) + 
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=.1) +
  geom_point() + 
  ylab(expression(paste("Mean " %+-% " 95% CI"))) + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

save_plot("../../SCALEFORES/CMR Spatial Modelling Talk/mtyrianthina_null.png", mt, base_width=6)

dh <- ggplot(res_dh, aes(x = Probabilities, y = mean)) + 
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=.1) +
  geom_point() + 
  ylab(expression(paste("Mean " %+-% " 95% CI"))) + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

save_plot("../../SCALEFORES/CMR Spatial Modelling Talk/dhumeralis_null.png", dh, base_width=6)
