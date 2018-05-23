# Load packages and functions --------------------------------------------------
require(RMark)
source("code/fnCleanBandingDat.R")
source("code/fnEncounterHistory.R")
source("code/fnSurvivalAnalysisMARK.R")

# Load and format banding data -------------------------------------------------
banding.dat.clean <- CleanBandingDat()

# Analysis for individual species ----------------------------------------------
species.list <- unique(banding.dat.clean$Specie.Name) 

# initialise list for storing results
res <- list()

for(species in species.list){
        # get records for species, get necessary columns & set the habitat type
        sp_dat <- filter(banding.dat.clean, Specie.Name == species) %>%
                select(Band.Number, session_new, habitat)
                
        # create species' encounter history
        sp_eh <- EncounterHistory(sp_dat[,c("Band.Number", "session_new")], 
                                  "session_new", 
                                  "Band.Number")
        
        # Process data for MARK ------------------------------------------------
        # get encounter history formatted correctly for MARK
        sp_mark <- sp_eh$eh.mark
        
        # get the habitat data onto the MARK formatted data
        sp_habitat <- select(sp_dat, Band.Number, habitat) %>%
                unique() %>%
                merge(sp_mark) %>%
                group_by(Band.Number, ch) %>%
                summarise(habitat = first(habitat))

        sp_habitat$habitat <- factor(sp_habitat$habitat)

        # process data for use in MARK
        sp_proc <- process.data(data.frame(sp_habitat), 
                                model = "CJS", 
                                group = "habitat")
        
        sp_ddl <- make.design.data(sp_proc)
        
        # create age field for TSM models
        sp_ddl=add.design.data(sp_proc,sp_ddl,parameter="Phi",type="age",
                            bins=c(0,.5,21),name="TSM")
        
        # run MARK analysis ----------------------------------------------------
        res[[species]] <- SurvivalAnalysisMARK(sp_proc, sp_ddl)
}

