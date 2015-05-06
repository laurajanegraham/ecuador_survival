# Load packages and functions --------------------------------------------------
require(Hmisc)
require(dplyr)
require(tidyr)
require(RMark)
require(R2WinBUGS)
source("code/fnEncounterHistory.R")

# Load and format banding data -------------------------------------------------

# read in the table from full_data.mdb and write out to banding_sheet.csv - this
# gets rid of column attribute issues
banding.dat <- mdb.get("data/full_data.mdb", tables="Banding Sheet")
write.csv(banding.dat, file="data/banding_sheet.csv")

# re-read in data
banding.dat <- read.csv("data/banding_sheet.csv", stringsAsFactors = FALSE)

# remove records with no species, session or band
banding.dat.clean <- filter(banding.dat, Band.Number != 99999, 
                            !(Specie.Name %in% c("", "U", "99999"))
                            , substr(Session, nchar(Session)-1, nchar(Session)) 
                            %in% c("06", "07", "08", "09", "10", "11", "12")
                            )

# create a lookup for the sessions so that they are ordered correctly
session.lookup <- data.frame(Session = unique(banding.dat.clean$Session))
session.name.length <- nchar(as.character(session.lookup$Session))
session.lookup$session_new <- paste0(
        "session_"
        , substr(session.lookup$Session, session.name.length - 1, session.name.length)
        , "_"
        , substr(session.lookup$Session, 1, 3))

# merge lookup with the dataset 
banding.dat.clean <- merge(banding.dat.clean, session.lookup)

# Create function to run analysis in MARK --------------------------------------
# This function can be adjusted to change the models analysed by MARK
survival_analysis <- function(){
        Phi.dot <- list(formula=~1)
        Phi.habitat <- list(formula=~habitat)
        Phi.TSM <- list(formula=~TSM)
        Phi.TSMhabitat <- list(formula=~TSM+habitat)
        p.dot <- list(formula=~1)
        #p.habitat <- list(formula=~habitat)
        cml <- create.model.list("CJS")
        mark.wrapper(cml,data=sp_proc,ddl=sp_ddl,output=FALSE)
}

# Analysis for individual species ----------------------------------------------
species.list <- unique(banding.dat.clean$Specie.Name) 

# initialise list for storing results
res <- list()

for(species in species.list){
        # get records for species, get necessary columns & set the habitat type
        sp_dat <- filter(banding.dat.clean, Specie.Name == species) %>%
                select(Band.Number, session_new, Location) %>%
                mutate(habitat = ifelse(Location %in% c("LLAV", "SANA"), "Scrub",
                                        ifelse(Location == "MASE", "Native", "Introduced")))

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
        res[[species]] <- survival_analysis()
}

