# Load packages and functions --------------------------------------------------
require(Hmisc)
require(dplyr)
require(tidyr)
require(RMark)
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

# Analysis for individual species ---------------------------------------------------------

# EncounterHistory() only requires the species of interest and columns for band
# number and session
species <- "COELIGENA IRIS"
sp_dat <- filter(banding.dat.clean, Specie.Name == species) %>%
        select(Band.Number, session_new, Location) %>%
        mutate(habitat = ifelse(Location %in% c("LLAV", "SANA"), "Scrub",
                                ifelse(Location == "MASE", "Native", "Introduced")))

sp_eh <- EncounterHistory(sp_dat[,c("Band.Number", "session_new")], "session_new", "Band.Number")

sp_mark <- sp_eh$eh.mark

sp_habitat <- select(c.iris, Band.Number, habitat) %>%
        unique() %>%
        merge(c.iris.mark) %>%
        group_by(Band.Number, ch) %>%
        summarise(habitat = first(habitat))

sp_habitat$habitat <- factor(sp_habitat$habitat)

# process data for use in MARK
sp_proc <- process.data(data.frame(sp_habitat), model = "CJS", group = "habitat")
sp_ddl <- make.design.data(sp_proc)

survival_analysis <- function(){
        Phi.dot <- list(formula=~1)
        Phi.habitat <- list(formula=~habitat)
        p.dot <- list(formula=~1)
        cml <- create.model.list("CJS")
        mark.wrapper(cml,data=sp_proc,ddl=sp_ddl,output=FALSE)
}

res <- survival_analysis()

res$model.table
