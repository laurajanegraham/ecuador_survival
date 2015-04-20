require(Hmisc)
require(dplyr)
require(tidyr)

# load the function to create encounter history
source("code/fnEncounterHistory.R")

# read in the table from full_data.mdb and write out to banding_sheet.csv - this
# gets rid of column attribute issues
banding.dat <- mdb.get("data/full_data.mdb", tables="Banding Sheet")
write.csv(banding.dat, file="data/banding_sheet.csv")

# re-read in data
banding.dat <- read.csv("data/banding_sheet.csv")

# remove records with no species, session or band
banding.dat.clean <- filter(banding.dat, Band.Number != 99999, 
                            !(Specie.Name %in% c("", "U", "99999")), 
                            Session != "")

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

# 
c.iris <- filter(banding.dat.clean, Specie.Name == "COELIGENA IRIS") %>%
        select(Band.Number, session_new)

