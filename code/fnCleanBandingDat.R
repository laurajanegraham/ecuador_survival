CleanBandingDat <- function() {
        require(Hmisc)
        require(dplyr)
        require(tidyr)
        # read in the table from full_data.mdb and write out to banding_sheet.csv - this
        # gets rid of column attribute issues
        banding.dat <- mdb.get("data/full_data.mdb", tables="Banding Sheet")
        write.csv(banding.dat, file="data/banding_sheet.csv")
        
        # re-read in data
        banding.dat <- read.csv("data/banding_sheet.csv", stringsAsFactors = FALSE)
        
        # remove records with no species, session or band (at the moment this
        # also limits the temporal extent to that used in the chapter by Boris)
        banding.dat.clean <- filter(banding.dat, Band.Number != 99999, 
                                    !(Specie.Name %in% c("", "U", "99999"))
                                    , substr(Session, nchar(Session)-1, nchar(Session)) 
                                    %in% c("06", "07", "08", "09", "10", "11", "12")
                                    ) %>%
                # also add on a habitat column based on location
                mutate(habitat = ifelse(Location %in% c("LLAV", "SANA"), "Scrub",
                                        ifelse(Location == "MASE", "Native", "Introduced")))
                
        
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
        return(banding.dat.clean)
}