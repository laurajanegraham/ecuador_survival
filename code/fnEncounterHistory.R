# function to create encounter histories (returns a list with a full encounter 
# history (one column per session), and one for input to mark (one field with 
# entire encounter history)). At the moment requires as input a data frame 
# (filtered to only include one species) with the band number and the session. 
# In the input session is the name of the session column and band.number is the 
# name of the unique individual identifier column

EncounterHistory <- function(data, session, band.number) {
        
        # rename session and band number columns for use in below code
        names(data)[which(names(data)==session)]  <- "session.id"
        names(data)[which(names(data)==band.number)]  <- "band.id"
        
        # will require unique identifier for each row
        data$id <- 1:nrow(data)
        
        # First part creates one line per record
        eh.full <- mutate(data, count=1) %>%
                spread(session.id, count) 
        
        # drop id column
        eh.full <- eh.full[, !names(eh.full) %in% "id"]
        
        # Set all NAs to zero
        eh.full[is.na(eh.full)] <- 0
                
        # This part summarises and creates one line per species/band.number combo
        eh.full <- group_by(eh.full, band.id) %>%
                summarise_each(funs(first))
        
        # concatenate the sessions for input to mark (i.e. as e.g. 100010100)
        eh.mark <- data.frame(band.id = select(eh.full, 1), 
                              ch = do.call(paste0, eh.full[,2:ncol(eh.full)]),
                              stringsAsFactors = FALSE)
        
        # rename band.number column back to what it should be
        names(eh.mark)[1] <- band.number
        
        eh <- list(eh.full = eh.full, eh.mark = eh.mark)
}