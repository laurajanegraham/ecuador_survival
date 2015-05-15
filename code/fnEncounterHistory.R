# function to create encounter histories (returns a list with a full encounter 
# history (one column per session), and one for input to mark (one field with 
# entire encounter history)). At the moment requires as input a data frame 
# (filtered to only include one species) with the band number and the session. 
# In the input session is the name of the session column and band.number is the 
# name of the unique individual identifier column

EncounterHistory <- function(data, session, band.number) {
        require(dplyr)
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
                summarise_each(funs(sum))
        
        # If the same individual was captured more than once in a season, it
        # will create a number > 1 here - this code corrects it for use in MARK
        eh.full[,2:ncol(eh.full)][eh.full[,2:ncol(eh.full)] > 1] <- 1
        
        # concatenate the sessions for input to mark (i.e. as e.g. 100010100)
        eh.mark <- data.frame(band.id = select(eh.full, 1), 
                              ch = do.call(paste0, eh.full[,2:ncol(eh.full)]),
                              stringsAsFactors = FALSE)
        
        
        # create m-array for input to WinBUGS 
        # TODO: add in a bit to do by group for habitat analysis
        CH <- eh.full[-1]
        nind <- nrow(CH)
        n.occasions <- ncol(CH)
        m.array <- matrix(data = 0, ncol = n.occasions + 1, nrow = n.occasions)
        m.array[, 1] <- colSums(CH)
        for (i in 1:nind){
                pos <- which(CH[i,]!=0)
                g <- length(pos)
                for (z in 1:(g - 1)) {
                        m.array[pos[z], pos[z + 1]] <- m.array[pos[z], pos[z + 1]] + 1
                }
        }
        
        for (t in 1:n.occasions) {
                m.array[t, n.occasions + 1] <- m.array[t, 1] - sum(m.array[t, 2:n.occasions])
        }
        
        m.array <- m.array[1:(n.occasions - 1), 2:(n.occasions + 1)]
        
        
        # rename band.number column back to what it should be
        names(eh.mark)[1] <- band.number
        names(eh.full)[1] <- band.number
        
        eh <- list(eh.full = eh.full, eh.mark = eh.mark, m.array = m.array)
        return(eh)
}