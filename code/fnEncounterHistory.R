# function to create encounter histories (returns a list with a full encounter 
# history (one column per session), and one for input to mark (one field with 
# entire encounter history)). At the moment requires as input a data frame 
# (filtered to only include one species) with the band number and the session. 
# In the input session is the name of the session column and band.number is the 
# name of the unique individual identifier column
fnMarray <- function(CH){
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
        return(m.array)
}

# For the TSD model, there needs to be two M-arrays, one with the first two 
# caputres (for phi[1]) and one minus the first capture (phi[2+]). The below two
# functions split the encounter histories into these
get.first <- function(x) {
    # identify first and second non-zero value (sometimes there is only a first)
    first <- na.omit(which(x!=0)[c(1,2)])
    # set all others to zero
    x[-first] <- 0
    return(x)
}

get.other <- function(x) {
    first <- min(which(x!=0))
    x[first] <- 0
    return(x)
}

EncounterHistory <- function(data, session, band.number, group = NULL, sessions) {
        require(dplyr)
        # rename session and band number columns for use in below code
        names(data)[which(names(data)==session)]  <- "session.id"
        names(data)[which(names(data)==band.number)]  <- "band.id"
        if(is.null(group)){
                data$group.id <- 1
        } else {
                names(data)[which(names(data)==group)] <- "group.id"
        }
                
        # will require unique identifier for each row
        data$id <- 1:nrow(data)
        
        # First part creates one line per record
        eh.full <- mutate(data, count=1) %>%
                spread(session.id, count) %>%
                select(-id)
        
        # Set all NAs to zero
        eh.full[is.na(eh.full)] <- 0
        
        # some species are not sampled at all sessions - this fixes it
        no.rec.sessions <- sessions[which(!sessions %in% colnames(eh.full))]
        zeroes <- matrix(data = 0, nrow = nrow(eh.full), ncol = length(no.rec.sessions))
        colnames(zeroes) <- no.rec.sessions
        eh.full <- cbind(eh.full, zeroes)
        eh.full <- eh.full[,order(names(eh.full))] 
        
        # This part summarises and creates one line per species/band.number combo
        eh.counts <- group_by(eh.full, band.id) %>%
                summarise_each(funs(sum), matches("session")) 
        
        eh.habitat <- group_by(eh.full, band.id) %>%
                summarise(group.id = first(group.id))
        
        eh.full <- merge(eh.habitat, eh.counts)
        
        # If the same individual was captured more than once in a season, it
        # will create a number > 1 here - this code corrects it for use in MARK
        eh.full[,3:ncol(eh.full)][eh.full[,3:ncol(eh.full)] > 1] <- 1
        
        # concatenate the sessions for input to mark (i.e. as e.g. 100010100)
        eh.mark <- data.frame(select(eh.full, 1, 2),
                              ch = do.call(paste0, eh.full[,3:ncol(eh.full)]),
                              stringsAsFactors = FALSE)
        
        
        # create the TSM encounter histories
        eh.first <- data.frame(t(apply(eh.full[3:ncol(eh.full)], 1, get.first)))
        eh.other <- data.frame(t(apply(eh.full[3:ncol(eh.full)], 1, get.other)))
        cap <- rowSums(eh.full[,3:ncol(eh.full)])
        not.cap <- eh.full[which(cap < 2),3:ncol(eh.full)]
        
        # create m-array for input to WinBUGS 
        CH <- eh.full[,-1]
        m.array <- fnMarray(CH[,-1])
        CH <- split(CH[,-1], f = CH$group.id)
        m.array.gp <- lapply(CH, fnMarray)
        m.array.TSM1 <- fnMarray(eh.first)
        m.array.TSM1[,ncol(m.array.TSM1)] <- 0 
        # final col is tot recaptures - need to adjust so it doesn't include
        # those recaptured in TSM2
        m.array.TSM1[,ncol(m.array.TSM1)] <- colSums(not.cap[,1:ncol(not.cap) - 1])
        m.array.TSM2 <- fnMarray(eh.other)
                
        # rename band.number column back to what it should be
        names(eh.mark)[1] <- band.number
        names(eh.full)[1] <- band.number
        if(is.null(group)){
                eh.mark <- select(eh.mark, -group.id)
                eh.full <- select(eh.full, -group.id)
        } else {
                names(eh.mark)[2] <- group
                names(eh.full)[2] <- group
        }
        
        # TODO: add in a bit to do by group for habitat analysis
        eh <- list(eh.full = eh.full, eh.mark = eh.mark, m.array = m.array, m.array.gp = m.array.gp, m.array.TSM1 = m.array.TSM1, m.array.TSM2 = m.array.TSM2)
        return(eh)
}