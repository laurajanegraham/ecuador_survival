# Functions for creating the initial values ------------------------------------ 
#this speeds the code up a bit because it takes into account where we know the
#survival of the individual


# by creating a matrix of known dates, we save computational power (and reduce 
# time) by minimising the number of times the model needs to estimate the state
known.state.cjs <- function(ch) {
        state <- ch
        for(i in 1:dim(ch)[1]) {
                n1 <- min(which(ch[i,]==1))
                n2 <- max(which(ch[i,]==1))
                state[i, n1:n2] <- 1
                state[i, n1] <- NA
        }
        state[state==0] <- NA
        return(state)
}

# function to create initial values for latent state z incorporating the known
# values (where we know the value, the init should be NA)
cjs.init.z <- function(ch, f) {
        for (i in 1:dim(ch)[1]) {
                if(sum(ch[i,]) ==1) next
                n2 <- max(which(ch[i,]==1))
                ch[i, f[i]:n2] <- NA
        }
        for (i in 1:dim(ch)[1]) {
                ch[i, 1:f[i]] <- NA
        }
        return(ch)
}
