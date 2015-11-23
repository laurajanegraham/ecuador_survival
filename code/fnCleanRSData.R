fnLocation <- function(df) {
    factor(df$Id, labels=c("SANA", "LLAV", "MASE", "MAIN"))
}

fnTimeWindowAllSites <- function(df, days, samp.dates) {
    df <- merge(df, samp.dates)
    df$date.interval <- interval(df$end - days, df$end)
    df$interval.boolean <- as.Date(df$date) %within% df$date.interval
    out <- select(df, site, session, measure, mean, date, interval.boolean) %>%
        filter(interval.boolean) %>%
        group_by(session, measure) %>%
        summarise(mean.val=mean(mean, na.rm=TRUE))
    return(out)
}

CleanRSData <- function(window) {
    # load the sample dates to create windows from
    samp.dates <- read.csv("data/sampling_dates_ec.csv") %>%
        mutate(session=substr(rownames(.), 1, 10)) %>%
        arrange(session) %>%
        mutate(beg=as.Date(beg), end=as.Date(end)) 

    # load RS data files
    evi <- do.call("rbind", lapply(list.files("data/rs-data/", pattern="EVI", full.names = TRUE), 
                                   function(x) read.csv(x, stringsAsFactors = FALSE))) %>%
        filter(complete.cases(.)==TRUE) %>%
        mutate(site=fnLocation(.),
               date=as.Date(gsub("_", "-", substr(system.index, 13, 22))),
               measure="EVI") %>%
        select(site, date, measure, mean) # 825 missing values (52%)
    
    ndvi <- do.call("rbind", lapply(list.files("data/rs-data/", pattern="NDVI", full.names = TRUE), 
                                    function(x) read.csv(x, stringsAsFactors = FALSE))) %>%
        filter(complete.cases(.)==TRUE) %>%
        mutate(site=fnLocation(.),
               date=as.Date(gsub("_", "-", substr(system.index, 13, 22))),
               measure="NDVI") %>%
        select(site, date, measure, mean) # 765 missing values (48%)
    
    temp <- do.call("rbind", lapply(list.files("data/rs-data/", pattern="temp", ignore.case=TRUE, full.names = TRUE), 
                                    function(x) read.csv(x, stringsAsFactors = FALSE))) %>%
        filter(complete.cases(.)==TRUE) %>%
        mutate(site=fnLocation(.),
               date=as.Date(paste(substr(system.index, 12, 15), 
                                  substr(system.index, 16, 17),
                                  substr(system.index, 18, 19), sep="-")),
               measure="Temperature") %>%
        select(site, date, measure, mean)
    
    # create means for the time window
    RS.dat <- fnTimeWindowAllSites(evi, window * 7, samp.dates)
    RS.dat$evi <- RS.dat$mean.val
    RS.dat$ndvi <- fnTimeWindowAllSites(ndvi, window * 7, samp.dates)$mean.val
    RS.dat$temp <- fnTimeWindowAllSites(temp, window * 7, samp.dates)$mean.val
    
    # rename the session so can be ordered
    RS.dat$session <- paste0(
        substr(RS.dat$session, nchar(RS.dat$session) - 1, nchar(RS.dat$session))
        , "_", substr(RS.dat$session, 1, 3))
    
    RS.dat <- arrange(RS.dat, session) %>%
        select(session, evi, ndvi, temp)
}
    
    