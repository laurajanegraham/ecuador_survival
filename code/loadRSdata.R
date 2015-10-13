require(dplyr)
require(lubridate)
require(ggplot2)

# functions required
fnLocation <- function(df) {
    factor(df$Id, labels=c("SANA", "LLAV", "MASE", "MAIN"))
}

fnTimeWindow <- function(df, days) {
    df <- merge(df, samp.dates)
    df$date.interval <- interval(df$end - days, df$end)
    df$interval.boolean <- as.Date(df$date) %within% df$date.interval
    out <- select(df, site, session, measure, mean, date, interval.boolean) %>%
        filter(interval.boolean) %>%
        group_by(site, session, measure) %>%
        summarise(mean.val=mean(mean, na.rm=TRUE),
                  SE.val=sd(mean, na.rm=TRUE)/sqrt(n())) %>%
        mutate(window.size=days/7)
    return(out)
}

fnTimeWindowAllSites <- function(df, days) {
    df <- merge(df, samp.dates)
    df$date.interval <- interval(df$end - days, df$end)
    df$interval.boolean <- as.Date(df$date) %within% df$date.interval
    out <- select(df, site, session, measure, mean, date, interval.boolean) %>%
        filter(interval.boolean) %>%
        group_by(session, measure) %>%
        summarise(mean.val=mean(mean, na.rm=TRUE),
                  SE.val=sd(mean, na.rm=TRUE)/sqrt(n())) %>%
        mutate(window.size=days/7)
    return(out)
}

# load the sample dates to create windows from
samp.dates <- read.csv("data/sampling_dates_ec.csv") %>%
    mutate(session=substr(rownames(.), 1, 10)) %>%
    arrange(session) %>%
    mutate(beg=as.Date(beg), end=as.Date(end)) 

# get the RS data files
files <- list.files("data/rs-data/", pattern=".csv", full.names = TRUE)

test <- 
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

# create a dataframe with the mean RS measurements for 2, 6 and 12 week windows
wk <- c(2, 4, 6, 8, 10, 12)
window.dat <- list()
window.dat.allsites <- list()
for(w in wk) {
    window.dat[[paste0("evi", w)]] <- fnTimeWindow(evi, w * 7)
    window.dat[[paste0("ndvi", w)]] <- fnTimeWindow(ndvi, w * 7)
    window.dat[[paste0("temp", w)]] <- fnTimeWindow(temp, w * 7)
    
    window.dat.allsites[[paste0("evi", w)]] <- fnTimeWindowAllSites(evi, w * 7)
    window.dat.allsites[[paste0("ndvi", w)]] <- fnTimeWindowAllSites(ndvi, w * 7)
    window.dat.allsites[[paste0("temp", w)]] <- fnTimeWindowAllSites(temp, w * 7)
}

window.dat <- do.call("rbind", window.dat)
window.dat.allsites <- do.call("rbind", window.dat.allsites)
# rename sessions so they display in the correct order
window.dat$session <- paste0(
    substr(window.dat$session, nchar(window.dat$session) - 1, nchar(window.dat$session))
    , "_", substr(window.dat$session, 1, 3))

window.dat.allsites$session <- paste0(
    substr(window.dat.allsites$session, nchar(window.dat.allsites$session) - 1, nchar(window.dat.allsites$session))
    , "_", substr(window.dat.allsites$session, 1, 3))

# create plots of mean RS measurements
window.dat.plot <- subset(window.dat, window.size %in% c(4, 8, 12))
ggplot(window.dat.plot, aes(x = session, y=mean.val, group=site, colour=site)) + 
    geom_line() + 
    geom_point() + 
    geom_errorbar(aes(ymin=mean.val-SE.val, ymax=mean.val+SE.val), width=.1) +
    facet_grid(measure~window.size, scales = "free_y") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("figures/RS_data_by_site.png", width=12, height=6)

ggplot(window.dat.allsites,
       aes(x = session, y=mean.val, group=as.factor(window.size), colour=as.factor(window.size))) + 
    geom_line() + 
    geom_point() + 
    geom_errorbar(aes(ymin=mean.val-SE.val, ymax=mean.val+SE.val), width=.1) +
    facet_grid(measure~., scales = "free_y") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("figures/RS_data_overall.png", width=12, height=6)