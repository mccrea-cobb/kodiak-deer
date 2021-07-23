
# Import and format replicate survey data for n-mixture model


import_format_array <- function(filedir = "./data/raw_data/survey/survey_data_reps.csv",
                                herd = "Tomales", 
                                startyear = 1978, 
                                endyear = 2020){
  library(tidyverse)
  
  dat <- read.csv(filedir)
 
   dat <- dat %>%
    filter(herd == herd) %>%
    select(-c(date, herd, total, cows, bulls)) %>%
    filter(year %in% startyear:endyear) %>%
     rename(occasion = replicate)
  
  d <- data.frame(year = rep(startyear:endyear, each = 4),
                  occasion = c(1:4))
  
  dat <- dat %>%
    group_by(year) %>%
    complete(year, nesting(occasion))
  
  dat <- merge(dat, d, all=T)
  dat <- dat %>%
    #select(-year) %>%
    select(c(year, occasion, calves, est_yrlng_cows, spikes, est_adult_cows, adult_bulls))
  
  dat <- dat %>%
    gather(calves, est_yrlng_cows, spikes, est_adult_cows, adult_bulls, key = age_class, value=n)
  
  dat <- dat %>%
    spread(occasion, n)
  
  dat$age_class <- factor(dat$age_class, levels = c("calves", "est_yrlng_cows", "spikes", "est_adult_cows", "adult_bulls"))
  
  dat <- dat %>%
    arrange(year, age_class)
  
  dat <- split(dat[, 3:6], as.factor(dat$year))  # Creates a list of dataframes, one for each year
  
  C <- array(as.numeric(unlist(dat)), dim = c(nrow(dat[[1]]), ncol(dat[[1]]), length(dat)))  # Converts the list to a 3-dimensional array (1:30 plots, 1:12 annual visits, 1:29 years)
  }




