library("dplyr")
library("anytime")
library("lubridate")
library("chron")
library("McMasterPandemic")
library("epigrowthfit")
library("ggplot2")

##Keep a copy of the data and ICU1.csv, the params, in my home folder so that this works.
allcases <- read.csv("mikelidata.csv")
##Set the start date to whatever we want it to be.
startdate <- anytime::anydate("2020-08-15")
##Province R0s from PHAC.
##Province population sizes to use in calibration
pops <- c("AB" = 4421876, "BC" = 5147712, "ON" = 14734000, "SK" = 1178681, "MB" = 1379263, "QC" = 8574571)
##Strip out the provinces we don't want, leaving only the ones Mike used.
toKeep <- c("AB", "BC", "ON", "SK", "MB", "QC")
allcases <- allcases[allcases$Province %in% toKeep,]
##re-use my code from A1.
splitintervalCases <- lapply(as.vector(unique(allcases$Province)), function(provinceName){
  ##For each province, select all the data corresponding to that province.
  casesdf <- allcases[as.vector(allcases$Province) == provinceName,]
  ##Get rid of missing values.
  casesdf <- casesdf[!is.na(casesdf$"confirmed_positive"),]
  ##Derive interval incidence by differencing (discarding the first reported entry to align all the columns and keep them the same length), keeping the column for each provience as a check to make sure we're doing everthing properly.
  intervalcasesdf <- bind_cols("Date" = casesdf$"Date",
                               "Province" =  casesdf$"Province",
                               "Note" = casesdf$"Note",
                               "intervalCases" = c(0, diff(casesdf$"confirmed_positive", lag = 1)))
  ##Some row of contain a negative value, so let's get rid of them
  intervalcasesdf <- intervalcasesdf[intervalcasesdf$intervalCases >= 0,]
  ##Make a list of sums of weekend case reports such that the ith element is the sum of the reported cases on Saturday and Sunday for the ith weekend.
  if (provinceName == "BC" || provinceName  == "AB"){
    weekendSums <- as.list(intervalcasesdf[ format(as.Date(anytime::anydate(intervalcasesdf$Date)), '%A') == "Saturday","intervalCases"])$intervalCases +
      as.list(intervalcasesdf[format(as.Date(anytime::anydate(intervalcasesdf$Date)), '%A') == "Sunday","intervalCases"])$intervalCases
    mondayIndices <- format(as.Date(anytime::anydate(intervalcasesdf$Date)), '%A') == "Monday"
    intervalcasesdf[mondayIndices,"intervalCases"] <- c(0, weekendSums[1:length(weekendSums)-1]) + intervalcasesdf[mondayIndices,"intervalCases"]
    ##Remove the weekends.
    intervalcasesdf <- intervalcasesdf[!is.weekend(anytime::anydate(intervalcasesdf$Date)),]
  }
  else{
  }
  ##Format the data so that it can be passed to calibrate.
  colnames(intervalcasesdf) <- c("date", "a", "b", "value")
  intervalcasesdf <- intervalcasesdf[,-c(2,3)]
  intervalcasesdf$var <- rep("report", nrow(intervalcasesdf))
  intervalcasesdf$date <- anytime::anydate(intervalcasesdf$date)
  ##Grab only the second wave, using the start date we set above.
  if (provinceName == "MB"){
    startdate <- anytime::anydate("2020-09-15")
  }
  intervalcasesdf <- intervalcasesdf[intervalcasesdf$date >= startdate,]
  return(intervalcasesdf)
})
names(splitintervalCases) <- as.vector(unique(allcases$Province))
initialvalswave2 <- lapply(names(splitintervalCases), function(provinceName){
  #For each province, select all the data corresponding to that province.
  casesdf <- allcases[as.vector(allcases$Province) == provinceName,]
  ##Get rid of missing values.
  casesdf <- casesdf[!is.na(casesdf$"confirmed_positive"),]
  return(casesdf[anytime::anydate(casesdf$Date) >= startdate, "confirmed_positive"][1])
})
names(initialvalswave2) <- names(splitintervalCases)
pars <- read_params("ICU1.csv")
##Calibrate each province based on that provinces reported cases and the same set of base parameters.
calibs <- lapply(names(splitintervalCases), function(provinceName){
  provincereport <-  splitintervalCases[[provinceName]]
  pars <- update(pars, c(N = pops[[provinceName]], E0 = initialvalswave2[[provinceName]]))
  ##E0 should be the number of cumulative reports at that date.
  ##Use epigrowthfit to estimate R0 and then fix the parameter file for each province to match that estimated rate from the data.
  data(covid_generation_interval)
  pars <- fix_pars(pars, target = c(R0 = compute_R0(egf(egf_init(date = provincereport$date, cases = provincereport$value)), breaks =  covid_generation_interval$breaks, probs = covid_generation_interval$probs), Gbar = 6))
  ##So we can see what's going on.
  print(paste0("now calibrating ", provinceName))
  calibrate(base_params = pars,
             data = provincereport,
             time_args = list(break_dates = NULL),
             opt_pars = list(params = c(beta0 = 0.1)))
})
names(calibs) <- names(splitintervalCases)
##Save a calibration so we don't have to run it again to get the same results.
savecalibs <- function(){
  saveRDS(calibs, "calibs.rds")
}
##Load a calibration that we saved using the filename abve.
loadcalibs <- function(){
  return(readRDS("calibs.rds"))
}
##Plot a calibrated simulation, changing the provinceName to whatever we want it to be..
###plot(calibs$ON, data = splitintervalCases$ON, predict_args=list(keep_vars=c("report")))