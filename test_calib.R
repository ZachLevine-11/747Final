##standard R packaes.
library("dplyr")
library("anytime")
library("lubridate")
library("chron")
library("ggplot2")
##McMaster -specific package
library("McMasterPandemic")
library("epigrowthfit")

##Keep a copy of the data and ICU1.csv, the params, in my home folder so that this works.
allcases <- read.csv("mikelidata.csv")
##Set the start date to whatever we want it to be.
startdate <- anytime::anydate("2020-08-15")
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

pars <- read_params("ICU1.csv")

##Mike Li method.
calibrate_good <- function(){
  goodcalibs <- lapply(names(splitintervalCases), function(provinceName){
    provincereport <-  splitintervalCases[[provinceName]]
    ##E0 should be the number of cumulative reports at that date.
    ##Use epigrowthfit to estimate R0 and then fix the parameter file for each province to match that estimated rate from the data.
    data(covid_generation_interval)
    pars <- fix_pars(pars, target = c(R0 = compute_R0(egf(egf_init(date = provincereport$date, cases = provincereport$value)), breaks =  covid_generation_interval$breaks, probs = covid_generation_interval$probs), Gbar = 6))
    pars <- update(pars, c(N = pops[[provinceName]]))
    ##So we can see what's going on.
    print(paste0("now calibrating ", provinceName))
    ##Rigged based on the calibration to Ontario in https://github.com/bbolker/McMasterPandemic/blob/master/ontario/Ontario_current.R
    optpars <- list(params = c(log_E0 = 2, log_beta0 = -1, logit_phi1 = -1), log_nb_disp=c(report=1, death=1, H=1)) 
    calibrate(base_params = pars,
              data = provincereport,
              debug_plot = FALSE,
              ##based on the calibration to ontario in https://github.com/bbolker/McMasterPandemic/blob/master/ontario/Ontario_current.R
              sim_args = list(ndt = 1, ratemat_args = list(testing_time = "report")),
              time_args = list(break_dates = NULL),
              opt_pars = optpars)
  })
  names(goodcalibs) <- names(splitintervalCases)
  return(goodcalibs)
}

##Save a calibration so we don't have to run it again to get the same results.

##Our method that is not as goood.
##Calibrate each province based on that provinces reported cases and the same set of base parameters.
calibrate_bad <- function(){
  badcalibs <- lapply(names(splitintervalCases), function(provinceName){
    provincereport <-  splitintervalCases[[provinceName]]
    pars <- update(pars, c(N = pops[[provinceName]]))
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
  names(badcalibs) <- names(splitintervalCases)
  return(badcalibs)
}

##Save a calibration so we don't have to run it again to get the same results.
savecalibs <- function(){
  saveRDS(goodcalibs, "calibs.rds")
}
##Load a calibration that we saved using the filename above. Make sure calibs.rds is in your working directory.
loadcalibs <- function(){
  return(readRDS("calibs.rds"))
}
##Plot a calibrated simulation, changing the provinceName to whatever we want it to be..
###plot(calibs$ON, data = splitintervalCases$ON, predict_args=list(keep_vars=c("report")))

##Do forecasts for each province.
##Options
##scenario =  1: status quo. Do nothing
##scenario = 2: Strict lockdown is imposed throughout the country on December 18th, 2020, and then relaxed six weeks lated
##scenario = 3: ICU's fill up on Jan 15th (a month into the simulation), and then clear exactly a month later.
forecast_province <- function(provinceName, sd = anytime::anydate("2020-08-01"), ed = anytime::anydate("2021-12-18"), scenario = 1, calibsList = goodcalibs){
  calib <- calibsList[[provinceName]]
  ff <- calib$forecast_args
  pars <- ff$base_params
  ##Do these manually I guess.
  #pars[["E0"]] <- coef(calib, "fitted")$params["E0"]
  #pars[["beta0"]] <- coef(calib, "fitted")$params["beta0"]
  #pars[["phi1"]] <- coef(calib, "fitted")$params["phi1"]
  ##The below line does the same thing. Lee, I'm leaving the above three in there so you can see.
  pars[names(coef(calib, "fitted")$params)] <- coef(calib, "fitted")$params
  ##Stil don't now why this didn't work.
  #pars <- update(pars, ff$opt_pars$params)
  ##The only thing that changes between scenarios is the time_pars.
  if (scenario == 1){
    ##Status quo
    time_pars <- data.frame(Date=c("2020-12-01"),
                            Symbol=c("beta0"),
                            Relative_value=c(1),
                            stringsAsFactors=FALSE)
  }
  else if (scenario == 2){
    ##Lockdown then relax six weeks after.
    ##Six weeks after is Jan 29th, 2021
    time_pars <- data.frame(Date=c("2020-12-18", "2021-01-29"),
                            Symbol=c("beta0", "beta0"),
                            Relative_value=c(0.5, 0.1),
                            stringsAsFactors=FALSE)
  }
  else if (scenario == 3){
    ##ICU's fill up and then clear exactly a month later.
    time_pars <- data.frame(Date=c("2021-01-15", "2021-02-15"),
                            Symbol=c("phi2", "beta0"),
                            Relative_value=c(0.1, 1),
                            stringsAsFactors=FALSE)
  }
  else{
  }
  sim <- run_sim(pars, start_date = sd, end_date = ed, params_timevar = time_pars)
  return(sim)
}

##Test a forecast by plotting it on the same graph as the observed report data and the calibrate to it.
test_forecast_plot <- function(provinceName, sim = forecast_province(provinceName)){
  plot(sim, drop_states = c("S", "R", "I", "cumRep", "E", "X", "incidence", "D", "H", "ICU")) +
    geom_point(data = splitintervalCases[[provinceName]],
               mapping = aes(x = date, y = value))

}