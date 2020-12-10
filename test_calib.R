#--------------------------------------------------------------------------
#   LOADING R PACKAGES
#--------------------------------------------------------------------------

## Loading Standard R Packages

library("dplyr") # For nice data frame manipulation

library("anytime") # Date/Time manipulation

library("lubridate") # Date/Time manipulation

library("chron") # Date/Time manipulation

library("ggplot2") # For fancy plots

## Loading Pandemic Modeling Packages

library("epigrowthfit") # For pandemic model parameter estimation

library("McMasterPandemic") # For COVID-19 Pandemic Simulations and Forecasting

#--------------------------------------------------------------------------
#   LOADING AND MANIPULATING CANADA COVID-19 DATA TO CALIBRATE
#--------------------------------------------------------------------------

## Keep a copy of the data (mikelidata.csv) and ICU1.csv (the params) in working directory

## Reading in provincial-level COVID-19 data

allcases <- read.csv("mikelidata.csv")

## We focus on provinces with high incidence and minimal noise

## Isolating the provinces of AB, BC, ON, SK, MB, QC

toKeep <- c("AB", "BC", "ON", "SK", "MB", "QC")

allcases <- allcases[allcases$Province %in% toKeep,]

## Setting provincial population sizes to use in calibration (source: StatCan)

pops <- c("AB" = 4421876, 
          "BC" = 5147712, 
          "ON" = 14734000, 
          "SK" = 1178681, 
          "MB" = 1379263, 
          "QC" = 8574571)

## Setting start date of calibration. Date signifies rough start of 2nd wave.

startdate <- anytime::anydate("2020-08-15")

## Creating a function to perform the following tasks:
# 1) Matches real data to model variables (report, deaths and hospitalization)
# 2) Creates non-negative interval incidence and death columns from cumulative data (hospitalization already recorded as daily count)
# 3) Removes missing values
# 4) Converts data to "long form" so that it can be passed to calibrate()
# 5) Filters data so that it starts at the beginning of the second wave

## Matching real data to model variables

splitinterval <- function(var){
  if (var == "report"){
    ## colselect grabs the data from the main data frame and var tells the calibration function what variable the data is.
    colselect <- "confirmed_positive"
    var <- "report"
  }
  else if (var == "deaths"){
    colselect <- "deceased"
    var <- "death"
  }
  else{
    colselect <- "Hospitalization"
    var <-  "H"
  }
  lapply(as.vector(unique(allcases$Province)), function(provinceName){
    ## Select all the data corresponding to each province
    casesdf <- allcases[as.vector(allcases$Province) == provinceName,]
    ## Remove missing values
    casesdf <- casesdf[!is.na(casesdf[,colselect]),]
    ## Create interval incidence column by differencing 
    # Discard the first reported entry to align all the columns and keep them the same length, keeping the column for each province as a check to make sure we're doing everything properly.
    intervalcasesdf <- bind_cols("Date" = casesdf$"Date",
                                 "Province" =  casesdf$"Province",
                                 "Note" = casesdf$"Note",
                                 "intervalCases" = c(0, diff(casesdf[,colselect], lag = 1)))
    ## Removal of negative values
    intervalcasesdf <- intervalcasesdf[intervalcasesdf$intervalCases >= 0,]
    ## Converting data to "long form" so that it can be passed to calibrate
    colnames(intervalcasesdf) <- c("date", "a", "b", "value")
    intervalcasesdf <- intervalcasesdf[,-c(2,3)]
    intervalcasesdf$var <- rep(var, nrow(intervalcasesdf))
    intervalcasesdf$date <- anytime::anydate(intervalcasesdf$date)
    ## Filtering data to begin at second wave (using startdate)
    if (provinceName == "MB"){
      startdate <- anytime::anydate("2020-09-15")
    }
    intervalcasesdf <- intervalcasesdf[intervalcasesdf$date >= startdate,]
    return(intervalcasesdf)
  })
}

## Create vectors of observed data.

splitintervalCases <- splitinterval(var = "report")

splitintervaldeaths <- splitinterval(var = "deaths")

splitintervalhosp <- splitinterval(var = "hosp")

## Setting the names in the function stops it from returning the list which is what we want, so we'll do it here.

names(splitintervalhosp) <- names(splitintervalCases) <- names(splitintervaldeaths) <- unique(allcases$Province)  

#--------------------------------------------------------------------------
#   CALIBRATING
#--------------------------------------------------------------------------

## Reading in preset parameter file from McMasterPandemic

pars <- read_params("ICU1.csv")

## Update iso_m and iso_s to be non-zero for social distancing scenario to work

pars <- update(pars, iso_m = 0.000001, iso_s = 0.000001)

## Creation of function "calibrate_good" which calibrates all provinces

calibrate_good <- function(){
  goodcalibs <- lapply(names(splitintervalCases), function(provinceName){
    dd <-  bind_rows(splitintervalCases[[provinceName]], splitintervaldeaths[[provinceName]])
    ## Order by date
    dd <- dd[order(anytime::anydate(dd$date)),]
    provincereport <- splitintervalCases[[provinceName]]
    ## Use epigrowthfit to estimate R0 and then fix the parameter file for each province to match the estimated rate from the data and setting Gbar = 6
    data(covid_generation_interval)
    pars <- fix_pars(pars, target = c(R0 = compute_R0(egf(egf_init(date = provincereport$date, cases = provincereport$value)), breaks =  covid_generation_interval$breaks, probs = covid_generation_interval$probs), Gbar = 6))
    ## Set population size for each province
    pars <- update(pars, c(N = pops[[provinceName]]))
    init_e0 <- provincereport$value[[16]]
    if (init_e0 == 0){
      loginit_e0 <- log(provincereport$value[[18]])
    }
    else if (provinceName == "SK"){
      loginit_e0 <-2
    }
    else{
      loginit_e0 <- log(provincereport$value[[16]])
    }
    ## Parameters to be optimized are based on the calibration to Ontario by Michael Li in 
    # https://github.com/bbolker/McMasterPandemic/blob/master/ontario/Ontario_current.R
    if (provinceName == "QC") {
      optpars <- list(params = c(log_E0 = loginit_e0, log_beta0 = -1, logit_phi1 = -1), logit_rel_beta0 = c(-0.223), log_nb_disp=c(report=1, death=1, H=1))
      bks <- c("2020-Oct-06")
    } else{
      optpars <- list(params = c(log_E0 = loginit_e0, log_beta0 = -1, logit_phi1 = -1), log_nb_disp=c(report=1, death=1, H=1))
      bks <- NULL
    }
    ## Calibration to SK deaths is not done due to noise.
    if (provinceName == "SK"){
      dd <- dd[dd$var == "report",]
    }
    else{
    }
    ## Print which province is currently being calibrated
    print(paste0("now calibrating ", provinceName))
    ## Calibration step
      calibrate(base_params = pars,
                data = dd,
                debug_plot = FALSE,
                ## sim_args below based on the calibration to Ontario in 
                # https://github.com/bbolker/McMasterPandemic/blob/master/ontario/Ontario_current.R
                sim_args = list(ndt = 1, ratemat_args = list(testing_time = "report")),
                time_args = list(break_dates = bks),
                opt_pars = optpars)
  })
  names(goodcalibs) <- names(splitintervalCases)
  return(goodcalibs)
}
##Our method that is not as good.
##Calibrate each province based on that provinces reported cases and the same set of base parameters.
calibrate_bad <- function(){
  badcalibs <- lapply(names(splitintervalCases), function(provinceName){
    dd <-  bind_rows(splitintervalCases[[provinceName]], splitintervaldeaths[[provinceName]])
    ##We need to order by date
    dd <- dd[order(anytime::anydate(dd$date)),]
    if (provinceName == "SK"){
      dd <- dd[dd$var == "report",]
    }
    else{
    }
    provincereport <- splitintervalCases[[provinceName]]
    ##Use epigrowthfit to estimate R0 and then fix the parameter file for each province to match that estimated rate from the data.
    data(covid_generation_interval)
    pars <- fix_pars(pars, target = c(R0 = compute_R0(egf(egf_init(date = provincereport$date, cases = provincereport$value)), breaks =  covid_generation_interval$breaks, probs = covid_generation_interval$probs), Gbar = 6))
    pars <- update(pars, c(N = pops[[provinceName]]))
    ##So we can see what's going on.
    print(paste0("now calibrating ", provinceName))
    calibrate(base_params = pars,
              data = dd,
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

#--------------------------------------------------------------------------
#   FORECAST FUNCTION WITH SCENARIOS
#--------------------------------------------------------------------------

## Creation of a forecasting function "forecast_province" with the following arguments
# 1) Name of province to forecast
# 2) Start date of forecast
# 3) End date of forecast
# 4) Scenario number (see comment below for details)
# 5) List of calibrations

## There are four scenarios to choose from
# scenario = 1: No government intervention (status quo is maintained)
# scenario = 2: Strict lockdown imposed in Canada on December 18th, 2020, and then relaxed six weeks later
# scenario = 3: ICUs fill up on Jan 15th (a month into the simulation), and then clear exactly a month later
# scenario = 4: Social distancing
# Note: initial iso_m and iso_s are very small (1/1mill). Large relative values are needed

forecast_province <- function(provinceName, 
                              sd = anytime::anydate("2020-08-01"), 
                              ed = anytime::anydate("2021-12-18"), 
                              scenario = 1, 
                              calibsList = goodcalibs,
                              lockdown_beta0 = 0.5,
                              lockdown_relax = 0.1,
                              phi2_1 = 0.1,
                              phi2_2 = 0.5,
                              isom_init = 500000,
                              isos_init = 700000){
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
    time_pars <- NULL
  }
  else if (scenario == 2){
    ##Lockdown then relax six weeks after.
    ##Six weeks after is Jan 29th, 2021
    time_pars <- data.frame(Date=c("2020-12-18", "2021-01-29"),
                            Symbol=c("beta0", "beta0"),
                            Relative_value=c(lockdown_beta0, lockdown_relax),
                            stringsAsFactors=FALSE)
  }
  else if (scenario == 3){
    ##ICU's fill up and then clear exactly a month later.
    time_pars <- data.frame(Date=c("2020-01-15", "2021-02-15"),
                            Symbol=c("phi2", "phi2"),
                            Relative_value=c(phi2_1, phi2_2),
                            stringsAsFactors=FALSE)
  }
  else if (scenario == 4){
    ##Social distancing protocols are implemented.
    time_pars <- data.frame(Date=c("2020-12-18", "2020-12-18", "2021-02-18", "2021-02-18"),
                            Symbol=c("iso_m", "iso_s", "iso_m", "iso_s"),
                            Relative_value=c(isom_init,isos_init, 1, 1),
                            stringsAsFactors=FALSE)
  }
  else{
  }
  sim <- run_sim(pars, start_date = sd, end_date = ed, params_timevar = time_pars)
  return(sim)
}
##Plot a calibrated simulation, changing the provinceName to whatever we want it to be..
test_calib_plot <- function(provinceName, calibslist = goodcalibs, reportlist = splitintervalCases, deathlist = splitintervaldeaths, hosplist = splitintervalhosp){
  dd <-  bind_rows(splitintervalCases[[provinceName]], splitintervaldeaths[[provinceName]])
  ##We might need to order by date, so it's safer to.
  dd <- dd[order(anytime::anydate(dd$date)),]
  plot(calibslist[[provinceName]], data = dd, predict_args=list(keep_vars=c("report", "death")))
}
##Test a forecast by plotting it on the same graph as the observed report data and the calibrate to it.
test_forecast_plot <- function(provinceName, sim = forecast_province(provinceName)){
  plot(sim, drop_states = c("S", "R", "I", "cumRep", "E", "X", "D", "incidence", "ICU", "H")) +
    ##Add observed reported cases to the plot
    geom_point(data = splitintervalCases[[provinceName]],
               mapping = aes(x = date, y = value)) +
    ##Add observed deaths to the plot
    geom_point(data = splitintervaldeaths[[provinceName]],
               mapping = aes(x = date, y = value)) 
  #+
  ##Add observed hospitalizations to the plot
  #geom_point(data = splitintervalhosp[[provinceName]],
  #           mapping = aes(x = date, y = value))
}

