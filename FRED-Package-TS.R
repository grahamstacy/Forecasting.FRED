library(dplyr)

get_fred <- function(fred_api_key, 
                     fred_series_id,
                     units = "pch",
                     aggregation_method = "avg",
                     frequency = NULL) {
  
  fredr::fredr_set_key(fred_api_key)
  
  tibble.get_fred <- fredr::fredr(series_id = fred_series_id, 
                                   frequency = frequency, 
                                   units = units, 
                                   aggregation_method = aggregation_method)
  
  tibble.get_fred.year <- as.numeric( substr(tibble.get_fred$date,1,4) )
  tibble.get_fred.month <- as.numeric( substr(tibble.get_fred$date,6,7) )
  tibble.get_fred.day <- as.numeric( substr(tibble.get_fred$date,9,10) )
  
  if ( all( unique(tibble.get_fred.day) %in% c(1) ) ) {
    invisible()
  } else {
    stop("Error: The FRED series selected is not monthly, quarterly, or annual.")
  }
  
  if ( all( unique(tibble.get_fred.month) %in% c(1) ) ) {
    
    ts.get_fred <- ts(tibble.get_fred$value, 
                      start = tibble.get_fred.year[1],
                      frequency = 1)
    
  } else if ( all( unique(tibble.get_fred.month) %in% c(1,4,7,10) ) ) {
    
    tibble.get_fred.month[tibble.get_fred.month==4] <- 2
    tibble.get_fred.month[tibble.get_fred.month==7] <- 3
    tibble.get_fred.month[tibble.get_fred.month==10] <- 4
    ts.get_fred <- ts(tibble.get_fred$value, 
                      start = c(tibble.get_fred.year[1],
                                tibble.get_fred.month[1]),
                      frequency = 4)
    
  } else if ( all( unique(tibble.get_fred.month) %in% c(1:12) ) ) {
    
    ts.get_fred <- ts(tibble.get_fred$value, 
                      start = c(tibble.get_fred.year[1],
                                tibble.get_fred.month[1]),
                      frequency = 12)
    
  } else {
    stop("Error: The FRED series selected is not monthly, quarterly, or annual.")
  }
  
  ts.get_fred <- tseries::na.remove(ts.get_fred)
  
  return(ts.get_fred)
  
}

fcast <- function(ts,
                  methods = c("a","e","m","n","sn"), 
                  h = 4,
                  nsim= 1000) {
  
  arima.fcast <- function(ts, h, nsim) {
    
    matrix.arima.fcast <- matrix(0:0, nrow=1,ncol=h)
    fit.arima.fcast <- forecast::auto.arima(ts,
                                            ic = "aicc")
    fitted.arima.fcast <- as.numeric(fit.arima.fcast$fitted)
    
    for (i in 1:nsim) {
      
      sim.arima.fcast <- stats::simulate(fit.arima.fcast, nsim = h)
      matrix.arima.fcast <- rbind(matrix.arima.fcast, sim.arima.fcast)
      
    }
    
    row.names(matrix.arima.fcast) <- c(0:nsim)
    
    return(list(matrix.arima.fcast, fitted.arima.fcast))
    
  }
  
  ets.fcast <- function(ts, h, nsim) {
    
    matrix.ets.fcast <- matrix(0:0, nrow=1,ncol=h)
    fit.ets.fcast <- forecast::ets(ts,
                                   allow.multiplicative.trend = TRUE,
                                   ic = "aicc")
    fitted.ets.fcast <- as.numeric(fit.ets.fcast$fitted)
    
    for (i in 1:nsim) {
      
      sim.ets.fcast <- stats::simulate(fit.ets.fcast, nsim = h)
      matrix.ets.fcast <- rbind(matrix.ets.fcast, sim.ets.fcast)
      
    }
    
    row.names(matrix.arima.fcast) <- c(0:nsim)
    
    return(list(matrix.arima.fcast, fitted.arima.fcast))
    
  }
  
  mean.fcast <- function(ts, h, nsim) {
    
    matrix.mean.fcast <- matrix(0:0, nrow=1,ncol=h)
    fitted.mean.fcast <- c()
    for (i in 1:(length(ts))) {
      fitted.mean.fcast <- c(fitted.mean.fcast, mean(as.numeric(ts)))
    }
    
    for (i in 1:nsim) {
      
      initial.mean.fcast <- rnorm(1, mean = mean(as.numeric(ts)), sd = sd(as.numeric(ts)))
      sim.mean.fcast <- c(initial.mean.fcast)
      total.mean.fcast <- c(as.numeric(ts), initial.mean.fcast)
      
      for (j in 1:(h-1)) {
        
        additional.mean.fcast <- rnorm(1, 
                                       mean = mean(total.mean.fcast),
                                       sd = sd(total.mean.fcast))
        sim.mean.fcast <- c(sim.mean.fcast, additional.mean.fcast)
        total.mean.fcast<- c(total.mean.fcast, additional.mean.fcast)
        
      }
      
      matrix.mean.fcast <- rbind(matrix.mean.fcast, sim.mean.fcast)
      
    }
    
    row.names(matrix.mean.fcast) <- c(0:nsim)
    
    return(list(matrix.mean.fcast, fitted.mean.fcast))
    
  }

  naive.fcast <- function(ts, h, nsim) {
    
    matrix.naive.fcast <- matrix(0:0, nrow=1,ncol=h)
    fit.naive.fcast <- forecast::Arima(ts, order = c(0,1,0), seasonal = c(0,0,0))
    fitted.naive.fcast <- as.numeric(fit.naive.fcast$fitted)
    
    for (i in 1:nsim) {
      
      sim.naive.fcast <- stats::simulate(fit.naive.fcast, nsim = h)
      matrix.naive.fcast <- rbind(matrix.naive.fcast, sim.naive.fcast)
      
    }
    
    row.names(matrix.naive.fcast) <- c(0:nsim)
    
    return(list(matrix.naive.fcast, fitted.naive.fcast))
    
  }
  
  snaive.fcast <- function(ts, h, nsim) {
    
    matrix.snaive.fcast <- matrix(0:0, nrow=1,ncol=h)
    fit.snaive.fcast <- forecast::Arima(ts, order = c(0,0,0), seasonal = c(0,1,0))
    fitted.snaive.fcast <- as.numeric( fit.snaive.fcast$fitted )
    
    for (i in 1:nsim) {
      
      sim.snaive.fcast <- stats::simulate(fit.snaive.fcast, nsim = h)
      matrix.snaive.fcast <- rbind(matrix.snaive.fcast, sim.snaive.fcast)
      
    }
    
    row.names(matrix.snaive.fcast) <- c(0:nsim)
    
    return(list(matrix.snaive.fcast, fitted.snaive.fcast))
    
  }
  
  rwf.fcast <- function(ts, h, nsim) {
    
    matrix.rwf.fcast <- matrix(0:0, nrow=1,ncol=h)
    fit.rwf.fcast <- forecast::Arima(ts, order = c(0,1,0), seasonal = c(0,0,0), include.drift = TRUE)
    fitted.rwf.fcast <- as.numeric(fit.rwf.fcast$fitted)
    
    for (i in 1:nsim) {
      
      sim.rwf.fcast <- stats::simulate(fit.rwf.fcast, nsim = h)
      matrix.rwf.fcast <- rbind(matrix.rwf.fcast, sim.rwf.fcast)
      
    }
    
    row.names(matrix.rwf.fcast) <- c(0:nsim)
    
    return(list(matrix.rwf.fcast, fitted.rwf.fcast))
    
  }
  
  nnetar.fcast <- function(ts, h, nsim) {
    
    matrix.nnetar.fcast <- matrix(0:0, nrow=1,ncol=h)
    fit.nnetar.fcast <- forecast::nnetar(ts)
    fitted.nnetar.fcast <- as.numeric(fit.nnetar.fcast$fitted)
    
    for (i in 1:nsim) {
      
      sim.nnetar.fcast <- stats::simulate(fit.nnetar.fcast, nsim = h)
      matrix.nnetar.fcast <- rbind(matrix.nnetar.fcast, sim.nnetar.fcast)
      
    }
    
    row.names(matrix.nnetar.fcast) <- c(0:nsim)
    
    return(list(matrix.nnetar.fcast, fitted.nnetar.fcast))
    
  }
  
}

x <- get_fred(fred_api_key = "0962c1d1b6316b7c4b6c68ced8211ac9", 
              fred_series_id = "GDP")