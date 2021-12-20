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
                  methods = c("a","e","m","n","s","rw","nn"), 
                  h = 4,
                  nsim= 1000,
                  conf = 0.90,
                  tscv = FALSE,
                  preds = 8,
                  boxcox = FALSE) {
  
  ts0 <- ts
  methods0 <- methods
  h0 <- h
  nsim0 <- nsim
  conf0 <- conf
  tscv0 <- tscv
  preds0 <- preds
  
  if (boxcox = TRUE) {
    lambda0 <- forecast::BoxCox.lambda(ts0)
    ts0 <- forecast::BoxCox(ts0, lambda0)
  }
  
  arima.fcast <- function(ts, h, nsim) {
    
    matrix.arima.fcast <- matrix(0:0, nrow=0,ncol=h)
    fit.arima.fcast <- forecast::auto.arima(ts,
                                            ic = "aicc")
    fitted.arima.fcast <- as.numeric(fit.arima.fcast$fitted)
    accuracy.arima.fcast <- accuracy(fit.arima.fcast)[1:6]
    
    for (i in 1:nsim) {
      
      sim.arima.fcast <- stats::simulate(fit.arima.fcast, nsim = h)
      matrix.arima.fcast <- rbind(matrix.arima.fcast, sim.arima.fcast)
      
    }
    
    row.names(matrix.arima.fcast) <- c(1:nsim)

    return(list(matrix.arima.fcast, fitted.arima.fcast, accuracy.arima.fcast))
    
  }
  
  ets.fcast <- function(ts, h, nsim) {
    
    matrix.ets.fcast <- matrix(0:0, nrow=0,ncol=h)
    fit.ets.fcast <- forecast::ets(ts,
                                   allow.multiplicative.trend = TRUE,
                                   ic = "aicc")
    fitted.ets.fcast <- as.numeric(fit.ets.fcast$fitted)
    accuracy.ets.fcast <- accuracy(fit.ets.fcast)[1:6]
    
    for (i in 1:nsim) {
      
      sim.ets.fcast <- stats::simulate(fit.ets.fcast, nsim = h)
      matrix.ets.fcast <- rbind(matrix.ets.fcast, sim.ets.fcast)
      
    }
    
    row.names(matrix.ets.fcast) <- c(1:nsim)
    
    return(list(matrix.ets.fcast, fitted.ets.fcast, accuracy.ets.fcast))
    
  }
  
  mean.fcast <- function(ts, h, nsim) {
    
    matrix.mean.fcast <- matrix(0:0, nrow=0,ncol=h)
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
    
    row.names(matrix.mean.fcast) <- c(1:nsim)
    
    return(list(matrix.mean.fcast, fitted.mean.fcast))
    
  }
  
  naive.fcast <- function(ts, h, nsim) {
    
    matrix.naive.fcast <- matrix(0:0, nrow=0,ncol=h)
    fit.naive.fcast <- forecast::Arima(ts, order = c(0,1,0), seasonal = c(0,0,0), include.constant = FALSE, include.drift = FALSE, include.mean = FALSE)
    fitted.naive.fcast <- as.numeric(fit.naive.fcast$fitted)
    accuracy.naive.fcast <- accuracy(fit.naive.fcast)[1:6]
    
    for (i in 1:nsim) {
      
      sim.naive.fcast <- stats::simulate(fit.naive.fcast, nsim = h)
      matrix.naive.fcast <- rbind(matrix.naive.fcast, sim.naive.fcast)
      
    }
    
    row.names(matrix.naive.fcast) <- c(1:nsim)
    
    return(list(matrix.naive.fcast, fitted.naive.fcast, accuracy.naive.fcast))
    
  }
  
  snaive.fcast <- function(ts, h, nsim) {
    
    matrix.snaive.fcast <- matrix(0:0, nrow=0,ncol=h)
    fit.snaive.fcast <- forecast::Arima(ts, order = c(0,0,0), seasonal = c(0,1,0), include.constant = FALSE, include.drift = FALSE, include.mean = FALSE)
    fitted.snaive.fcast <- as.numeric( fit.snaive.fcast$fitted )
    accuracy.snaive.fcast <- accuracy(fit.snaive.fcast)[1:6]
    
    for (i in 1:nsim) {
      
      sim.snaive.fcast <- stats::simulate(fit.snaive.fcast, nsim = h)
      matrix.snaive.fcast <- rbind(matrix.snaive.fcast, sim.snaive.fcast)
      
    }
    
    row.names(matrix.snaive.fcast) <- c(1:nsim)
    
    return(list(matrix.snaive.fcast, fitted.snaive.fcast, accuracy.snaive.fcast))
    
  }
  
  rwd.fcast <- function(ts, h, nsim) {
    
    matrix.rwd.fcast <- matrix(0:0, nrow=0,ncol=h)
    fit.rwd.fcast <- forecast::Arima(ts, order = c(0,1,0), seasonal = c(0,0,0), include.constant = TRUE)
    fitted.rwd.fcast <- as.numeric(fit.rwd.fcast$fitted)
    accuracy.rwd.fcast <- accuracy(fit.rwd.fcast)[1:6]
    
    for (i in 1:nsim) {
      
      sim.rwd.fcast <- stats::simulate(fit.rwd.fcast, nsim = h)
      matrix.rwd.fcast <- rbind(matrix.rwd.fcast, sim.rwd.fcast)
      
    }
    
    row.names(matrix.rwd.fcast) <- c(1:nsim)
    
    return(list(matrix.rwd.fcast, fitted.rwd.fcast, accuracy.rwd.fcast))
    
  }
  
  nnetar.fcast <- function(ts, h, nsim) {
    
    matrix.nnetar.fcast <- matrix(0:0, nrow=0,ncol=h)
    fit.nnetar.fcast <- forecast::nnetar(ts)
    fitted.nnetar.fcast <- as.numeric(fit.nnetar.fcast$fitted)
    accuracy.nnetar.fcast <- accuracy(fit.nnetar.fcast)[1:6]
    
    for (i in 1:nsim) {
      
      sim.nnetar.fcast <- stats::simulate(fit.nnetar.fcast, nsim = h)
      matrix.nnetar.fcast <- rbind(matrix.nnetar.fcast, sim.nnetar.fcast)
      
    }
    
    row.names(matrix.nnetar.fcast) <- c(1:nsim)
    
    return(list(matrix.nnetar.fcast, fitted.nnetar.fcast, accuracy.nnetar.fcast))
    
  }
  
  matrix.results.fcast <- matrix(0:0, nrow=nsim0,ncol=h0)
  
  fitted.results.fcast <- matrix(0:0, nrow=length(ts0),ncol=1)
  
  matrix.tscv.actual <- matrix(0:0, nrow = 0, ncol = h0)
  for (i in 1:preds0) {
    l <- length(ts0)
    matrix.tscv.actual <- rbind(
      matrix.tscv.actual,
      ts0[(l-preds0+i-h0+1):(l-preds0+i)]
    )
  }
  
  matrix.tscv.results <- matrix(0:0, nrow=0, ncol=1)
  colnames(matrix.tscv.results) <- "Mean Squared Diff Score"
  
  accuracy.results.fcast <- matrix(0:0, nrow=0,ncol=6)
  colnames(accuracy.results.fcast) <- c("ME","RMSE","MAE","MPE","MAPE","MASE")
  
  if ("a" %in% methods0) {
    
    arima.results <- arima.fcast(ts = ts0, h = h0, nsim = nsim0)
    matrix.results.fcast <- matrix.results.fcast + arima.results[[1]]
    fitted.results.fcast <- fitted.results.fcast + arima.results[[2]]
    arima.accuracy <- arima.results[[3]]
    accuracy.results.fcast <- rbind(accuracy.results.fcast, 
                                    arima.accuracy)
    
    if (tscv0 == TRUE) {
      
      matrix.tscv.arima <- matrix(0:0, nrow=0, ncol=h0)
      
      for (i in 1:preds0) {
        
        increment.arima <- preds0-i+h0
        window_end.arima <- tsp(ts0)[2] - increment.arima/(tsp(ts0)[3]) 
        tscv_result.arima <- window(ts0, end=window_end.arima, freq = tsp(3)[4]) %>%
          arima.fcast(h = h0, nsim = nsim0) %>%
          '[['(1) %>%
          apply(MARGIN = 2, FUN = mean) 
        matrix.tscv.arima <- rbind(matrix.tscv.arima, tscv_result.arima)
        
      }
      
      results.tscv.arima <- (matrix.tscv.arima - matrix.tscv.actual)^2 %>%
        colMeans() %>%
        mean()
      
      matrix.tscv.results <- rbind(matrix.tscv.results,
                                   arima = results.tscv.arima)
      
    }
    
  }
  
  if ("e" %in% methods0) {
    
    ets.results <- ets.fcast(ts = ts0, h = h0, nsim = nsim0)
    matrix.results.fcast <- matrix.results.fcast + ets.results[[1]]
    fitted.results.fcast <- fitted.results.fcast + ets.results[[2]]
    ets.accuracy <- ets.results[[3]]
    accuracy.results.fcast <- rbind(accuracy.results.fcast, 
                                    ets.accuracy)
    
    if (tscv0 == TRUE) {
      
      matrix.tscv.ets <- matrix(0:0, nrow=0, ncol=h0)
      
      for (i in 1:preds0) {
        
        increment.ets <- preds0-i+h0
        window_end.ets <- tsp(ts0)[2] - increment.ets/(tsp(ts0)[3]) 
        tscv_result.ets <- window(ts0, end=window_end.ets, freq = tsp(3)[4]) %>%
          ets.fcast(h = h0, nsim = nsim0) %>%
          '[['(1) %>%
          apply(MARGIN = 2, FUN = mean) 
        matrix.tscv.ets <- rbind(matrix.tscv.ets, tscv_result.ets)
        
      }
      
      results.tscv.ets <- (matrix.tscv.ets - matrix.tscv.actual)^2 %>%
        colMeans() %>%
        mean()
      
      matrix.tscv.results <- rbind(matrix.tscv.results,
                                   ets = results.tscv.ets)
      
    }
    
  }
  
  if ("m" %in% methods0) {
    
    mean.results <- mean.fcast(ts = ts0, h = h0, nsim = nsim0)
    matrix.results.fcast <- matrix.results.fcast + mean.results[[1]]
    fitted.results.fcast <- fitted.results.fcast + mean.results[[2]]
    
    if (tscv0 == TRUE) {
      
      matrix.tscv.mean <- matrix(0:0, nrow=0, ncol=h0)
      
      for (i in 1:preds0) {
        
        increment.mean <- preds0-i+h0
        window_end.mean <- tsp(ts0)[2] - increment.mean/(tsp(ts0)[3]) 
        tscv_result.mean <- window(ts0, end=window_end.mean, freq = tsp(3)[4]) %>%
          mean.fcast(h = h0, nsim = nsim0) %>%
          '[['(1) %>%
          apply(MARGIN = 2, FUN = mean) 
        matrix.tscv.mean <- rbind(matrix.tscv.mean, tscv_result.mean)
        
      }
      
      results.tscv.mean <- (matrix.tscv.mean - matrix.tscv.actual)^2 %>%
        colMeans() %>%
        mean()
      
      matrix.tscv.results <- rbind(matrix.tscv.results,
                                   mean = results.tscv.mean)
      
    }
    
  }
  
  if ("n" %in% methods0) {
    
    naive.results <- naive.fcast(ts = ts0, h = h0, nsim = nsim0)
    matrix.results.fcast <- matrix.results.fcast + naive.results[[1]]
    fitted.results.fcast <- fitted.results.fcast + naive.results[[2]]
    naive.accuracy <- naive.results[[3]]
    accuracy.results.fcast <- rbind(accuracy.results.fcast, 
                                    naive.accuracy)
    
    if (tscv0 == TRUE) {
      
      matrix.tscv.naive <- matrix(0:0, nrow=0, ncol=h0)
      
      for (i in 1:preds0) {
        
        increment.naive <- preds0-i+h0
        window_end.naive <- tsp(ts0)[2] - increment.naive/(tsp(ts0)[3]) 
        tscv_result.naive <- window(ts0, end=window_end.naive, freq = tsp(3)[4]) %>%
          naive.fcast(h = h0, nsim = nsim0) %>%
          '[['(1) %>%
          apply(MARGIN = 2, FUN = mean) 
        matrix.tscv.naive <- rbind(matrix.tscv.naive, tscv_result.naive)
        
      }
      
      results.tscv.naive <- (matrix.tscv.naive - matrix.tscv.actual)^2 %>%
        colMeans() %>%
        mean()
      
      matrix.tscv.results <- rbind(matrix.tscv.results,
                                   naive = results.tscv.naive)
      
    }
    
  }
  
  if ("s" %in% methods0) {
    
    snaive.results <- snaive.fcast(ts = ts0, h = h0, nsim = nsim0)
    matrix.results.fcast <- matrix.results.fcast + snaive.results[[1]]
    fitted.results.fcast <- fitted.results.fcast + snaive.results[[2]]
    snaive.accuracy <- snaive.results[[3]]
    accuracy.results.fcast <- rbind(accuracy.results.fcast, 
                                    snaive.accuracy)
    
    if (tscv0 == TRUE) {
      
      matrix.tscv.snaive <- matrix(0:0, nrow=0, ncol=h0)
      
      for (i in 1:preds0) {
        
        increment.snaive <- preds0-i+h0
        window_end.snaive <- tsp(ts0)[2] - increment.snaive/(tsp(ts0)[3]) 
        tscv_result.snaive <- window(ts0, end=window_end.snaive, freq = tsp(3)[4]) %>%
          snaive.fcast(h = h0, nsim = nsim0) %>%
          '[['(1) %>%
          apply(MARGIN = 2, FUN = mean) 
        matrix.tscv.snaive <- rbind(matrix.tscv.snaive, tscv_result.snaive)
        
      }
      
      results.tscv.snaive <- (matrix.tscv.snaive - matrix.tscv.actual)^2 %>%
        colMeans() %>%
        mean()
      
      matrix.tscv.results <- rbind(matrix.tscv.results,
                                   snaive = results.tscv.snaive)
      
    }
    
  }
  
  if ("rw" %in% methods0) {
    
    rwd.results <- rwd.fcast(ts = ts0, h = h0, nsim = nsim0)
    matrix.results.fcast <- matrix.results.fcast + rwd.results[[1]]
    fitted.results.fcast <- fitted.results.fcast + rwd.results[[2]]
    rwd.accuracy <- rwd.results[[3]]
    accuracy.results.fcast <- rbind(accuracy.results.fcast, 
                                    rwd.accuracy)
    
    if (tscv0 == TRUE) {
      
      matrix.tscv.rwd <- matrix(0:0, nrow=0, ncol=h0)
      
      for (i in 1:preds0) {
        
        increment.rwd <- preds0-i+h0
        window_end.rwd <- tsp(ts0)[2] - increment.rwd/(tsp(ts0)[3]) 
        tscv_result.rwd <- window(ts0, end=window_end.rwd, freq = tsp(3)[4]) %>%
          rwd.fcast(h = h0, nsim = nsim0) %>%
          '[['(1) %>%
          apply(MARGIN = 2, FUN = mean) 
        matrix.tscv.rwd <- rbind(matrix.tscv.rwd, tscv_result.rwd)
        
      }
      
      results.tscv.rwf <- (matrix.tscv.rwd - matrix.tscv.actual)^2 %>%
        colMeans() %>%
        mean()
      
      matrix.tscv.results <- rbind(matrix.tscv.results,
                                   rwd = results.tscv.rwd)
      
    }
    
  }
  
  if ("nn" %in% methods0) {
    
    nnetar.results <- nnetar.fcast(ts = ts0, h = h0, nsim = nsim0)
    matrix.results.fcast <- matrix.results.fcast + nnetar.results[[1]]
    fitted.results.fcast <- fitted.results.fcast + nnetar.results[[2]]
    nnetar.accuracy <- nnetar.results[[3]]
    accuracy.results.fcast <- rbind(accuracy.results.fcast, 
                                    nnetar.accuracy)
    
    if (tscv0 == TRUE) {
      
      matrix.tscv.nnetar <- matrix(0:0, nrow=0, ncol=h0)
      
      for (i in 1:preds0) {
        
        increment.nnetar <- preds0-i+h0
        window_end.nnetar <- tsp(ts0)[2] - increment.nnetar/(tsp(ts0)[3]) 
        tscv_result.nnetar <- window(ts0, end=window_end.nnetar, freq = tsp(3)[4]) %>%
          nnetar.fcast(h = h0, nsim = nsim0) %>%
          '[['(1) %>%
          apply(MARGIN = 2, FUN = mean) 
        matrix.tscv.nnetar <- rbind(matrix.tscv.nnetar, tscv_result.nnetar)
        
      }
      
      results.tscv.nnetar <- (matrix.tscv.nnetar - matrix.tscv.actual)^2 %>%
        colMeans() %>%
        mean()
      
      matrix.tscv.results <- rbind(matrix.tscv.results,
                                   nnetar = results.tscv.nnetar)
      
    }
    
  }
  
  matrix.results.fcast <- (1/length(methods0)) * matrix.results.fcast
  fitted.results.fcast <- (1/length(methods0)) * fitted.results.fcast
  
  ts.freq <- tsp(ts0)[3]
  ts.end_date <- tsp(ts0)[2]
  ts.start_date <- tsp(ts0)[1]
  
  upper.quantile <- apply(matrix.results.fcast, 
                          MARGIN = 2, 
                          FUN = quantile, 
                          probs = 1 - (1-conf0)/2)
  lower.quantile <- apply(matrix.results.fcast, 
                          MARGIN = 2, 
                          FUN = quantile, 
                          probs = (1-conf0)/2)
  middle.quantile <- apply(matrix.results.fcast, 
                  MARGIN = 2, 
                  FUN = mean)
  
  upper.fcast <- ts(data = upper.quantile,
                    start = ts.end_date + (1/ts.freq),
                    freq = ts.freq)
  lower.fcast <- ts(data = lower.quantile,
                    start = ts.end_date + (1/ts.freq),
                    freq = ts.freq)
  middle.fcast <- ts(data = middle.quantile,
                     start = ts.end_date + (1/ts.freq),
                     freq = ts.freq)
  fitted.fcast <- ts(data = fitted.results.fcast,
                     start = ts.start_date,
                     freq = ts.freq)
  
  if (boxcox = TRUE) {
  
    object.fcast <- list(x = forecast::InvBoxCox(ts0, lambda=lambda0),
                       mean = forecast::InvBoxCox(middle.fcast, lambda=lambda0),
                       lower = forecast::InvBoxCox(lower.fcast, lambda=lambda0),
                       upper = forecast::InvBoxCox(upper.fcast, lambda=lambda0),
                       fitted = forecast::InvBoxCox(fitted.fcast, lambda=lambda0),
                       level = (conf0*100)) %>% structure(class="forecast")
  
  } else {
    
    object.fcast <- list(x = ts0,
                         mean = middle.fcast,
                         lower = lower.fcast,
                         upper = upper.fcast,
                         fitted = fitted.fcast, lambda=lambda0,
                         level = (conf0*100)) %>% structure(class="forecast")
    
  }
  
#  structure.fcast <- structure(list.fcast, class="forecast")
  
  if (tscv0 == TRUE) {
    
    return(list(object.fcast, matrix.tscv.results, accuracy.results.fcast))
    
  } else {
    
    return(list(object.fcast, accuracy.results.fcast))
    
  }
  
}

y <- get_fred(
  frequency = "q",
  fred_api_key = "0962c1d1b6316b7c4b6c68ced8211ac9",
  fred_series_id = "UNRATE",
  units = NULL) %>%
  fcast(methods = c("m","n","s","a","e","nn"),
        tscv = TRUE,
        h = 6)
