# FILE: /EmpStudyProject/EmpStudyProject/MethodCode/PELTtrendARp.R
# This file contains the implementation of the PELT (Pruned Exact Linear Time) algorithm for detecting change points in time series data.

PELT.trendARp <- function(Y, p, pen, minseglen) {
  # Implementation of the PELT algorithm for detecting change points
  # Y: time series data
  # p: order of the autoregressive model
  # pen: penalty term for the number of change points
  # minseglen: minimum segment length
  
  n <- length(Y)
  cost <- numeric(n + 1)
  change_points <- list()
  
  # Initialize cost for the first segment
  cost[1] <- 0
  
  for (t in 2:(n + 1)) {
    cost[t] <- Inf
    for (s in max(1, t - minseglen): (t - 1)) {
      segment <- Y[s:(t - 1)]
      segment_fit <- trendARpsegfit_arimax(segment, p)
      segment_cost <- sum((segment - segment_fit$trend_fit)^2) + pen
      if (cost[s] + segment_cost < cost[t]) {
        cost[t] <- cost[s] + segment_cost
        change_points[[t]] <- s - 1
      }
    }
  }
  
  # Extract change points
  cp <- unlist(change_points)
  return(cp)
}

trendARpsegfit_arimax <- function(y, p) {
  # Fit an ARIMA model to the segment
  fit <- arima(y, order = c(p, 0, 0), include.mean = TRUE)
  trend_fit <- fitted(fit)
  
  return(list(trend_fit = trend_fit))
}