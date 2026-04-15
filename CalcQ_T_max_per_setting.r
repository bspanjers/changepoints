##########################
###### Vantage Simulations
##########################

library(parallel)
library(doParallel)
library(foreach)
source('./MethodCode/PELTtrendARpJOIN.R')

# Load data
load("./data/temperature_anomalies.RData")
ANOM <- Tanom_annual_df[,]
source('./MethodCode/PELTtrendARpJOIN.R')

seed <- 12345
set.seed(seed) # for reproducibility


########### Continuous model AR(1) 
trendarjoin=list()
mresiduals=list()
fits=list()
dates=list()

setting_to_col <- c(
  "NASA" = 2,
  "HadCRUT" = 4,
  "NOAA" = 5,
  "Berkeley" = 6
)


# Define the simulation function
simulateone <- function(index = 1, n, interc, slope, sd, phi) {
  tryCatch({
    seg <- 1:n
    y <- interc + slope * seg + arima.sim(n = n, sd = sd, model = list(ar = phi))
    min_len <- max(round(0.1 * n), 2)
    max_len <- min(n - round(0.1 * n), n - 2)
    seglen <- max_len:min_len
    stats <- numeric(length(seglen))
    
    for (i in seq_along(seglen)) {
      seg2 <- (n - seglen[i] + 1):n
      seg1 <- 1:(n - seglen[i])
      X <- cbind(
        rep(1, n),
        c(seg1, rep(seg1[length(seg1)], seglen[i])),
        c(rep(0, (n - seglen[i])), seg2 - seg2[1] + 1)
      )
      vec <- c(0, 0, -1, 1)

      armafit <- tryCatch({
        arima(y, xreg = X, order = c(1, 0, 0), include.mean = FALSE)
      }, error = function(e) {
        return(NULL)
      })

      if (is.null(armafit)) {
        stats[i] <- NA
        next
      }

      sddiff <- sqrt(t(vec) %*% armafit$var.coef %*% vec)
      stats[i] <- (armafit$coef[3] - armafit$coef[4]) / sddiff
    }

    return(max(abs(stats), na.rm = TRUE))

  }, error = function(e) {
    return(NA_real_)
  })
}

gettmax <- function(n, interc, slope, sd, phi, n_sim) {
  set.seed(n)
  hold <- sapply(1:n_sim, simulateone, n = n, interc = interc, slope = slope, sd = sd, phi = phi)
  return(quantile(hold, 0.95))
}

gettmax_chunked <- function(n, interc, slope, sd, phi, n_sim, chunk_size = 1000) {
  n_chunks <- ceiling(n_sim / chunk_size)
  hold <- numeric(n_sim)
  pos <- 1

  for (k in seq_len(n_chunks)) {
    m <- min(chunk_size, n_sim - pos + 1)
    hold[pos:(pos + m - 1)] <- replicate(
      m,
      simulateone(n = n, interc = interc, slope = slope, sd = sd, phi = phi)
    )
    pos <- pos + m
  }

  quantile(hold, 0.95)
}


ns <- 33:78
n_sim <- 100000
n_cores <-46 # detectCores() - 1
settings_to_run <- c("NOAA")

for (setting in settings_to_run) {
  col_index <- unname(setting_to_col[setting])
  if (is.na(col_index)) stop(sprintf("Invalid setting specified: %s", setting))

  data <- Tanom_annual_df[, c(1, col_index)][!is.na(Tanom_annual_df[, col_index]), ]
  y_full <- data[, 2]
  years <- data[, 1]
  n <- nrow(data)

  # Model fit
  itrendarjoin <- PELT.trendARpJOIN(y_full, p = 1, pen = 4 * log(n), minseglen = 10)
  fittrend <- fit.trendARpJOIN_with_se(
    y_full, itrendarjoin, p = 1, dates = years, plot = F, add.ar = F, fit = T,
    title = names(Tanom_annual_df[col_index]), pred = F
  )
  coeffs_matrix <- fittrend$coeffs

  y_after_1970 <- y_full[121:length(y_full)] # beaulieu fix change point at 1970 (index 121)
  times <- 1:length(y_after_1970)
  shortfit <- arima(y_after_1970, order = c(1, 0, 0), xreg = times)

  # Extract point estimates
  phi_hat <- shortfit$coef["ar1"]
  intercept_hat <- shortfit$coef["intercept"]
  slope_hat <- shortfit$coef["times"]
  sigma2_hat <- shortfit$sigma2
  sd_hat <- sqrt(sigma2_hat)

  shortfit_params <- list(
    phi = as.numeric(phi_hat),
    intercept = as.numeric(intercept_hat),
    slope = as.numeric(slope_hat),
    sd = sd_hat
  )

  cat(sprintf("\n=== Setting %s - Fitted ARIMA(1,0,0) with Trend (1970-2023) ===\n", setting))
  cat(sprintf("Point Estimates:\n"))
  cat(sprintf("  Intercept = %.6f \n", intercept_hat))
  cat(sprintf("  Slope = %.6f \n", slope_hat))
  cat(sprintf("  AR1 = %.6f \n", phi_hat))
  cat(sprintf("  SD = %.6f \n", sd_hat))
  cat(sprintf("Using %d cores\n", n_cores))

  # Register parallel backend
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
  cl <- makeCluster(n_cores, outfile = "")
  registerDoParallel(cl)
  clusterSetRNGStream(cl, iseed = seed + col_index)

  # Export required functions and data to workers
  clusterExport(cl, c("simulateone", "gettmax", "gettmax_chunked", "shortfit_params", "ns"))

  # Create a task grid per setting and skip completed tasks
  tasks <- expand.grid(n = ns)
  tasks$save_file <- sprintf("./Results/Tmaxquant%s_n%d.Rdata", setting, tasks$n)
  tasks <- tasks[!file.exists(tasks$save_file), ]

  if (nrow(tasks) > 0) {
    results <- foreach(task_idx = 1:nrow(tasks), .packages = c("stats"), .errorhandling = "pass") %dopar% {
      task <- tasks[task_idx, ]
      n_val <- task$n
      save_file <- sprintf("./Results/Tmaxquant%s_n%d.Rdata", setting, n_val)

      if (file.exists(save_file)) {
        cat(sprintf("Setting %s, n=%d already completed. Skipping...\n", setting, n_val))
        return("skipped")
      }

      QN <- gettmax_chunked(
        n = n_val,
        interc = shortfit_params$intercept,
        slope = shortfit_params$slope,
        sd = shortfit_params$sd,
        phi = shortfit_params$phi,
        n_sim = n_sim,
        chunk_size = 1000
      )
      save(QN, file = save_file)
      cat(sprintf("Saved result for setting %s, n=%d: %.6f\n", setting, n_val, QN))
      "completed"
    }
  } else {
    cat(sprintf("All n values already completed for setting %s.\n", setting))
  }

  stopCluster(cl)

  # Aggregate all n results into one per-setting file
  setting_files <- list.files(
    "./Results",
    pattern = sprintf("Tmaxquant%s_n[0-9]+\\.Rdata", setting),
    full.names = TRUE
  )

  if (length(setting_files) > 0) {
    setting_results_list <- lapply(setting_files, function(file) {
      env <- new.env()
      load(file, envir = env)
      n_val <- as.integer(sub(".*_n([0-9]+)\\.Rdata$", "\\1", basename(file)))
      data.frame(ns = n_val, QNs = env$QN)
    })

    setting_results <- do.call(rbind, setting_results_list)
    setting_results <- setting_results[order(setting_results$ns), ]
    save(setting_results, file = sprintf("./Results/Tmaxquant%s.Rdata", setting))
    cat(sprintf("Saved aggregated results for %s: ./Results/Tmaxquant%s.Rdata\n", setting, setting))
  } else {
    cat(sprintf("No per-n result files found to aggregate for setting %s\n", setting))
  }
}

cat("Task completed for all settings!\n")
