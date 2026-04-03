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


########### Continuous model AR(1)
trendarjoin=list()
mresiduals=list()
fits=list()
dates=list()

data=Tanom_annual_df[,c(1,4)][!is.na(Tanom_annual_df[,4]),]
y_full = data[,2]
years = data[,1]
n=nrow(data)

#Model fit
itrendarjoin=PELT.trendARpJOIN(y_full, p=1,pen=4*log(n),minseglen=10)
fittrend = fit.trendARpJOIN_with_se(y_full, itrendarjoin,p=1,dates=years,plot=F,add.ar=F,fit=T,
                            title=names(Tanom_annual_df[4]),pred=F)# get fit without AR - to visualize trend segments
fits = fittrend$fit
dates = fittrend$dates
fittrendAR = fit.trendARpJOIN_with_se(y_full,itrendarjoin,p=1,dates=years,plot=F,add.ar=T,fit=T,
                            title=names(Tanom_annual_df[4])) #get fit with AR - to compute residuals below
yeartrendarjoin=years[itrendarjoin]  # put cpts in terms of year
cat(paste(names(Tanom_annual_df),':TrendAR1join \n'))


# Extract coefficients from the fitted trend model with 2 segments
coeffs_matrix <- fittrend$coeffs
se_matrix <- fittrend$se_coeffs

cat("\nFitted coefficients (2 segments):\n")
print(coeffs_matrix)

# Get residuals from AR fit to estimate SD
residuals_ar <- y_full - fittrendAR$fit
residuals_ar <- residuals_ar[!is.na(residuals_ar)]
sd_residuals <- sd(residuals_ar)

# Extract coefficients for 2-segment model
carry_level <- coeffs_matrix[1, 1]      # CarryLevel (intercept)
beta_seg1 <- coeffs_matrix[1, 2]        # Beta segment 1
beta_seg2 <- coeffs_matrix[2, 2]        # Beta segment 2
ar1_seg1 <- coeffs_matrix[1, 3]         # AR1 segment 1
ar1_seg2 <- coeffs_matrix[2, 3]         # AR1 segment 2

# Extract standard errors
se_carry_level <- se_matrix[1, 1]
se_beta_seg1 <- se_matrix[1, 2]
se_beta_seg2 <- se_matrix[2, 2]
se_ar1_seg1 <- se_matrix[1, 3]
se_ar1_seg2 <- se_matrix[2, 3]
se_sd <- sd_residuals * 0.1  # Approximate SE for SD

# 95% CI bounds (z = 1.96)
z <- 1.96

# Define 5 settings within 95% confidence region
# Vary the slopes and AR1s while keeping intercept and SD at center
settings_full <- data.frame(
  setting = 1:5,
  intercept = rep(carry_level, 5),
  
  slope_seg1 = c(
    beta_seg1,                           # 1: center
    beta_seg1 + 1 * z * se_beta_seg1,  # 2: upper-middle
    beta_seg1 - 0.5 * z * se_beta_seg1,  # 3: lower-middle
    beta_seg1 + 0.5 * z * se_beta_seg1, # 4: upper
    beta_seg1 - 1 * z * se_beta_seg1  # 5: lower
  ),
  
  slope_seg2 = c(
    beta_seg2,                           # 1: center
    beta_seg2 + 1 * z * se_beta_seg2,  # 2: upper-middle
    beta_seg2 - 1 * z * se_beta_seg2,  # 3: lower-middle
    beta_seg2 - 0.5 * z * se_beta_seg2, # 4: lower
    beta_seg2 + .5 * z * se_beta_seg2  # 5: upper
  ),
  
  phi_seg1 = c(
    ar1_seg1,                           # 1: center
    ar1_seg1 + 1 * z * se_ar1_seg1,   # 2: upper-middle
    ar1_seg1 - 0.5 * z * se_ar1_seg1,   # 3: lower-middle
    ar1_seg1 + 0.5 * z * se_ar1_seg1,  # 4: slight upper
    ar1_seg1 - 1 * z * se_ar1_seg1   # 5: lower
  ),
  
  phi_seg2 = c(
    ar1_seg2,                           # 1: center
    ar1_seg2 - 1 * z * se_ar1_seg2,   # 2: lower-middle
    ar1_seg2 + 0.5 * z * se_ar1_seg2,   # 3: upper-middle
    ar1_seg2 + 1 * z * se_ar1_seg2,  # 4: slight upper
    ar1_seg2 - 0.5 * z * se_ar1_seg2   # 5: slight lower
  ),
  
  sd = rep(sd_residuals, 5)
)

# Ensure phi values are within valid range (-1, 1)
settings_full$phi_seg1 <- pmin(pmax(settings_full$phi_seg1, -0.99), 0.99)
settings_full$phi_seg2 <- pmin(pmax(settings_full$phi_seg2, -0.99), 0.99)

# Ensure sd is positive (lower bound)
settings_full$sd <- pmax(settings_full$sd, 0.01)

# Print settings information
cat("\n=== Settings from 2-Segment Fitted Model ===\n")
cat(sprintf("Point Estimates:\n"))
cat(sprintf("  CarryLevel (intercept) = %.6f (SE: %.6f)\n", carry_level, se_carry_level))
cat(sprintf("  Beta Segment 1 = %.6f (SE: %.6f)\n", beta_seg1, se_beta_seg1))
cat(sprintf("  Beta Segment 2 = %.6f (SE: %.6f)\n", beta_seg2, se_beta_seg2))
cat(sprintf("  AR1 Segment 1 = %.6f (SE: %.6f)\n", ar1_seg1, se_ar1_seg1))
cat(sprintf("  AR1 Segment 2 = %.6f (SE: %.6f)\n", ar1_seg2, se_ar1_seg2))
cat(sprintf("  SD (residuals) = %.6f (SE: %.6f)\n\n", sd_residuals, se_sd))

cat("95% Confidence Intervals:\n")
cat(sprintf("  Beta Segment 1: [%.6f, %.6f]\n", 
            beta_seg1 - z*se_beta_seg1, beta_seg1 + z*se_beta_seg1))
cat(sprintf("  Beta Segment 2: [%.6f, %.6f]\n", 
            beta_seg2 - z*se_beta_seg2, beta_seg2 + z*se_beta_seg2))
cat(sprintf("  AR1 Segment 1: [%.6f, %.6f]\n", 
            ar1_seg1 - z*se_ar1_seg1, ar1_seg1 + z*se_ar1_seg1))
cat(sprintf("  AR1 Segment 2: [%.6f, %.6f]\n\n", 
            ar1_seg2 - z*se_ar1_seg2, ar1_seg2 + z*se_ar1_seg2))

cat("5 Settings (all parameters within 95% CI):\n")
print(settings_full)

# Validation
cat("\nValidation - All settings within bounds:\n")
for (i in 1:nrow(settings_full)) {
  in_bounds <- (
    settings_full$slope_seg1[i] >= beta_seg1 - z*se_beta_seg1 && 
    settings_full$slope_seg1[i] <= beta_seg1 + z*se_beta_seg1 &&
    settings_full$slope_seg2[i] >= beta_seg2 - z*se_beta_seg2 && 
    settings_full$slope_seg2[i] <= beta_seg2 + z*se_beta_seg2 &&
    settings_full$phi_seg1[i] >= ar1_seg1 - z*se_ar1_seg1 && 
    settings_full$phi_seg1[i] <= ar1_seg1 + z*se_ar1_seg1 &&
    settings_full$phi_seg2[i] >= ar1_seg2 - z*se_ar1_seg2 && 
    settings_full$phi_seg2[i] <= ar1_seg2 + z*se_ar1_seg2
  )
  cat(sprintf("Setting %d within bounds: %s\n", i, ifelse(in_bounds, "YES", "NO")))
}

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

# Create fittrend objects for each setting
fittrend_list <- list()

for (setting_id in 1:nrow(settings_full)) {
  setting_row <- settings_full[setting_id, ]
  
  # Create coefficients matrix matching the structure of fittrend$coeffs
  # Row 1: CarryLevel, Beta_seg1, AR1_seg1
  # Row 2: NA, Beta_seg2, AR1_seg2
  coeffs_setting <- matrix(NA_real_, nrow = 2, ncol = 3)
  colnames(coeffs_setting) <- c("CarryLevel", "Beta", "AR1")
  
  coeffs_setting[1, 1] <- setting_row$intercept
  coeffs_setting[1, 2] <- setting_row$slope_seg1
  coeffs_setting[1, 3] <- setting_row$phi_seg1
  coeffs_setting[2, 2] <- setting_row$slope_seg2
  coeffs_setting[2, 3] <- setting_row$phi_seg2
  
  # Create a fittrend-like object for this setting
  fittrend_list[[setting_id]] <- list(
    coeffs = coeffs_setting,
    se_coeffs = se_matrix,  # Use the SEs from original fit
    fit = fittrend$fit,      # Reuse the original fit
    dates = fittrend$dates
  )
  
  # Assign to global environment with names fittrend1, fittrend2, etc.
  assign(sprintf("fittrend%d", setting_id), fittrend_list[[setting_id]], envir = .GlobalEnv)
}

# Now you can access them as:
# fittrend1$coeffs, fittrend2$coeffs, etc.

cat("\nCreated fittrend objects for all 5 settings:\n")
for (i in 1:5) {
  cat(sprintf("\nfittrend%d$coeffs:\n", i))
  print(get(sprintf("fittrend%d", i))$coeffs)
}

# for each setting simulate one path via simulate_trendARpJOIN (not for the first setting as we have the original data). 
# we use numeric(n) just to get the size of the data, the actual values will be generated based on the coefficients in fittrend objects
seed <- 12345
ysetting1 = y_full
set.seed(seed) # for reproducibility
ysetting2 = simulate_trendARpJOIN(numeric(n), fittrend2, 121)
ysetting3 = simulate_trendARpJOIN(numeric(n), fittrend3, 121)
ysetting4 = simulate_trendARpJOIN(numeric(n), fittrend4, 121)
ysetting5 = simulate_trendARpJOIN(numeric(n), fittrend5, 121)

y1_after_1970 = ysetting1[121:length(ysetting1)]
y2_after_1970 = ysetting2[121:length(ysetting2)]
y3_after_1970 = ysetting3[121:length(ysetting3)]
y4_after_1970 = ysetting4[121:length(ysetting4)]
y5_after_1970 = ysetting5[121:length(ysetting5)]

times <- 1:length(y1_after_1970) # Create a time index for the ARIMA model (starting from 1 for the post-1970 data)

shortfit1 <- arima(y1_after_1970, order = c(1, 0, 0), xreg = times)
shortfit2 <- arima(y2_after_1970, order = c(1, 0, 0), xreg = times)
shortfit3 <- arima(y3_after_1970, order = c(1, 0, 0), xreg = times)
shortfit4 <- arima(y4_after_1970, order = c(1, 0, 0), xreg = times)
shortfit5 <- arima(y5_after_1970, order = c(1, 0, 0), xreg = times)

# Extract parameters and SEs for all 5 settings
shortfits <- list(shortfit1, shortfit2, shortfit3, shortfit4, shortfit5)
shortfit_params <- list()

for (i in 1:5) {
  shortfit <- shortfits[[i]]
  
  # Extract point estimates
  phi_hat <- shortfit$coef["ar1"]
  intercept_hat <- shortfit$coef["intercept"]
  slope_hat <- shortfit$coef["times"]
  sigma2_hat <- shortfit$sigma2

  # Point estimate for SD
  sd_hat <- sqrt(sigma2_hat)
  
  # Store results in a list
  shortfit_params[[i]] <- list(
    setting = i,
    phi = as.numeric(phi_hat),
    intercept = as.numeric(intercept_hat),
    slope = as.numeric(slope_hat),
    sd = sd_hat
  )
  
  # Print summary for this setting
  cat(sprintf("\n=== Setting %d - Fitted ARIMA(1,0,0) with Trend (1970-2023) ===\n", i))
  cat(sprintf("Point Estimates:\n"))
  cat(sprintf("  Intercept = %.6f \n", intercept_hat))
  cat(sprintf("  Slope = %.6f \n", slope_hat))
  cat(sprintf("  AR1 = %.6f \n", phi_hat))
  cat(sprintf("  SD = %.6f \n", sd_hat))
  
}

# Convert to data frame for easy viewing
settings <- do.call(rbind, lapply(shortfit_params, function(x) {
  data.frame(
    Setting = x$setting,
    Intercept = x$intercept,
    Slope = x$slope,
    AR1 = x$phi,
    SD = x$sd
  )
}))

print(settings)





# Parallel processing
ns <- 33:71#33:71
n_sim <- 7
n_cores <- 7#detectCores() - 1
cat(sprintf("Using %d cores\n", n_cores))

# Register parallel backend
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
cl <- makeCluster(n_cores, outfile="")
registerDoParallel(cl)

# Export required functions and data to workers
clusterExport(cl, c("simulateone", "gettmax", "settings", "ns"))

# Create a task grid: all combinations of settings and ns
tasks <- expand.grid(setting_id = 1:nrow(settings), n = 33:71)
tasks$save_file <- sprintf("./Results/TmaxquantHad_setting%d_n%d.Rdata",
                           tasks$setting_id, tasks$n)

tasks <- tasks[!file.exists(tasks$save_file), ]

# Run all tasks in parallel
results <- foreach(task_idx = 1:nrow(tasks), .packages = c("stats"), .errorhandling = "pass") %dopar% {
  task <- tasks[task_idx, ]
  setting_id <- task$setting_id
  n_val <- task$n
  setting_row <- settings[setting_id, ]
  
  # Save file for this task
  save_file <- sprintf("./Results/TmaxquantHad_setting%d_n%d.Rdata", setting_id, n_val)
  
  # Skip if already completed
  if (file.exists(save_file)) {
    cat(sprintf("Setting %d, n=%d already completed. Skipping...\n", setting_id, n_val))
    return("skipped")
  }
  
  # Run the simulation
  QN <- gettmax_chunked(
    n = n_val,
    interc = setting_row$Intercept,
    slope = setting_row$Slope,
    sd = setting_row$SD,
    phi = setting_row$AR1,
    n_sim = n_sim,
    chunk_size = 1000
  )
  # Save result to disk immediately
  save(QN, file = save_file)
  cat(sprintf("Saved result for setting %d, n=%d: %.6f\n", setting_id, n_val, QN))
  
  # Return only a simple status indicator, not the full result
  "completed"
}

# Stop cluster
stopCluster(cl)

# Aggregate results per setting
for (setting_id in 1:nrow(settings)) {
  setting_files <- list.files("./Results", 
                              pattern = sprintf("TmaxquantHad_setting%d_n[0-9]+\\.Rdata", setting_id),
                              full.names = TRUE)
  
  if (length(setting_files) > 0) {
    setting_results_list <- lapply(setting_files, function(file) {
      env <- new.env()
      load(file, envir = env)
      QN <- env$QN
      n_val <- as.numeric(gsub(".*_n([0-9]+)\\.Rdata", "\\1", file))
      data.frame(n = n_val, QN = QN)
    })
    
    setting_results <- do.call(rbind, setting_results_list)
    setting_results <- setting_results[order(setting_results$n), ]
    
    # Save aggregated results for the setting
    save(setting_results, file = sprintf("./Results/TmaxquantHad_setting%d.Rdata", setting_id))
    cat(sprintf("Saved aggregated results for setting %d\n", setting_id))
  }
}

cat("All tasks completed!\n")
