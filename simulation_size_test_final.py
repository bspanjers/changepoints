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




# ...existing code...

##############################################################################
###     BOOTSTRAP SIZE TEST FOR JOIN MODEL - ALL 5 SETTINGS (PARALLEL)     ###
##############################################################################

library(foreach)
library(doParallel)
library(parallel)

source('./MethodCode/PELTtrendARpJOIN.R')


### --- PARAMETERS --- ###
nsim       <- 20000
refyear    <- 1973
lag_thresh <- 15
tau_star   <- 124  # Index for 1973 (changepoint location)
penalty    <- 4 * log(n)

### --- SINGLE SIMULATION FUNCTION --- ###
run_one_sim <- function(i, y_setting, fittrend_setting, tau_star, years, QNs, penalty, refyear, lag_thresh) {
  
  # Simulate a path under the null (no additional changepoint after tau_star)
  simu <- simulate_trendARpJOIN(y_setting, fittrend_setting, tau_star)
  
  # Fit PELT to detect changepoints

  bp <- suppressWarnings(tryCatch(
    PELT.trendARpJOIN(simu, p = 1, pen = penalty, minseglen = 10),
    error = function(e) NULL
  ))
  
  # Classify result

if (is.null(bp) || length(bp) == 0) return(list(status = "zero", bp = NA, year = NA, Tmax = NA, thr = NA))
  if (length(bp) >= 2) return(list(status = "2ormore", bp = bp, year = years[bp[1]], Tmax = NA, thr = NA))
  
  # Single breakpoint case
  year_est <- years[bp]
  if (abs(year_est - refyear) > lag_thresh) {
    return(list(status = "far", bp = bp, year = year_est, Tmax = NA, thr = NA))
  }
  
  # Valid 1-break case: compute Tmax statistic
  n_after <- length(years) - bp + 1
  thr <- QNs[QNs$n == n_after, "QN"]
  if (length(thr) == 0) thr <- NA
  
  yafter <- simu[bp:length(simu)]
  n2     <- length(yafter)
  min_len <- max(round(0.1 * n2), 2)
  max_len <- min(n2 - round(0.1 * n2), n2 - 2)
  seglen <- max_len:min_len
  stats  <- numeric(length(seglen))
  
  for (j in seq_along(seglen)) {
    seg2 <- (n2 - seglen[j] + 1):n2
    seg1 <- 1:(n2 - seglen[j])
    
    Xreg <- cbind(
      rep(1, n2),
      c(seg1, rep(seg1[length(seg1)], seglen[j])),
      c(rep(0, (n2 - seglen[j])), seg2 - seg2[1] + 1)
    )
    
    vec <- c(0, 0, -1, 1)
    fit <- tryCatch(
      arima(yafter, xreg = Xreg, order = c(1, 0, 0), include.mean = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      stats[j] <- NA
      next
    }
    
    sdd <- sqrt(t(vec) %*% fit$var.coef %*% vec)
    stats[j] <- (fit$coef[3] - fit$coef[4]) / sdd
  }
  
  Tmax <- max(abs(stats), na.rm = TRUE)
  
  list(status = "valid", bp = bp, year = year_est, Tmax = Tmax, thr = thr)
}

### --- RUN BOOTSTRAP FOR ONE SETTING --- ###
run_bootstrap_for_setting <- function(setting_id, nsim, refyear, lag_thresh, penalty, tau_star, years, n) {
  
  cat("\n##############################################\n")
  cat(sprintf("  Running bootstrap for Setting %d\n", setting_id))
  cat("##############################################\n")
  
  # Get the fittrend object and y data for this setting
  fittrend_setting <- get(sprintf("fittrend%d", setting_id))
  y_setting <- get(sprintf("ysetting%d", setting_id))
  
  # Load quantile table for this setting
  qn_file <- sprintf("./Results/TmaxquantHad_setting%d.Rdata", setting_id)
  if (!file.exists(qn_file)) {
    stop(sprintf("Quantile file not found: %s", qn_file))
  }
  load(qn_file)  # Loads 'setting_results'
  QNs <- setting_results
  names(QNs) <- c("n", "QN")
  
  cat(sprintf("Loaded quantiles from %s\n", qn_file))
  cat(sprintf("Running %d simulations...\n", nsim))
  
  # Setup parallel backend
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export necessary objects to workers
  clusterExport(cl, c("y_setting", "fittrend_setting", "tau_star", "years", "QNs", 
                      "penalty", "refyear", "lag_thresh",
                      "simulate_trendARpJOIN", "PELT.trendARpJOIN", 
                      "fit.trendARpJOIN", "run_one_sim"),
                envir = environment())
  clusterEvalQ(cl, source('./MethodCode/PELTtrendARpJOIN.R'))
  
  # Run simulations in parallel
  results_list <- foreach(i = 1:nsim, 
                          .packages = c("stats"),
                          .errorhandling = "pass") %dopar% {
    run_one_sim(i, y_setting, fittrend_setting, tau_star, years, QNs, 
                penalty, refyear, lag_thresh)
  }
  
  stopCluster(cl)
  
  # Process results
  statuses <- sapply(results_list, function(x) {
    if (inherits(x, "error")) return("error")
    x$status
  })
  
  # Count outcomes
  n_zero     <- sum(statuses == "zero")
  n_2ormore  <- sum(statuses == "2ormore")
  n_far      <- sum(statuses == "far")
  n_valid    <- sum(statuses == "valid")
  n_error    <- sum(statuses == "error")
  
  # For valid cases, compute rejection rate
  valid_results <- results_list[statuses == "valid"]
  if (length(valid_results) > 0) {
    Tmax_vals <- sapply(valid_results, function(x) x$Tmax)
    thr_vals  <- sapply(valid_results, function(x) x$thr)
    rejections <- Tmax_vals > thr_vals
    rejection_rate <- mean(rejections, na.rm = TRUE)
    
    # Year distribution of detected changepoints
    years_detected <- sapply(valid_results, function(x) x$year)
    year_table <- table(years_detected)
  } else {
    rejection_rate <- NA
    year_table <- NULL
  }
  
  # Summary
  summary_df <- data.frame(
    Setting = setting_id,
    N_sim = nsim,
    N_zero_bp = n_zero,
    N_2ormore_bp = n_2ormore,
    N_far_from_ref = n_far,
    N_valid = n_valid,
    N_error = n_error,
    Rejection_rate = rejection_rate
  )
  
  cat("\n--- Results Summary ---\n")
  print(summary_df)
  
  if (!is.null(year_table)) {
    cat("\nDetected changepoint years (valid cases):\n")
    print(year_table)
  }
  
  list(
    setting_id = setting_id,
    summary = summary_df,
    year_distribution = year_table,
    all_results = results_list
  )
}

### --- RUN FOR ALL 5 SETTINGS --- ###
all_results <- list()

for (setting_id in 1:5) {
  result <- run_bootstrap_for_setting(
    setting_id = setting_id,
    nsim = nsim,
    refyear = refyear,
    lag_thresh = lag_thresh,
    penalty = penalty,
    tau_star = tau_star,
    years = years,
    n = n
  )
  
  all_results[[setting_id]] <- result
  
  # Save intermediate results
  saveRDS(result, file = sprintf("./Results/bootstrap_size_test_setting%d.Rds", setting_id))
  cat(sprintf("\nSaved results for setting %d\n", setting_id))
}

### --- AGGREGATE AND DISPLAY FINAL RESULTS --- ###
cat("\n\n##############################################\n")
cat("         FINAL SUMMARY - ALL SETTINGS\n")
cat("##############################################\n")

summary_all <- do.call(rbind, lapply(all_results, function(x) x$summary))
print(summary_all)

# Save combined results
saveRDS(all_results, file = "./Results/bootstrap_size_test_all_settings.Rds")
cat("\nSaved combined results to ./Results/bootstrap_size_test_all_settings.Rds\n")

# Print settings parameters for reference
cat("\n--- Settings Parameters ---\n")
print(settings)