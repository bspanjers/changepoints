##########################
###### Vantage Simulations
##########################

library(WeightedPortTest)
library(parallel)
library(doParallel)
library(foreach)

# Load data
load("./data/temperature_anomalies.RData")
ANOM <- Tanom_annual_df[,]

# Define settings
ns <- 33:71
n_sim <- 100000
n_cores <- min(128, detectCores() - 1)  # Use up to 128 cores
cat(sprintf("Using %d cores\n", n_cores))

# Intermediate save directory
save_dir <- "./Results"
if (!dir.exists(save_dir)) dir.create(save_dir)

# Define the simulation function
simulateone <- function(index = 1, n, interc, slope, sd, phi) {
  seg <- 1:n
  y <- interc + slope * seg + arima.sim(n = n, sd = sd, model = list(ar = phi))
  min_len <- max(round(0.1 * n), 2)
  max_len <- min(n - round(0.1 * n), n - 2)
  seglen <- max_len:min_len
  stats <- numeric(length(seglen))
  
  for (i in seq_along(seglen)) {
    seg2 <- (n - seglen[i] + 1):n
    seg1 <- 1:(n - seglen[i])
    X <- cbind(rep(1, n), c(seg1, rep(seg1[length(seg1)], seglen[i])), c(rep(0, (n - seglen[i])), seg2 - seg2[1] + 1))
    vec <- c(0, 0, -1, 1)
    armafit <- arima(y, xreg = X, order = c(1, 0, 0), include.mean = FALSE)
    sddiff <- sqrt(t(vec) %*% armafit$var.coef %*% vec)
    stats[i] <- (armafit$coef[3] - armafit$coef[4]) / sddiff
  }
  
  return(max(abs(stats)))
}

gettmax <- function(n, interc, slope, sd, phi, n_sim) {
  set.seed(n)
  hold <- sapply(1:n_sim, simulateone, n = n, interc = interc, slope = slope, sd = sd, phi = phi)
  return(quantile(hold, 0.95))
}

# Parallel processing for all settings
run_simulation <- function(settings, ns, n_sim, save_dir) {
  # Register parallel backend
  cl <- makeCluster(n_cores, type = "FORK")  # Use FORK for efficiency on Linux
  registerDoParallel(cl)
  
  # Process each setting
  foreach(i = 1:nrow(settings), .packages = c("stats"), .errorhandling = "pass") %dopar% {
    setting_row <- settings[i, ]
    setting_id <- setting_row$setting
    
    # Check if results already exist
    save_file <- sprintf("%s/TmaxquantHad_setting%d.Rdata", save_dir, setting_id)
    if (file.exists(save_file)) {
      cat(sprintf("Setting %d already completed. Skipping...\n", setting_id))
      return(NULL)
    }
    
    # Run simulations for all ns
    QNs <- sapply(ns, function(n) {
      cat(sprintf("Processing setting %d, n=%d...\n", setting_id, n))  # Print progress for each n
      tryCatch({
        gettmax(n = n, interc = setting_row$intercept, slope = setting_row$slope, sd = setting_row$sd, phi = setting_row$phi, n_sim = n_sim)
      }, error = function(e) {
        cat(sprintf("Error for setting %d, n=%d: %s\n", setting_id, n, e$message))
        NA_real_
      })
    })
    
    # Save results
    TmaxquantHad_setting <- data.frame(ns = ns, QNs = QNs)
    save(TmaxquantHad_setting, file = save_file)
    cat(sprintf("Saved results for setting %d\n", setting_id))
  }
  # Stop cluster
  stopCluster(cl)
}

# Define settings (already provided in your script)
settings <- data.frame(
  setting = 1:5,
  phi = c(phi_hat, phi_hat + z * se_phi, phi_hat + z * se_phi, phi_hat - z * se_phi, phi_hat - z * se_phi),
  intercept = c(intercept_hat, intercept_hat + z * se_intercept, intercept_hat - z * se_intercept, intercept_hat + z * se_intercept, intercept_hat - z * se_intercept),
  slope = c(slope_hat, slope_hat + z * se_slope, slope_hat - z * se_slope, slope_hat - z * se_slope, slope_hat + z * se_slope),
  sd = c(sd_hat, sd_hat + z * se_sd, sd_hat - z * se_sd, sd_hat - z * se_sd, sd_hat + z * se_sd)
)

# Run the simulation
run_simulation(settings, ns, n_sim, save_dir)

cat("All settings completed!\n")