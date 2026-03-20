##########################
###### Vantage Simulations
##########################

library(parallel)
library(doParallel)
library(foreach)

# Load data
load("./data/temperature_anomalies.RData")
ANOM <- Tanom_annual_df[,]

# Fit the null (no change) distribution to HadCRUT (1970:2023)
y <- ANOM[121:174, 4]
times <- 1:length(y)
shortfitHad <- arima(y, order = c(1, 0, 0), xreg = times)

# Extract point estimates
phi_hat <- shortfitHad$coef["ar1"]
intercept_hat <- shortfitHad$coef["intercept"]
slope_hat <- shortfitHad$coef["times"]
sigma2_hat <- shortfitHad$sigma2

# Extract variance-covariance matrix for regression coefficients
vcov_mat <- shortfitHad$var.coef

# Standard errors for regression coefficients
se_phi <- sqrt(vcov_mat["ar1", "ar1"])
se_intercept <- sqrt(vcov_mat["intercept", "intercept"])
se_slope <- sqrt(vcov_mat["times", "times"])

# Standard error for sigma2 (asymptotic: Var(sigma2_hat) ≈ 2*sigma4/n)
n_obs <- length(y)
se_sigma2 <- sqrt(2 * sigma2_hat^2 / n_obs)
se_sd <- se_sigma2 / (2 * sqrt(sigma2_hat))  # Delta method: se(sd) = se(sigma2) / (2*sd)

# Point estimates
sd_hat <- sqrt(sigma2_hat)

# 95% CI bounds (z = 1.96)
z <- 1.96

# Define 5 settings
settings <- data.frame(
  setting = 1:5,
  phi = c(
    phi_hat,                          # 1: center
    phi_hat + z * se_phi,             # 2: high phi
    phi_hat + z * se_phi,             # 3: high phi
    phi_hat - z * se_phi,             # 4: low phi
    phi_hat - z * se_phi              # 5: low phi
  ),
  intercept = c(
    intercept_hat,                    # 1: center
    intercept_hat + z * se_intercept, # 2: high intercept
    intercept_hat - z * se_intercept, # 3: low intercept
    intercept_hat + z * se_intercept, # 4: high intercept
    intercept_hat - z * se_intercept  # 5: low intercept
  ),
  slope = c(
    slope_hat,                        # 1: center
    slope_hat + z * se_slope,         # 2: high slope
    slope_hat - z * se_slope,         # 3: low slope
    slope_hat - z * se_slope,         # 4: low slope
    slope_hat + z * se_slope          # 5: high slope
  ),
  sd = c(
    sd_hat,                           # 1: center
    sd_hat + z * se_sd,               # 2: high sd
    sd_hat - z * se_sd,               # 3: low sd
    sd_hat - z * se_sd,               # 4: low sd
    sd_hat + z * se_sd                # 5: high sd
  )
)

# Ensure sd is positive (lower bound)
settings$sd <- pmax(settings$sd, 0.01)

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

# Parallel processing
ns <- 33:71
n_sim <- 2
n_cores <- detectCores() - 1
cat(sprintf("Using %d cores\n", n_cores))

# Register parallel backend
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export required functions and data to workers
clusterExport(cl, c("simulateone", "gettmax", "settings", "ns"))

# Create a task grid: all combinations of settings and ns
tasks <- expand.grid(setting_id = 1:nrow(settings), n = ns)

# Run all tasks in parallel
results <- foreach(task_idx = 1:nrow(tasks), .packages = c("stats"), .errorhandling = "pass") %dopar% {
  task <- tasks[task_idx, ]
  setting_id <- task$setting_id
  n_val <- task$n
  setting_row <- settings[setting_id, ]
  
  # Run the simulation
  QN <- tryCatch({
    gettmax(n = n_val, interc = setting_row$intercept, slope = setting_row$slope, sd = setting_row$sd, phi = setting_row$phi, n_sim = n_sim)
  }, error = function(e) {
    cat(sprintf("Error for setting %d, n=%d: %s\n", setting_id, n_val, e$message))
    NA_real_
  })
  
  list(setting_id = setting_id, n = n_val, QN = QN)
}

# Stop cluster
stopCluster(cl)

# Aggregate results per setting
for (setting_id in 1:nrow(settings)) {
  setting_results <- do.call(rbind, lapply(results, function(res) {
    if (res$setting_id == setting_id) {
      data.frame(n = res$n, QN = res$QN)
    }
  }))
  
  # Save results for the setting
  save(setting_results, file = sprintf("./Results/TmaxquantHad_setting%d.Rdata", setting_id))
  cat(sprintf("Saved results for setting %d\n", setting_id))
}

cat("All tasks completed!\n")