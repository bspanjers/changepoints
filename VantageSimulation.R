##########################
###### Vantage Simulations
##########################

library(WeightedPortTest)
library(shape) # for colour palette
library(lattice) # for levelplot
library(parallel)
library(doParallel)
library(foreach)
#load('Results/resultstrend.Rdata')
#load('Results/resultstrendar.Rdata')
#load('Results/resultstrendFIXar4.Rdata')
#load('Results/resultsjointrend.Rdata')

load("./data/temperature_anomalies.RData")
ANOM <- Tanom_annual_df[,]



#Fit the null (no change) distribution to HadCRUT (1970:2023)
y=ANOM[121:174,4]
times=1:length(y)
shortfitHad=arima(y,order=c(1,0,0),xreg=times)
#Test for remaining autocorrelation up to lag 10; pvalue>.05 indicates AR(1) is adequate
#Weighted.Box.test(shortfitHad$residuals,lag=10,fitdf=1,type="Ljung")

# -------------------------
# Generate 5 settings within 95% confidence region
# -------------------------

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

# Create 5 settings: center + 4 corners varying phi, slope, and sd
# Setting 1: Point estimates (center)
# Setting 2: High phi, high slope, high sd
# Setting 3: High phi, low slope, low sd
# Setting 4: Low phi, high slope, low sd
# Setting 5: Low phi, low slope, high sd

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

# Print the settings
cat("\n========================================\n")
cat("Point Estimates:\n")
cat(sprintf("  phi (AR1):     %.4f (SE: %.4f)\n", phi_hat, se_phi))
cat(sprintf("  intercept:     %.4f (SE: %.4f)\n", intercept_hat, se_intercept))
cat(sprintf("  slope:         %.4f (SE: %.4f)\n", slope_hat, se_slope))
cat(sprintf("  sigma (sd):    %.4f (SE: %.4f)\n", sd_hat, se_sd))
cat("========================================\n\n")

cat("95% Confidence Intervals:\n")
cat(sprintf("  phi:       [%.4f, %.4f]\n", phi_hat - z*se_phi, phi_hat + z*se_phi))
cat(sprintf("  intercept: [%.4f, %.4f]\n", intercept_hat - z*se_intercept, intercept_hat + z*se_intercept))
cat(sprintf("  slope:     [%.4f, %.4f]\n", slope_hat - z*se_slope, slope_hat + z*se_slope))
cat(sprintf("  sd:        [%.4f, %.4f]\n", max(0.01, sd_hat - z*se_sd), sd_hat + z*se_sd))
cat("========================================\n\n")

cat("5 Settings within 95% Confidence Region:\n")
print(settings, digits = 4)


simulateone=function(index=1,n,interc=-0.0662,slope=.019,sd=.09,phi=0.178){
  #Code to simulate one time series and one Tmax value
  seg=1:n
  y=interc+slope*(seg)+arima.sim(n=n,sd=sd,model=list(ar=phi))
  min=max(round(.1*n),2)
  max=min(n-round(.1*n),n-2)
  seglen=max:min
  stats=seglen
  sddiff=stats
  for(i in 1:length(seglen)){
    seg2=(n-seglen[i]+1):n
    seg1=1:(n-seglen[i])
    X=cbind(rep(1,n),c(seg1,rep(seg1[length(seg1)],seglen[i])),c(rep(0,(n-seglen[i])),seg2-seg2[1]+1))
    vec=c(0,0,-1,1)
    armafit=arima(y,xreg=X,order=c(1,0,0),include.mean=FALSE)
    sddiff=sqrt(t(vec)%*%armafit$var.coef%*%vec)
    stats[i]=(armafit$coef[3]-armafit$coef[4])/sddiff
  }
  return((max(abs(stats))))
}

gettmax=function(n,interc=-0.0662,slope=.019,sd=.09,phi=0.178, n_sim=100000){
  #This function will simulate 100000 Tmax values and get .95 quantile QN
  set.seed(n)
  hold=sapply(1:n_sim,simulateone,n=n,interc=interc,slope=slope,sd=sd,phi=phi)
  return(quantile(hold,.95))			
}

#N values for AMOC changepoint detection 2023 through 2040
#ns=(2023-1969):(2040-1969)#from beaulieu
#ns2= (2023-1990):(2023-1970)
#The following line takes more than 24 hours to run; It simulates 0.95 quantiles of TMAX null distribution
#QNs=sapply(ns,gettmax,interc=shortfitHad$coef[2],slope=shortfitHad$coef[3],phi=shortfitHad$coef[1])
#QNs2=sapply(ns2,gettmax,interc=shortfitHad$coef[2],slope=shortfitHad$coef[3],phi=shortfitHad$coef[1])
#TmaxquantHad=as.data.frame(cbind(ns,QNs))
#TmaxquantHad2=as.data.frame(cbind(ns2,QNs2))
#save(TmaxquantHad,file="./Results/TmaxquantHad.Rdata")
#save(TmaxquantHad2,file="./Results/TmaxquantHad2.Rdata")
#load("./Results/TmaxquantHad.Rdata")
#QNs=TmaxquantHad$QNs


ns <- 33:71
n_sim <- 2
# Number of cores to use
n_cores <- detectCores() - 1
cat(sprintf("Using %d cores\n", n_cores))

# Register parallel backend
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export required functions and data to workers
clusterExport(cl, c("simulateone", "gettmax", "settings", "ns"))

# Run all settings in parallel
results_list <- foreach(
  i = 1:nrow(settings),
  .packages = c("stats"),
  .errorhandling = "pass"
) %dopar% {
  
  setting_row <- settings[i, ]
  
  QNs <- sapply(ns, function(n) {
    gettmax(
      n = n,
      interc = setting_row$intercept,
      slope = setting_row$slope,
      sd = setting_row$sd,
      phi = setting_row$phi,
      n_sim=n_sim
    )
  })
  
  TmaxquantHad_setting <- as.data.frame(cbind(ns = ns, QNs = QNs))
  
  list(
    setting_id = i,
    setting_params = setting_row,
    TmaxquantHad = TmaxquantHad_setting
  )
}

# Stop cluster
stopCluster(cl)

# Save results
for (i in seq_along(results_list)) {
  res <- results_list[[i]]
  if (!inherits(res, "error")) {
    TmaxquantHad_setting <- res$TmaxquantHad
    save(TmaxquantHad_setting, 
         file = sprintf("./Results/TmaxquantHad_setting%d.Rdata", res$setting_id))
    
    setting_info <- list(
      setting = res$setting_params,
      TmaxquantHad = res$TmaxquantHad
    )
    save(setting_info, 
         file = sprintf("./Results/TmaxquantHad_setting%d_full.Rdata", res$setting_id))
    
    cat(sprintf("Saved setting %d\n", res$setting_id))
  } else {
    cat(sprintf("Error in setting %d: %s\n", i, res$message))
  }
}

cat("All settings completed!\n")


