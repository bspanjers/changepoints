
library(WeightedPortTest) # for the portmanteau test
library(EnvCpt) # for changepoint detection
library(changepoint) # for changepoint detection
library(parallel)

# load data
load('./data/temperature_anomalies.RData')
# use Tanom_annual_df, matrix, 6 columns, first year
# at the time of analysis, the Japan Met dataset had not been updated and is not included in the paper

# setup for PORTMANTEAU results - Table S1
names = c("NASA","Japan Met","HadCRUT","NOAA","Berkeley")
FGstatdisc = matrix(data=NA,nrow=4,ncol=5,dimnames=list(c("IID","GlobalAR(1)","GlobalAR(4)","ChangingAR(1)"),names))
FGstatcont = matrix(data=NA,nrow=4,ncol=5,dimnames=list(c("IID","GlobalAR(1)","GlobalAR(4)","ChangingAR(1)"),names))
LAG=20

source('./MethodCode/PELTtrendARpJOIN.R')


########### Continuous model AR(1)
trendarjoin=list()
mresiduals=list()
fits=list()
dates=list()

data=Tanom_annual_df[,c(1,4)][!is.na(Tanom_annual_df[,4]),]
y = data[,2]
years = data[,1]
n=nrow(data)

#Model fit
itrendarjoin=PELT.trendARpJOIN(y, p=1,pen=4*log(n),minseglen=10)
fittrend = fit.trendARpJOIN(y, itrendarjoin,p=1,dates=years,plot=T,add.ar=F,fit=T,
                            title=names(Tanom_annual_df[4]),pred=F)# get fit without AR - to visualize trend segments
fits = fittrend$fit
dates = fittrend$dates
fittrendAR = fit.trendARpJOIN(y,itrendarjoin,p=1,dates=years,plot=F,add.ar=T,fit=T,
                            title=names(Tanom_annual_df[4])) #get fit with AR - to compute residuals below
yeartrendarjoin=years[itrendarjoin]  # put cpts in terms of year
cat(paste(names(Tanom_annual_df),':TrendAR1join \n'))
print(fittrend$coeffs)




#simulate path
simu = simulate_trendARpJOIN(y, fittrendAR, itrendarjoin)

#Model fit
itrendarjoin_simu = PELT.trendARpJOIN(simu, p=1,pen=4*log(n),minseglen=10)
fittrend_simu = fit.trendARpJOIN(simu, itrendarjoin_simu,p=1,dates=years,plot=T,add.ar=F,fit=T,
                            title=names(Tanom_annual_df[4]),pred=F)# get fit without AR - to visualize trend segments
fits = fittrend_simu$fit
dates = fittrend_simu$dates
fittrendAR_simu = fit.trendARpJOIN(simu,itrendarjoin_simu,p=1,dates=years,plot=F,add.ar=T,fit=T,
                            title=names(Tanom_annual_df[4])) #get fit with AR - to compute residuals below
yeartrendarjoin=years[itrendarjoin_simu]  # put cpts in terms of year
cat(paste(names(Tanom_annual_df),':TrendAR1join \n'))
print(fittrend_simu$coeffs)


simu = simulate_trendARpJOIN(y, fittrendAR, itrendarjoin)


yafterbreak=simu[121: length(simu)]
n=length(yafterbreak)
seg=1:n
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
  armafit=arima(yafterbreak,xreg=X,order=c(1,0,0),include.mean=FALSE)
  sddiff=sqrt(t(vec)%*%armafit$var.coef%*%vec)
  stats[i]=(armafit$coef[3]-armafit$coef[4])/sddiff
}
#TMAX for 1970-2023
cat((max(abs(stats))),"\n")	
i=order(abs(stats),decreasing=TRUE)[1] 
cat("End of segment 1=",simu[121:174][n-seglen[i]] ,"\n")



# 1. true changepoint (e.g. 1973)
# 2. estimated changepoint (from PELT)
# 3. fixed changepoint (1970)


########## "TRUE" CHANGEPOINT: 1973 #################
nsim <- 2000
threshold <- 3.146627   # TMAX threshold recalculated for 1973
count_exceed <- 0  # counter

for(rep in 1:nsim){

  ### ---- simulate from fitted model ----
  simu <- simulate_trendARpJOIN(y, fittrendAR, itrendarjoin)

  ### ---- extract part AFTER 1973 ----
  yafterbreak <- simu[124:length(simu)]

  n <- length(yafterbreak)
  seg <- 1:n
  min <- max(round(.1*n),2)
  max <- min(n-round(.1*n),n-2)
  seglen <- max:min
  stats <- seglen

  ### ---- compute TMAX ----
  for(i in 1:length(seglen)){
    seg2 <- (n-seglen[i]+1):n
    seg1 <- 1:(n-seglen[i])

    X <- cbind(rep(1,n),
               c(seg1,rep(seg1[length(seg1)],seglen[i])),
               c(rep(0,(n-seglen[i])),seg2-seg2[1]+1))

    vec <- c(0,0,-1,1)
    armafit <- arima(yafterbreak, xreg=X, order=c(1,0,0), include.mean=FALSE)
    sddiff <- sqrt(t(vec) %*% armafit$var.coef %*% vec)
    stats[i] <- (armafit$coef[3] - armafit$coef[4]) / sddiff
  }

  Tmax <- max(abs(stats))

  ### ---- count exceedance ----
  if(Tmax > threshold) count_exceed <- count_exceed + 1
}

### ---- report result ----
cat("\n----------------------------------------\n")
cat("Out of", nsim, "simulations:\n")
cat("TMAX >", threshold, " occurred", count_exceed, "times.\n")
cat("Estimated probability =", count_exceed/nsim, "\n")
cat("----------------------------------------\n")



# ...existing code...

##############################################################################
###               BOOTSTRAP SIZE TEST FOR JOIN MODEL (PARALLEL)            ###
##############################################################################

library(foreach)
library(doParallel)
library(progress)

### --- Load quantile table used for thresholding --- ###
load("./Results/TmaxquantHad_combined.Rdata")
QNs <- merged

### --- PARAMETERS --- ###
nsim       <- 20000
refyear    <- 1973
lag_thresh <- 15
penalty    <- 4 * log(n)

### --- SETUP PARALLEL BACKEND --- ###
n_cores <- detectCores() - 1  # leave one core free
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export necessary objects to workers
clusterExport(cl, c("y", "fittrendAR", "itrendarjoin", "years", "QNs", 
                    "penalty", "refyear", "lag_thresh",
                    "simulate_trendARpJOIN", "PELT.trendARpJOIN"))
clusterEvalQ(cl, source('./MethodCode/PELTtrendARpJOIN.R'))

cat("Running", nsim, "simulations on", n_cores, "cores...\n")

### --- SINGLE SIMULATION FUNCTION --- ###
run_one_sim <- function(i) {
  
  simu <- simulate_trendARpJOIN(y, fittrendAR, itrendarjoin)
  
  bp <- suppressWarnings(tryCatch(
    PELT.trendARpJOIN(simu, p=1, pen=penalty, minseglen=10),
    error = function(e) NULL,
    warning = function(w) NULL
  ))
  
  # Return early for failures
  if (is.null(bp)) return(list(status = "pelt_fail"))
  if (length(bp) == 0) return(list(status = "zero"))
  if (length(bp) >= 2) return(list(status = "2ormore"))
  
  year_est <- years[bp]
  if (abs(year_est - refyear) > lag_thresh) return(list(status = "far"))
  
  # Valid 1-break case
  thr <- QNs[QNs$ns == 2023 - year_est + 1, ]$QNs
  
  yafter <- simu[bp:length(simu)]
  n2     <- length(yafter)
  seglen <- max(round(.1*n2), 2) : min(n2 - round(.1*n2), n2 - 2)
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
      arima(yafter, xreg = Xreg, order = c(1,0,0), include.mean = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    sdd <- sqrt(t(vec) %*% fit$var.coef %*% vec)
    stats[j] <- (fit$coef[3] - fit$coef[4]) / sdd
  }
  
  Tmax <- max(abs(stats), na.rm = TRUE)
  
  list(status = "valid", year = year_est, Tmax = Tmax, thr = thr)
}

### --- RUN PARALLEL WITH PROGRESS --- ###
# Use pbapply for progress bar with parallel
library(pbapply)
pboptions(type = "timer")  # shows ETA

results_list <- pblapply(1:nsim, run_one_sim, cl = cl)

stopCluster(cl)

### --- AGGREGATE RESULTS --- ###
count_pelt_fail <- sum(sapply(results_list, function(x) x$status == "pelt_fail"))
count_zero      <- sum(sapply(results_list, function(x) x$status == "zero"))
count_2ormore    <- sum(sapply(results_list, function(x) x$status == "2ormore"))
count_far       <- sum(sapply(results_list, function(x) x$status == "far"))

valid_results <- Filter(function(x) x$status == "valid", results_list)

# Build year-wise counts
valid_1break  <- list()
reject_1break <- list()

for (res in valid_results) {
  yr <- as.character(res$year)
  if (!yr %in% names(valid_1break)) {
    valid_1break[[yr]]  <- 0
    reject_1break[[yr]] <- 0
  }
  valid_1break[[yr]] <- valid_1break[[yr]] + 1
  if (res$Tmax > res$thr) reject_1break[[yr]] <- reject_1break[[yr]] + 1
}

##############################################################################
###                             SUMMARY OUTPUT                             ###
##############################################################################

rej_1break <- sapply(sort(names(valid_1break)), function(yr) {
  reject_1break[[yr]] / valid_1break[[yr]]
})

results <- list(
  no_break          = count_zero / nsim,
  more_than_two     = count_2ormore / nsim,
  far_from_1973     = count_far / nsim,
  rejection_by_year = rej_1break,
  pelt_failures     = count_pelt_fail,
  counts = list(
    zero         = count_zero,
    count_2ormore        = count_2ormore,
    far          = count_far,
    valid_1break = valid_1break[sort(names(valid_1break))]
  )
)

saveRDS(results, file = "simulation_results_beaulieu.Rds")

cat("\n=========== SUMMARY ===========\n")
cat("> PELT failed (error/warning):", count_pelt_fail, "\n")
cat("> No breaks found            :", count_zero, "\n")
cat("> >=2 breaks found           :", count_2ormore, "\n")
cat("> 1 break but too far        :", count_far, "\n")
cat("> Valid 1-break cases        :", length(valid_results), "\n")
cat("================================\n")

print(results$rejection_by_year)