library(strucchange)
library(zoo)
source('./MethodCode/PELTtrendARp.R')

library(WeightedPortTest) # for the portmanteau test

robust_change_point_detection <- function(Y, X=NULL, block_size=5, alpha = 0.05, candidate_method = "WBS", 
                                         p_n = NULL, eta_n = 30, kernel_type = "sign") {
  
  # Input validation
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  n <- nrow(Y)
  d <- ncol(Y)
  
  if (n %% 2 != 0) {
    warning("Sample size n is odd. Removing last observation to make n even.")
    Y <- Y[1:(n-1), , drop = FALSE]
    n <- n - 1
  }
  
  m <- n / 2
  if (is.null(p_n)) p_n <- floor(2 * n^(2/5))
  

  # Split data into odd and even indices
  Y_indices <- seq(1, n)
  #O_indices <- seq(1, n, by = 2)
  #E_indices <- seq(2, n, by = 2)
  

  idx <- seq_len(n)
  block_id <- ceiling(idx / block_size)

  O_indices <- idx[block_id %% 2 == 1]
  E_indices <- idx[block_id %% 2 == 0]

  YO_data <- Y[O_indices, , drop = FALSE]
  YE_data <- Y[E_indices, , drop = FALSE]  

  if (!is.null(X)) {
    # If X is a vector (no columns!)
    if (is.vector(X)) {
      XO_data <- matrix(X[O_indices], ncol = 1)
      XE_data <- matrix(X[E_indices], ncol = 1)

    # If X is already a matrix/data.frame with ≥ 1 column
    } else {
      XO_data <- X[O_indices, , drop = FALSE]
      XE_data <- X[E_indices, , drop = FALSE]
    }
  } else {
    XO_data <- NULL
    XE_data <- NULL
  }

  # Step 1: Generate candidate change-points using first half of data
  candidate_cps <- generate_candidate_change_points(YO_data, XO_data, method = candidate_method, 
                                                   p_n = p_n, eta_n = eta_n)
  p_n_actual <- length(candidate_cps)
  
  if (p_n_actual == 0) {
    return(list(change_points = numeric(0), statistics = numeric(0), 
               threshold = Inf, FDP_estimate = numeric(0)))
  }
  
  # Add boundaries
  candidate_cps <- c(0, candidate_cps, m)
  p_n_actual <- length(candidate_cps) - 2

  
  # Step 2: Compute test statistics T_j
  #T_stats <- compute_test_statistics_new(candidate_cps, Y, X,O_indices, E_indices, YO_data, YE_data, XO_data, XE_data, kernel_type)
  T_stats <- compute_test_statistics_old(candidate_cps, Y, X, Y_indices, X_indices, O_indices, E_indices, YO_data, YE_data, XO_data, XE_data, kernel_type)

  # Step 3: Select threshold controlling FDR
  threshold_result <- select_fdr_threshold(T_stats, alpha)
  
  # Step 4: Identify significant change-points
  significant_cps <- which(T_stats >= threshold_result$threshold)
  detected_cps <- candidate_cps[significant_cps + 1]  # +1 to skip initial 0
  
  # Map back to original indices (approximately)
  original_indices <- Y_indices[detected_cps * 2]  # Rough mapping
  
  return(list(
    change_points = original_indices,
    statistics = T_stats,
    threshold = threshold_result$threshold,
    FDP_estimate = threshold_result$FDP_curve,
    candidate_points = candidate_cps[-c(1, length(candidate_cps))] * 2  # Remove boundaries
  ))
}




#' Compute test statistics T_j for all candidate change-points
compute_test_statistics_old <- function(candidate_cps, Y, X, Y_indices, X_indices, O_indices, E_indices, YO_data, YE_data, XO_data, XE_data, kernel_type) {
  p_n <- length(candidate_cps) - 2
  T_stats <- numeric(p_n)
  
  for (j in 1:p_n) {
    # Segment boundaries
    start <- candidate_cps[j] + 1 # +1 because candidate_cps[1] = 0
    cp <- candidate_cps[j + 1]
    end <- candidate_cps[j + 2]
    
    n_j <- cp - start + 1
    n_j1 <- end - cp
    
    if (n_j < 1 || n_j1 < 1) {
      T_stats[j] <- -Inf
      next
    }


    # Map segment indices back to full data
    start_full <- 2 * start - 1  # because of odd/even split we need to multiply by 2 for detrending the full sample
    cp_full <- 2 * cp       # we subtract 1 from start as the first index is odd since the breakpoints are defined on odd indices
    end_full <- 2 * end        # we don't subtract 1 from cp or end as otherwise we never obtain the last index of the full data

    # Full data segment indices
    idx1 <- start_full:cp_full
    idx2 <- (cp_full + 1):end_full

    # Check if indices are in O or E
    in_O_seg1 <- idx1 %in% O_indices
    in_E_seg1 <- idx1 %in% E_indices

    in_O_seg2 <- idx2 %in% O_indices
    in_E_seg2 <- idx2 %in% E_indices

    # Extract segments from full data
    segment1Y <- Y[idx1, drop = FALSE]
    segment2Y <- Y[idx2, drop = FALSE]

    segment1X <- X[idx1, , drop = FALSE]
    segment2X <- X[idx2, , drop = FALSE]

    y1O <- segment1Y[in_O_seg1, drop = FALSE]
    y2O <- segment2Y[in_O_seg2, drop = FALSE]
    X1O <- segment1X[in_O_seg1, , drop = FALSE]
    X2O <- segment2X[in_O_seg2, , drop = FALSE]

    y1E <- segment1Y[in_E_seg1, drop = FALSE]
    y2E <- segment2Y[in_E_seg2, drop = FALSE]
    X1E <- segment1X[in_E_seg1, , drop = FALSE]
    X2E <- segment2X[in_E_seg2, , drop = FALSE]

    # Fit trends on odd and even segments 
    result2O = trendARpsegfit_arimax(y1O, X1O, p=1)
    beta1_hatO <- as.numeric(result2O$coef[1:2])
    phiO <- as.numeric(result2O$coef["ar1"])
    #print(result2O$coef)

    result2E = trendARpsegfit_arimax(y1E, X1E, p=1)
    beta1_hatE <- as.numeric(result2E$coef[1:2])
    phiE <- as.numeric(result2E$coef["ar1"])
    #print(result2E$coef)

    # obtain detrended residuals using SAME beta1_hatO
    residuals1star_betaO <- segment1Y - segment1X %*% beta1_hatO
    residuals2star_betaO <- segment2Y - segment2X %*% beta1_hatO

    residuals1star_betaE <- segment1Y - segment1X %*% beta1_hatE
    residuals2star_betaE <- segment2Y - segment2X %*% beta1_hatE

    u1O <- as.numeric(residuals1star_betaO)  
    u2O <- as.numeric(residuals2star_betaO)
    u1E <- as.numeric(residuals1star_betaE)  
    u2E <- as.numeric(residuals2star_betaE)

    residuals1star_betaO_white <- numeric(length(u1O))
    residuals2star_betaO_white <- numeric(length(u2O))

    residuals1star_betaE_white <- numeric(length(u1E))
    residuals2star_betaE_white <- numeric(length(u2E))
    
    residuals1star_betaO_white[-1] <- u1O[-1] - phiO * u1O[-length(u1O)]       # within seg1
    residuals2star_betaO_white[1] <- u2O[1] - phiO * u1O[length(u1O)]          # boundary step
    residuals2star_betaO_white[-1] <- u2O[-1] - phiO * u2O[-length(u2O)]       # within seg2

    residuals1star_betaE_white[-1] <- u1E[-1] - phiE * u1E[-length(u1E)]       # within seg1
    residuals2star_betaE_white[1] <- u2E[1] - phiE * u1E[length(u1E)]          # boundary step
    residuals2star_betaE_white[-1] <- u2E[-1] - phiE * u2E[-length(u2E)]       # within seg2

    residuals1_O = residuals1star_betaO_white[in_O_seg1, drop = FALSE]
    residuals1_E = residuals1star_betaE_white[in_E_seg1, drop = FALSE]

    residuals2_O = residuals2star_betaO_white[in_O_seg2, drop = FALSE]
    residuals2_E = residuals2star_betaE_white[in_E_seg2, drop = FALSE]

    # scores (up to sign; consistent across segments)
    S1O <- X1O * as.numeric(residuals1_O)                    # (n1-1) x d
    S2O <- X2O * as.numeric(residuals2_O)                    # (n2-1) x d

    S1E <- X1E * as.numeric(residuals1_E)                    # (n1-1) x d
    S2E <- X2E * as.numeric(residuals2_E)                    # (n2-1) x d

    # Compute L_j^E and L_j^O
    result_E <- compute_U_statistic_linear_trend_old(S1E, S2E)
    result_O <- compute_U_statistic_linear_trend_old(S1O, S2O)
    
    L_j_E <- result_E$U_stat
    L_j_O <- result_O$U_stat
    S1 <- result_O$S1
    S2 <- result_O$S2 

    # Estimate covariance matrix
    Sigma_hat <- moving_range_Omega_inv(S1, S2)

    # Compute test statistic
    if (rcond(Sigma_hat) < 1e-10) {
      # Add regularization if matrix is near singular
      Sigma_hat <- Sigma_hat + diag(ncol(XO_data)) * 1e-6
    }

    scale_factor <- (n_j * n_j1) / (n_j + n_j1)
    T_stats[j] <- scale_factor * t(L_j_O) %*% solve(Sigma_hat) %*% L_j_E
  }
  
  return(T_stats)
}

#' Compute U-statistic for linear trend change-point
compute_U_statistic_linear_trend_old <- function(segment1, segment2, include_phi = False, kernel_type = "moment") {
  U_stat <- colMeans(segment1, na.rm = TRUE) - colMeans(segment2, na.rm = TRUE)
  return(list(U_stat=U_stat, S1=segment1, S2=segment2))
}

trendARpsegfit_arimax <- function(y, X, p = 1) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)

  fit <- arima(
    y,
    order = c(p, 0, 0),
    xreg = X,
    include.mean = FALSE,
    method = "ML"
  )

  coef_raw <- coef(fit)

  beta_names <- colnames(X)
  ar_names   <- paste0("ar", seq_len(p))

  coef_all <- c(
    coef_raw[beta_names],
    coef_raw[ar_names]
  )

  beta_hat <- coef_all[beta_names]
  phi_hat  <- coef_all[ar_names]

  # trend fit + residuals
  trend_fit <- as.numeric(X %*% beta_hat)
  trend_res <- y - trend_fit

  # whitened residuals (innovations)
  e <- trend_res
  if (p >= 1) {
    for (k in 1:p) {
      e[(k + 1):n] <- e[(k + 1):n] - phi_hat[k] * trend_res[1:(n - k)]
    }
    e[1:p] <- NA_real_
  }

  list(
    coef       = coef_all,      # Intercept, Trend, ar1
    beta_hat   = beta_hat,
    phi_hat    = phi_hat,
    trend_fit  = trend_fit,
    trend_res  = trend_res,
    ar_fit     = e,             # whitened residuals
    ar_res     = e,
    fit        = fit
  )
}

#' Generate candidate change-points using PELT
generate_candidate_change_points <- function(Y, X, method = "WBS", p_n, eta_n) {
  n <- nrow(data)
    if (method == "PELT") {
      candidates = PELT.trendARp(Y, p=1, pen=0, minseglen=eta_n)
      if (length(candidates) > p_n) {
        print("More candidates than p_n found.")
      }

      return(candidates)
    } else {
      stop("Unknown candidate method for linear_trend-type CPs. Use 'PELT'.")
    }
  }
# Moving-range / difference-based estimator from Zou-Wang-Li style
# Returns Omega_inv_hat (estimate of Omega^{-1})
moving_range_Omega_inv <- function(S1, S2) {
  S1 <- as.matrix(S1)
  S2 <- as.matrix(S2)
  stopifnot(ncol(S1) == ncol(S2))

  S <- rbind(S1, S2)          # m x d
  m <- nrow(S)
  d <- ncol(S)

  if (m < 5) stop("Need at least 5 rows total to use 2-step differences (2i-1, 2i-3).")

  # indices: 2i-1 for i=2..floor((m+1)/2)  => j = 3,5,7,...
  j_max <- if (m %% 2 == 1) m else (m - 1)
  j_seq <- seq(3, j_max, by = 2)  # 3,5,...,j_max
  n_terms <- length(j_seq)
  if (n_terms < 1) stop("Not enough terms for the moving-range estimator.")

  Omega_inv_hat <- matrix(0, d, d)

  for (j in j_seq) {
    diff <- S[j, ] - S[j - 2, ]      # (S_{2i-1} - S_{2i-3})
    Omega_inv_hat <- Omega_inv_hat + tcrossprod(diff)  # diff %*% t(diff)
  }

  Omega_inv_hat <- Omega_inv_hat / (2 * (n_terms))  # 2(m-1) in their notation; here n_terms = m-1 over odd subseq
  Omega_inv_hat
}



#' Select threshold controlling FDR using symmetry property
select_fdr_threshold <- function(T_stats, alpha) {
  # Remove infinite values
  finite_idx <- is.finite(T_stats)
  T_finite <- T_stats[finite_idx]
  
  if (length(T_finite) == 0) {
    return(list(threshold = Inf, FDP_curve = numeric(0)))
  }
  
  # Generate threshold candidates
  t_candidates <- sort(unique(c(0, abs(T_finite))))
  if (length(t_candidates) == 0) {
    return(list(threshold = Inf, FDP_curve = numeric(0)))
  }
  
  FDP_estimates <- numeric(length(t_candidates))
  
  for (i in 1:length(t_candidates)) {
    t <- t_candidates[i]
    
    # Count positive statistics above threshold
    R_plus <- sum(T_finite >= t)
    
    # Count negative statistics below -t (using symmetry)
    R_minus <- sum(T_finite <= -t)
    
    # Estimated FDP
    #print("FDP .1 instead of 1")
    FDP_estimates[i] <- (.05 + R_minus) / max(R_plus, 1)
  }
  
  # Find the largest threshold that controls FDR
  valid_thresholds <- t_candidates[FDP_estimates <= alpha]
  
  if (length(valid_thresholds) > 0) {
    threshold <- min(valid_thresholds)
  } else {
    threshold <- Inf
  }
  
  return(list(
    threshold = threshold,
    FDP_curve = data.frame(threshold = t_candidates, FDP = FDP_estimates)
  ))
}

simulate_piecewise_linear_discontinuous <- function(
  n = 2000, K_n = 10,
  error = "gaussian", sigma = 12,
  AR_phi = 0.4
) {
  ## 1) Generate true change-points
  tau <- floor((1:K_n) * n / (K_n + 1)) + runif(K_n, -n^(1/4), n^(1/4))
  tau <- sort(pmax(2, pmin(n - 1, round(tau))))
  tau <- unique(tau)
  K_n <- length(tau)

  breaks <- c(tau, n)
  starts <- c(1, tau + 1)

  ## 2) Slopes + intercepts per segment (DIScontinuous)
  beta1_grid <- c(0.01, -0.02, 0.03, -0.015, 0.025,
                  -0.01, 0.05, -0.04, 0.02, -0.03, 0.01)

  beta0 <- beta1 <- numeric(n)

  # segment-wise trend used ONLY to generate Y
  Xseg <- numeric(n)

  jump_sd <- 1
  b0_grid <- rnorm(K_n + 1, mean = 0, sd = jump_sd)

  for (k in 1:(K_n + 1)) {
    s <- starts[k]
    e <- breaks[k]

    idx1 <- ((k - 1) %% length(beta1_grid)) + 1
    b1_k <- sample(beta1_grid, size = 1)
    b0_k <- b0_grid[k]

    # reset-to-1 trend within segment
    Xk <- seq_len(e - s + 1)
    Xseg[s:e] <- Xk

    beta0[s:e] <- b0_k
    beta1[s:e] <- b1_k
  }

  ## 3) i.i.d innovations u_t
  u <- switch(error,
    "gaussian" = rnorm(n),
    "t3"       = rt(n, df = 3) / sqrt(3),
    "mix"      = 0.8 * rnorm(n) + 0.2 * rnorm(n, sd = 3),
    "chi2"     = (rchisq(n, df = 3) - 3) / sqrt(6),
    "cauchy"   = rcauchy(n) * 0.3,
    stop("Unknown error option.")
  )

  ## 4) AR(1) filtering (for ALL error types)
  err <- numeric(n)
  err[1] <- u[1]
  for (i in 2:n) err[i] <- AR_phi * err[i - 1] + u[i]
  err <- err / sd(err)  # optional standardization

  ## 5) Construct response using segment-wise trend
  Y <- beta0 + beta1 * Xseg + sigma * err

  ## 6) Return GLOBAL trend as covariate
  X <- 1:n
  Xmat <- cbind(Intercept = 1, Trend = X)

  list(
    Y = Y,
    X = Xmat,        # <-- global trend
    Xseg = Xseg,     # <-- segment-wise trend used in DGP (kept for diagnostics)
    beta0 = beta0,
    beta1 = beta1,
    tau_true = tau
  )
}

dist_cp_mat <- function(T_hat, T_star) {
  T_hat  <- as.numeric(T_hat)
  T_star <- as.numeric(T_star)

  if (length(T_star) == 0) return(0)
  if (length(T_hat)  == 0) return(Inf)

  D <- abs(outer(T_hat, T_star, "-"))   # rows: T_hat, cols: T_star
  mean(apply(D, 2, min))
}

refine<-function(th_L,Location_e){
###### output: refined MOPS change-point Location and informative points #######        
        Location_p<-numeric()
        for(i in 1:length(Location_e)){
        dista=abs(setdiff(th_L,Location_p)-(Location_e)[i])
        d_min=min(dista)             
        Location_p[length(Location_p)+1]=setdiff(th_L,Location_p)[which.min(dista)]            
              } 
        infor_tau<-sapply(1:K, function(t){
                        infor <- Location_p[which.min(abs(Location_p-Location_e[t]))]
                        return(infor)                          
                         }              
                  )
        Re<-list(Location_p=Location_p,infor_tau=infor_tau)          
        return(Re)         
      }
  


meas_index<-function(Location,infor_tau,tau){
##### output: FDR,TPR,Pa,estimate number,dist ##########
        numb_dif<-length(Location)
        trueD<-intersect(Location,infor_tau)                        
        fdr<-(numb_dif-length(trueD))/max(numb_dif,1)
        tdr<-length(trueD)/K
        Pa<-ifelse(tdr==1,1,0)  
        if(length(Location)>0)             
         {temp_dist<-sapply(1:K, function(t){
                          min(abs(Location-tau[t])) 
                          }
                   )
          dist<-mean(temp_dist)
         } else {dist<-Inf}
        return(c(fdr,tdr,Pa,numb_dif,dist))
}




# load raw monthly data
#file_raw_monthly = 'ObservedData.csv'
#df_raw_monthly = read.csv(file_raw_monthly)#load(file)
#df_raw_monthly = df_raw_monthly#[241:length(df_raw_monthly$t),]
#Y_raw_monthly = df_raw_monthly$NOAA
#dates_raw_monthly = df_raw_monthly$t

# load adjusted data
#file_adj_monthly = 'temperature_anomalies_adj_monthly.RData'
#load(file_adj_monthly)
#df_adj_monthly <- subset(Tanom_annual_df, time >= 1880)
#names(df_adj_monthly)[1] <- "t"
#Y_adj_monthly <- df_adj_monthly$Berk.adj
#dates_adj = df_adj_monthly$t


# load adj yearly data
file_adj_yearly = './data/temperature_anomalies_adj.RData'
load(file_adj_yearly)
df_adj_yearly <- subset(Tanom_annual_df, Year >= 1880)
names(df_adj_yearly)[1] <- "t"
Y_adj_yearly <- df_adj_yearly$NASA.adj
dates_adj_yearly = df_adj_yearly$t

# For adj: NASA / NOAA blocksize=3, Berk/HadCRU blocksize=2


# load raw yearly data
file_raw_yearly = './data/temperature_anomalies.RData'
load(file_raw_yearly)
df_raw_yearly <- subset(Tanom_annual_df, year >= 1880)
names(df_raw_yearly)[1] <- "t"
Y_raw_yearly <- df_raw_yearly$NASA
dates_raw_yearly = df_raw_yearly$t

# For raw: HadCRU/NOAA/NASA blocksize=2, Berk blocksize=4

Y = Y_raw_yearly### Y_adj or Y_raw

dates = dates_raw_yearly ### fix at dates_raw
if (length(dates) %% 2 == 1) {
  dates <- dates[-length(dates)]
}
if (length(Y) %% 2 == 1) {
  Y <- Y[-length(Y)]
}
n <- length(Y)
X <- cbind(
  Intercept = rep(1, length(Y)),
    X_year = (1:length(Y)) / 12
)
eta_n <- 9#min(60, floor(n^(1/2)))

# Compare with moment kernel
result_moment <- robust_change_point_detection(
  Y,
  X, 
  block_size = 2,
  alpha = .5, 
  candidate_method = "PELT",
  kernel_type = "moment", eta_n=eta_n,
)
Location_e <- result_moment$change_points # detected change-points from moment kernel T
fittrend_e = fit.trendARp(Y,Location_e,p=1,plot=T,add.ar=F,fit=T,
                            dates=dates)#
K <- length(Location_e)
candidate_cps_refined <- generate_candidate_change_points(Y, X, method = "PELT", 
                                                   p_n = floor(2 * n^(2/5)), eta_n = eta_n) # A(Z)


R_MOPS <- refine(candidate_cps_refined, Location_e)
Location_p <- R_MOPS$Location_p #Refined MOPS change-point Location T(Z)

fittrend_p = fit.trendARp(Y, Location_p, p=1,plot=T,add.ar=F,fit=T,
                            dates=dates)#
LAG=20
resid <- Y - fittrend_p$fit_total
Weighted.Box.test(resid,lag=LAG,type="Ljung",fitdf=1)$p.value



