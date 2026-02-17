library(strucchange)
library(zoo)
source('./MethodCode/PELTtrendARp.R')

library(WeightedPortTest) # for the portmanteau test

robust_change_point_detection <- function(Y, X=NULL, p=1, block_size=5, alpha = 0.05, candidate_method = "WBS", 
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
  candidate_cps <- generate_candidate_change_points(YO_data, XO_data, p=p, method = candidate_method, 
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
  T_stats <- compute_test_statistics(candidate_cps, Y, X,
                                    O_indices, E_indices,
                                    kernel_type, p = 1)
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



compute_test_statistics <- function(candidate_cps, Y, X,
                                    O_indices, E_indices,
                                    kernel_type = NULL, p = 1) {

  # -------- helpers --------
  safe_beta2 <- function(coefs) {
    if (is.null(coefs) || length(coefs) < 2) return(c(0, 0))
    as.numeric(coefs[1:2])
  }
  safe_phi1 <- function(coefs, p) {
    if (p < 1) return(0)
    if (is.null(coefs) || is.null(names(coefs))) return(0)
    if (!("ar1" %in% names(coefs))) return(0)
    as.numeric(coefs["ar1"])
  }

  # Build AR(1) innovations eps_t and corresponding "whitened" regressors Xtilde_t
  # for a concatenation of (segment1, segment2), with boundary step using last u1/x1.
  make_eps_and_Xtilde_ar1 <- function(u1, X1, u2, X2, phi) {
    n1 <- length(u1); n2 <- length(u2)

    eps1 <- rep(NA_real_, n1)
    eps2 <- rep(NA_real_, n2)

    Xtilde1 <- matrix(NA_real_, nrow = n1, ncol = ncol(X1))
    Xtilde2 <- matrix(NA_real_, nrow = n2, ncol = ncol(X2))

    # within seg1: t = 2..n1
    if (n1 >= 2) {
      eps1[2:n1] <- u1[2:n1] - phi * u1[1:(n1 - 1)]
      Xtilde1[2:n1, ] <- X1[2:n1, , drop = FALSE] - phi * X1[1:(n1 - 1), , drop = FALSE]
    }

    # boundary for first point in seg2 (if exists)
    if (n2 >= 1 && n1 >= 1) {
      eps2[1] <- u2[1] - phi * u1[n1]
      Xtilde2[1, ] <- X2[1, , drop = FALSE] - phi * X1[n1, , drop = FALSE]
    }

    # within seg2: t = 2..n2
    if (n2 >= 2) {
      eps2[2:n2] <- u2[2:n2] - phi * u2[1:(n2 - 1)]
      Xtilde2[2:n2, ] <- X2[2:n2, , drop = FALSE] - phi * X2[1:(n2 - 1), , drop = FALSE]
    }

    list(eps1 = eps1, eps2 = eps2, Xtilde1 = Xtilde1, Xtilde2 = Xtilde2)
  }

  # -------- main loop --------
  p_n <- length(candidate_cps) - 2
  T_stats <- rep(-Inf, p_n)

  for (j in 1:p_n) {

    # Segment boundaries in "candidate index space"
    start <- candidate_cps[j] + 1
    cp    <- candidate_cps[j + 1]
    end   <- candidate_cps[j + 2]

    n_j  <- cp - start + 1
    n_j1 <- end - cp

    if (n_j < 1 || n_j1 < 1) {
      T_stats[j] <- -Inf
      next
    }

    # Map back to full-data indices (your original mapping)
    start_full <- 2 * start - 1
    cp_full    <- 2 * cp
    end_full   <- 2 * end

    idx1 <- start_full:cp_full
    idx2 <- (cp_full + 1):end_full

    if (length(idx1) < 2 || length(idx2) < 1) {
      T_stats[j] <- -Inf
      next
    }

    # Split masks
    in_O_seg1 <- idx1 %in% O_indices
    in_E_seg1 <- idx1 %in% E_indices
    in_O_seg2 <- idx2 %in% O_indices
    in_E_seg2 <- idx2 %in% E_indices

    # Extract full segments
    segment1Y <- Y[idx1]
    segment2Y <- Y[idx2]
    segment1X <- X[idx1, , drop = FALSE]
    segment2X <- X[idx2, , drop = FALSE]

    # Subsamples used to fit parameters (odd/even)
    y1O <- segment1Y[in_O_seg1]
    X1O <- segment1X[in_O_seg1, , drop = FALSE]
    y1E <- segment1Y[in_E_seg1]
    X1E <- segment1X[in_E_seg1, , drop = FALSE]

    # Need enough points to fit; if not, skip
    if (length(y1O) < 2 || length(y1E) < 2) {
      T_stats[j] <- -Inf
      next
    }

    # Fit on O and E (once each)
    fitO  <- trendARpsegfit_arimax(y1O, X1O, p = p)
    fitE  <- trendARpsegfit_arimax(y1E, X1E, p = p)
    coefO <- fitO$coef
    coefE <- fitE$coef

    beta_hatO <- safe_beta2(coefO)
    beta_hatE <- safe_beta2(coefE)
    phiO      <- safe_phi1(coefO, p)
    phiE      <- safe_phi1(coefE, p)

    # Residuals for FULL segments using SAME beta (as you intended)
    u1O_full <- as.numeric(segment1Y - segment1X %*% beta_hatO)
    u2O_full <- as.numeric(segment2Y - segment2X %*% beta_hatO)

    u1E_full <- as.numeric(segment1Y - segment1X %*% beta_hatE)
    u2E_full <- as.numeric(segment2Y - segment2X %*% beta_hatE)

    # Build innovations + whitened regressors for O and E
    wO <- make_eps_and_Xtilde_ar1(u1O_full, segment1X, u2O_full, segment2X, phiO)
    wE <- make_eps_and_Xtilde_ar1(u1E_full, segment1X, u2E_full, segment2X, phiE)

    # Now restrict to O/E indices AND drop NA rows (first obs of seg1 is NA)
    # Segment 1
    keep1O <- in_O_seg1 & !is.na(wO$eps1)
    keep1E <- in_E_seg1 & !is.na(wE$eps1)
    # Segment 2
    keep2O <- in_O_seg2 & !is.na(wO$eps2)
    keep2E <- in_E_seg2 & !is.na(wE$eps2)

    eps1O <- wO$eps1[keep1O]
    eps2O <- wO$eps2[keep2O]
    X1O_t <- wO$Xtilde1[keep1O, , drop = FALSE]
    X2O_t <- wO$Xtilde2[keep2O, , drop = FALSE]

    eps1E <- wE$eps1[keep1E]
    eps2E <- wE$eps2[keep2E]
    X1E_t <- wE$Xtilde1[keep1E, , drop = FALSE]
    X2E_t <- wE$Xtilde2[keep2E, , drop = FALSE]

    # Need non-empty for U-statistic
    if (length(eps1O) < 1 || length(eps2O) < 1 || length(eps1E) < 1 || length(eps2E) < 1) {
      T_stats[j] <- -Inf
      next
    }

    # Scores: s_t = Xtilde_t * eps_t  (no sigma scaling)
    S1O <- X1O_t * as.numeric(eps1O)
    S2O <- X2O_t * as.numeric(eps2O)

    S1E <- X1E_t * as.numeric(eps1E)
    S2E <- X2E_t * as.numeric(eps2E)

    # Compute L_j^E and L_j^O
    result_E <- compute_U_statistic_linear_trend(S1E, S2E)
    result_O <- compute_U_statistic_linear_trend(S1O, S2O)

    L_j_E <- result_E$U_stat
    L_j_O <- result_O$U_stat
    S1    <- result_O$S1
    S2    <- result_O$S2

    # Covariance estimate
    Sigma_hat <- moving_range_Omega_inv(S1, S2)

    if (rcond(Sigma_hat) < 1e-10) {
      Sigma_hat <- Sigma_hat + diag(ncol(X)) * 1e-6
    }

    scale_factor <- (n_j * n_j1) / (n_j + n_j1)
    T_stats[j] <- as.numeric(scale_factor * t(L_j_O) %*% solve(Sigma_hat) %*% L_j_E)
  }

  T_stats
}

#' Compute U-statistic for linear trend change-point
compute_U_statistic_linear_trend <- function(segment1, segment2, kernel_type = "moment") {
  U_stat <- colMeans(segment1, na.rm = TRUE) - colMeans(segment2, na.rm = TRUE)
  return(list(U_stat=U_stat, S1=segment1, S2=segment2))
}

trendARpsegfit_arimax <- function(y, X, p = p) {
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

  # helper: safely extract by name, fill missing with 0
  safe_extract_vec <- function(coefs, names_vec) {
    out <- numeric(length(names_vec))
    names(out) <- names_vec
    
    if (is.null(coefs)) return(out)
    
    present <- intersect(names_vec, names(coefs))
    out[present] <- coefs[present]
    
    out
  }

  beta_hat <- safe_extract_vec(coef_raw, beta_names)
  phi_hat  <- safe_extract_vec(coef_raw, ar_names)

  coef_all <- c(beta_hat, phi_hat)

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
generate_candidate_change_points <- function(Y, X, p=p, method = "WBS", p_n, eta_n) {
  n <- nrow(data)
    if (method == "PELT") {
      candidates = PELT.trendARp(Y, p=p, pen=0, minseglen=eta_n)
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
  AR_phi = 0.4,
  tau_manual = NULL
) {
  ## 1) Generate true change-points
  if (is.null(tau_manual)) {
    tau <- floor((1:K_n) * n / (K_n + 1)) + runif(K_n, -n^(1/4), n^(1/4))
    tau <- sort(pmax(2, pmin(n - 1, round(tau))))
    tau <- unique(tau)
  } else {
    tau <- sort(unique(as.integer(round(tau_manual))))
    if (length(tau) == 0) stop("tau_manual must contain at least one change-point.")
    if (any(!is.finite(tau))) stop("tau_manual must contain only finite values.")
    if (any(tau < 2 | tau > (n - 1))) {
      stop("tau_manual must be between 2 and n-1.")
    }
  }
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
                        infor <- Location_p[which.min(abs(Location_p-tau_star[t]))]
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




n <- 2000
K_n <- 10
sigma <- 1
AR_phi <- 0.4
error <- "gaussian"
#tau_manual <- c(91,135)
sim <- simulate_piecewise_linear_discontinuous(n=n, K_n=K_n, error = error, AR_phi = AR_phi, sigma=sigma)
Y <- sim$Y
X <- sim$X
block_size <- 8
alpha <- .2

tau_star <- sim$tau_true
K <- length(tau_star)
eta_n <- min(60, floor(n^(1/2)))
# tau_star

# Compare with moment kernel
result_moment <- robust_change_point_detection(
  Y,
  X=X, block_size=block_size,
  alpha = alpha, 
  candidate_method = "PELT",
  kernel_type = "moment", eta_n=eta_n,
)

th_L <- result_moment$candidate_points # candidate change-points from moment kernel A(O)
Location_e <- result_moment$change_points # detected change-points from moment kernel T

candidate_cps_refined <- generate_candidate_change_points(Y, X, method = "PELT", p=1,
                                                   p_n = floor(2 * n^(2/5)), eta_n = eta_n) # A(Z)

infor_tau<-sapply(1:K, function(t){             #informative points from A(O)
                infor<-th_L[which.min(abs(th_L-floor(tau_star[t])))]
                return(infor)                          
                  }              
          )

meas_index(Location_e, infor_tau, tau_star)

R_MOPS <- refine(candidate_cps_refined, Location_e)
Location_p <- R_MOPS$Location_p #Refined MOPS change-point Location T(Z)
R_infor_tau <- R_MOPS$infor_tau #Informative points A(Z)

meas_index(Location_p, R_infor_tau, tau_star)


fittrend = fit.trendARp(Y,Location_p,p=1,plot=T,add.ar=F,fit=T)#




# ---- Simulation study ----




n <- 2000
K_n <- 10
n_sim <- 200
sigma <- 1
AR_phi <- 0.4

error <- "gaussian"
block_size <- 8
alpha <- .2
eta_n <- min(60, floor(n^(1/2))) # set to 35 for n=1500

# pre-allocate
MOPS_results  <- matrix(NA_real_, nrow = n_sim, ncol = 5)
RMOPS_results <- matrix(NA_real_, nrow = n_sim, ncol = 5)

# define what "bad" means
is_bad_row <- function(x) {
  x <- as.numeric(x)
  any(!is.finite(x)) || all(x == 0) || identical(x, c(0, 0, 0, 0, Inf))
  # You can adjust: e.g., all(x[1:4] == 0) && is.infinite(x[5])
}

s <- 1
attempt <- 0

while (s <= n_sim) {
  attempt <- attempt + 1

  # (optional) safety break to avoid infinite loops
  # if (attempt > 10 * n_sim) stop("Too many invalid simulations; check settings.")

  # --- Simulate data ---
  sim <- simulate_piecewise_linear_discontinuous(n=n, K_n=K_n, error = error, AR_phi = AR_phi,  sigma=sigma)


  Y <- sim$Y
  X <- sim$X
  tau_star <- sim$tau_true
  K <- length(tau_star)


  # Wrap the whole iteration so failures also get skipped
  ok <- tryCatch({

    # --- Run MOPS ---
    result_moment <- robust_change_point_detection(
      Y,
      X=X, block_size=block_size,
      alpha = alpha, 
      candidate_method = "PELT",
      kernel_type = "moment", eta_n=eta_n,
    )


    th_L       <- result_moment$candidate_points
    Location_e <- result_moment$change_points

    infor_tau <- sapply(1:K, function(t){
      th_L[which.min(abs(th_L - floor(tau_star[t])))]
    })

    mops_row <- meas_index(Location_e, infor_tau, tau_star)

    # --- Run R-MOPS ---
    candidate_cps_refined <- generate_candidate_change_points(
      Y, X, method = "PELT", p=1,
      p_n = floor(2 * n^(2/5)),
      eta_n = eta_n
    )

    R_MOPS <- refine(candidate_cps_refined, Location_e)

    location_p  <- R_MOPS$Location_p
    R_infor_tau <- R_MOPS$infor_tau

    rmops_row <- meas_index(location_p, R_infor_tau, tau_star)

    # decide if we keep this iteration
    if (is_bad_row(mops_row) || is_bad_row(rmops_row)) {
      FALSE
    } else {
      MOPS_results[s, ]  <- as.numeric(mops_row)
      RMOPS_results[s, ] <- as.numeric(rmops_row)
      print(RMOPS_results[s,])
      TRUE
    }

  }, error = function(e) {
    FALSE
  })

  if (ok) {
    if (s %% 50 == 0) cat("Finished valid simulation", s, "(attempt", attempt, ")\n")
    s <- s + 1
  } else {
    # invalid attempt: do not increment s
    if (attempt %% 50 == 0) cat("Skipped invalid attempt", attempt, "\n")
  }
}

MOPS_df  <- as.data.frame(MOPS_results)
RMOPS_df <- as.data.frame(RMOPS_results)

summary_MOPS  <- colMeans(MOPS_df, na.rm = TRUE)
summary_RMOPS <- colMeans(RMOPS_df, na.rm = TRUE)

summary_table <- rbind(
  MOPS  = summary_MOPS,
  RMOPS = summary_RMOPS
)

print(summary_table)





# ---- grids you want to sweep ----
alpha_grid <- c(0.5)
sigma_grid <- c(.1, .5, 1)
block_grid <- c(2,3,4)
eta_grid <- c(5,6,7,8,9)

# ---- fixed settings ----
n <- 144
K_n <- 2
tau_manual <- c(91,135)
n_sim <- 1000
AR_phi <- 0.4
error <- "gaussian"
eta_n <- min(60, floor(n^(1/2)))

# "bad" result detector
is_bad_row <- function(x) {
  x <- as.numeric(x)
  any(!is.finite(x)) || all(x == 0) || identical(x, c(0, 0, 0, 0, Inf))
}

# one setting runner (your loop, parameterized)
run_one <- function(alpha, sigma, block_size, eta_n) {

  MOPS_results  <- matrix(NA_real_, nrow = n_sim, ncol = 5)
  RMOPS_results <- matrix(NA_real_, nrow = n_sim, ncol = 5)

  s <- 1
  attempt <- 0
  max_attempts <- 50 * n_sim

  while (s <= n_sim) {
    attempt <- attempt + 1
    if (attempt > max_attempts) {
      warning(sprintf(
        "Max attempts reached for alpha=%s sigma=%s block=%s (kept %s/%s).",
        alpha, sigma, block_size, s - 1, n_sim
      ))
      break
    }

    sim <- simulate_piecewise_linear_discontinuous(
      n = n, K_n = K_n, error = error, AR_phi = AR_phi, tau_manual=tau_manual,sigma = sigma
    )

    Y <- sim$Y
    X <- sim$X
    tau_star <- sim$tau_true
    K <- length(tau_star)

    ok <- tryCatch({

      # --- MOPS ---
      result_moment <- robust_change_point_detection(
        Y,
        X = X,
        block_size = block_size,
        alpha = alpha,
        candidate_method = "PELT",
        kernel_type = "moment",
        eta_n = eta_n
      )

      th_L       <- result_moment$candidate_points
      Location_e <- result_moment$change_points

      infor_tau <- sapply(seq_len(K), function(t) {
        th_L[which.min(abs(th_L - floor(tau_star[t])))]
      })

      mops_row <- meas_index(Location_e, infor_tau, tau_star)

      # --- R-MOPS ---
      candidate_cps_refined <- generate_candidate_change_points(
        Y, X, method = "PELT", p=1,
        p_n = floor(2 * n^(2/5)),
        eta_n = eta_n
      )

      R_MOPS <- refine(candidate_cps_refined, Location_e)
      location_p  <- R_MOPS$Location_p
      R_infor_tau <- R_MOPS$infor_tau

      rmops_row <- meas_index(location_p, R_infor_tau, tau_star)

      if (is_bad_row(mops_row) || is_bad_row(rmops_row)) {
        FALSE
      } else {
        MOPS_results[s, ]  <- as.numeric(mops_row)
        RMOPS_results[s, ] <- as.numeric(rmops_row)
        print(RMOPS_results[s,])
        TRUE
      }
    }, error = function(e) FALSE)

    if (ok) {
      s <- s + 1
    }
  }

  # summarize
  summary_MOPS  <- colMeans(as.data.frame(MOPS_results),  na.rm = TRUE)
  summary_RMOPS <- colMeans(as.data.frame(RMOPS_results), na.rm = TRUE)

  list(
    params = data.frame(alpha = alpha, sigma = sigma, block_size = block_size),
    summary_table = rbind(MOPS = summary_MOPS, RMOPS = summary_RMOPS),
    n_valid = sum(complete.cases(MOPS_results) & complete.cases(RMOPS_results)),
    attempts = attempt
  )
}

# ---- run all combinations ----
grid <- expand.grid(
  alpha = alpha_grid,
  sigma = sigma_grid,
  block_size = block_grid,
  eta_n = eta_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  cat(sprintf("Running %d/%d: alpha=%s sigma=%s block=%s eta_n=%s\n",
              i, nrow(grid), grid$alpha[i], grid$sigma[i], grid$block_size[i], grid$eta_n[i]))
  results[[i]] <- run_one(alpha = grid$alpha[i],
                          sigma = grid$sigma[i],
                          block_size = grid$block_size[i],
                          eta_n = grid$eta_n[i])
}

# ---- tidy output ----
out <- do.call(rbind, lapply(seq_along(results), function(i) {
  res <- results[[i]]
  if (is.null(res) || is.null(res$summary_table) || nrow(res$summary_table) == 0) return(NULL)

  tab <- res$summary_table
  eta_val <- if (!is.null(res$params$eta_n) && length(res$params$eta_n) == 1) {
    res$params$eta_n
  } else {
    grid$eta_n[i]
  }

  data.frame(
    alpha = res$params$alpha,
    sigma = res$params$sigma,
    block_size = res$params$block_size,
    eta_n = eta_val,
    method = rownames(tab),
    n_valid = res$n_valid,
    attempts = res$attempts,
    tab,
    row.names = NULL,
    check.names = FALSE
  )
}))

print(out)
saveRDS(out, "simulation_results_alpha_sigma_block_n144.rds")
