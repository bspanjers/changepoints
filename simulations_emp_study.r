library(strucchange)
library(zoo)
library(WeightedPortTest) # for the portmanteau test
source('./MethodCode/PELTtrendARp.R')
source('./MethodCode/MainCode.R')

n <- 144
K_n <-2
sigma2 <- .01
AR_phi <- 0.4
error <- "gaussian"
tau_manual <-c(91,135)
coef_matrix <- matrix(c(
  -0.1005093, 0.006467981,  0.6531327,
  -0.9315647, 0.017054831, -0.1998985,
  -3.3092150, 0.037958700,  0.1967182), nrow = 3, ncol = 3, byrow = TRUE)
sim <- simulate_piecewise_linear_discontinuous(n=n, K_n=K_n, error = error, tau_manual=tau_manual, AR_phi = AR_phi, coef_matrix=coef_matrix, sigma2=sigma2)
Y <- sim$Y
X <- sim$X
block_size <- 2
alpha <- .2
q <- 0
tau_star <- sim$tau_true
K <- length(tau_star)
eta_n <- 9#9min(60, floor(n^(1/2)))


# tau_star

# Compare with moment kernel
result_moment <- robust_change_point_detection(
  Y,q = q,
  X=X, block_size=block_size,
  alpha = alpha, 
  candidate_method = "PELT",
  kernel_type = "moment", eta_n=eta_n,
)

th_L <- result_moment$candidate_points # candidate change-points from moment kernel A(O)
Location_e <- result_moment$change_points # detected change-points from moment kernel T

candidate_cps_refined <- generate_candidate_change_points(Y, X, method = "PELT", p=1, q=q,
                                                   p_n = floor(2 * n^(2/5)), eta_n = eta_n) # A(Z)

infor_tau<-sapply(1:K, function(t){             #informative points from A(O)
                infor<-th_L[which.min(abs(th_L-floor(tau_star[t])))]
                return(infor)                          
                  }              
          )

meas_index(Location_e, infor_tau, tau_star, K)

R_MOPS <- refine(candidate_cps_refined, Location_e, K, tau_star) # Refine to get R-MOPS change-points T(Z) and informative points A(Z)
Location_p <- R_MOPS$Location_p #Refined MOPS change-point Location T(Z)
R_infor_tau <- R_MOPS$infor_tau #Informative points A(Z)

meas_index(Location_p, R_infor_tau, tau_star, K)


fittrend = fit.trendARp(Y,Location_p,p=1,plot=T,add.ar=F,fit=T)#




# ---- Simulation study ----



K <- 2 #

n <- 2000
K_n <- 10
sigma2 <- .1
AR_phi <- 0.4
n_sim <- 100
tau_manual <-NULL#c(91,135)

error <- "gaussian"
block_size <- 5
alpha <- .2
eta_n <- min(60, floor(n^(1/2))) 

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
  sim <- simulate_piecewise_linear_discontinuous(n=n, K_n=K_n, error = error, tau_manual=tau_manual, AR_phi = AR_phi,  sigma2=sigma2)


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

    mops_row <- meas_index(Location_e, infor_tau, tau_star, K_n)

    # --- Run R-MOPS ---
    candidate_cps_refined <- generate_candidate_change_points(
      Y, X, method = "PELT", p=1,q=0,
      p_n = floor(2 * n^(2/5)),
      eta_n = eta_n
    )

    R_MOPS <- refine(candidate_cps_refined, Location_e, K_n, tau_star)

    location_p  <- R_MOPS$Location_p
    R_infor_tau <- R_MOPS$infor_tau

    rmops_row <- meas_index(location_p, R_infor_tau, tau_star, K_n)

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
sigma2_grid <- c(.01, .25, .5)
block_grid <- c(2,3,4)
eta_grid <- c(6,7,8,9)

# ---- fixed settings ----
n <- 144
K_n <- 2
tau_manual <- c(91,135)
coef_matrix <- matrix(c(
  -0.1005093, 0.006467981,  0.6531327,
  -0.9315647, 0.017054831, -0.1998985,
  -3.3092150, 0.037958700,  0.1967182), nrow = 3, ncol = 3, byrow = TRUE)
n_sim <- 1000
AR_phi <- 0.4
error <- "gaussian"

# "bad" result detector
is_bad_row <- function(x) {
  x <- as.numeric(x)
  any(!is.finite(x)) || all(x == 0) || identical(x, c(0, 0, 0, 0, Inf))
}

# one setting runner (your loop, parameterized)
run_one <- function(alpha, sigma2, block_size, eta_n) {

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
      n = n, K_n = K_n, error = error, AR_phi = AR_phi, coef_matrix=coef_matrix,tau_manual=tau_manual,sigma2 = sigma2
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
        Y, X, method = "PELT", p=1,q=q,
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
    params = data.frame(alpha = alpha, sigma2 = sigma2, block_size = block_size),
    summary_table = rbind(MOPS = summary_MOPS, RMOPS = summary_RMOPS),
    n_valid = sum(complete.cases(MOPS_results) & complete.cases(RMOPS_results)),
    attempts = attempt
  )
}

# ---- run all combinations ----
grid <- expand.grid(
  alpha = alpha_grid,
  sigma2 = sigma2_grid,
  block_size = block_grid,
  eta_n = eta_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  cat(sprintf("Running %d/%d: alpha=%s sigma2=%s block=%s eta_n=%s\n",
              i, nrow(grid), grid$alpha[i], grid$sigma2[i], grid$block_size[i], grid$eta_n[i]))
  results[[i]] <- run_one(alpha = grid$alpha[i],
                          sigma2 = grid$sigma2[i],
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
    sigma2 = res$params$sigma2,
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
saveRDS(out, "simulation_results_alpha_sigma2_block_n144.rds")
