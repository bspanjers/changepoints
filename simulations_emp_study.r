library(strucchange)
library(zoo)
library(WeightedPortTest) # for the portmanteau test
source('./MethodCode/PELTtrendARp.R')
source('./MethodCode/MainCode.R')

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
  -3.3092150, 0.037958700,  0.1967182), nrow = 3, ncol = 3, byrow = TRUE) # setting B
coef_matrix <- NULL # setting A, 

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

      mops_row <- meas_index(Location_e, infor_tau, tau_star, K)

      # --- R-MOPS ---
      candidate_cps_refined <- generate_candidate_change_points(
        Y, X, method = "PELT", p=1,q=q,
        p_n = floor(2 * n^(2/5)),
        eta_n = eta_n
      )

      R_MOPS <- refine(candidate_cps_refined, Location_e, K, tau_star)
      location_p  <- R_MOPS$Location_p
      R_infor_tau <- R_MOPS$infor_tau

      rmops_row <- meas_index(location_p, R_infor_tau, tau_star, K)

      if (is_bad_row(mops_row) || is_bad_row(rmops_row)) {
        FALSE
      } else {
        MOPS_results[s, ]  <- as.numeric(mops_row)
        RMOPS_results[s, ] <- as.numeric(rmops_row)
        #print(RMOPS_results[s,])
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
saveRDS(out, "./Results/simulation_results_alpha_sigma2_block_n144_settingA.rds")
