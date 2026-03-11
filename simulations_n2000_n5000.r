library(strucchange)
library(zoo)
library(WeightedPortTest)
library(parallel)
library(foreach)
library(doParallel)

source('./MethodCode/PELTtrendARp.R')
source('./MethodCode/MainCode.R')

# ---- grids you want to sweep ----
alpha_grid <- c(.1, .2)
sigma2_grid <- c(.01, .1, .25, .5)
eta_grid <- c(34, 39, 44, 49, 54)

# ---- block_size grid depends on sigma2 ----
get_block_grid <- function(sigma2) {
  if (sigma2 %in% c(0.01, 0.1)) {
    return(1:5)
  } else if (sigma2 %in% c(0.25, 0.5)) {
    return(5:10)
  } else {
    return(1:5)  # default
  }
}

# ---- fixed settings ----
n <- 2000
K_n <- 10
n_sim <- 1000
AR_phi <- 0.4
error <- "gaussian"
tau_manual <- NULL

# settings needed for refine() and meas_index() to run, but not used in evaluation
K <- K_n # doesnt matter but needed for refine() to run
tau_star <- c(rep(0, K)) # dummy for refine() to run; not used in our evaluation

# "bad" result detector
is_bad_row <- function(x) {
  x <- as.numeric(x)
  any(!is.finite(x)) || all(x == 0) || identical(x, c(0, 0, 0, 0, Inf))
}

# Single simulation run (self-contained, no global dependencies)
run_single_sim <- function(sim_id, alpha, sigma2, block_size, eta_n, n, K_n, error, AR_phi, tau_manual) {
  
  max_attempts <- 50
  attempt <- 0
  
  while (attempt < max_attempts) {
    attempt <- attempt + 1
    
    result <- tryCatch({
      sim <- simulate_piecewise_linear_discontinuous(
        n = n, K_n = K_n, error = error, AR_phi = AR_phi, 
        tau_manual = tau_manual, sigma2 = sigma2
      )
      
      Y <- sim$Y
      X <- sim$X
      tau_star <- sim$tau_true
      K <- length(tau_star)
      
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
      
      th_L <- result_moment$candidate_points
      Location_e <- result_moment$change_points
      
      if (length(th_L) == 0) {
        return(NULL)
      }
      
      infor_tau <- sapply(seq_len(K), function(t) {
        th_L[which.min(abs(th_L - floor(tau_star[t])))]
      })
      
      # Pass K explicitly to meas_index
      mops_row <- meas_index(Location_e, infor_tau, tau_star, K)
      
      # --- R-MOPS ---
      candidate_cps_refined <- generate_candidate_change_points(
        Y, X, method = "PELT", p = 1, q=0,
        p_n = floor(2 * n^(2 / 5)),
        eta_n = eta_n
      )
      
      # Pass K and tau_star explicitly to refine
      R_MOPS <- refine(candidate_cps_refined, Location_e, K, tau_star)
      location_p <- R_MOPS$Location_p
      R_infor_tau <- R_MOPS$infor_tau
      
      # Pass K explicitly to meas_index
      rmops_row <- meas_index(location_p, R_infor_tau, tau_star, K)

      if (is_bad_row(mops_row) || is_bad_row(rmops_row)) {
        NULL
      } else {
        list(mops = as.numeric(mops_row), rmops = as.numeric(rmops_row))
      }
    }, error = function(e) {
      # Uncomment for debugging:
      # cat(sprintf("Error in sim %d: %s\n", sim_id, e$message))
      NULL
    })
    
    if (!is.null(result)) {
      return(result)
    }
  }
  
  # Return NA if all attempts failed
  return(list(mops = rep(NA_real_, 5), rmops = rep(NA_real_, 5)))
}

# One setting runner using foreach (more efficient)
run_one_parallel <- function(alpha, sigma2, block_size, eta_n, n_cores = NULL) {
  
  if (is.null(n_cores)) {
    n_cores <- max(1, detectCores() - 2)
  }
  
  cat(sprintf("  Using %d cores for %d simulations\n", n_cores, n_sim))
  
  # Register parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Ensure cleanup on exit
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
  })
  
  # Load packages on workers FIRST
  clusterEvalQ(cl, {
    library(strucchange)
    library(zoo)
    library(WeightedPortTest)
    library(forecast)
    source('/Users/barend/Documents/changepoints/MethodCode/PELTtrendARp.R')
    source('/Users/barend/Documents/changepoints/MethodCode/MainCode.R')
  })
  
  # Export the current function parameters to workers
  local_params <- list(
    alpha = alpha,
    sigma2 = sigma2,
    block_size = block_size,
    eta_n = eta_n,
    n = n,
    K_n = K_n,
    error = error,
    AR_phi = AR_phi,
    tau_manual = tau_manual
  )
  clusterExport(cl, c("local_params", "is_bad_row", "run_single_sim", "n_sim"), envir = environment())

  # Run using foreach
  results_list <- foreach(
    i = 1:n_sim,
    .errorhandling = "pass",
    .packages = c("strucchange", "zoo", "WeightedPortTest", "forecast")
  ) %dopar% {
    run_single_sim(
      sim_id = i,
      alpha = local_params$alpha,
      sigma2 = local_params$sigma2,
      block_size = local_params$block_size,
      eta_n = local_params$eta_n,
      n = local_params$n,
      K_n = local_params$K_n,
      error = local_params$error,
      AR_phi = local_params$AR_phi,
      tau_manual = local_params$tau_manual
    )
  }
  
  # Collect results
  MOPS_results <- matrix(NA_real_, nrow = n_sim, ncol = 5)
  RMOPS_results <- matrix(NA_real_, nrow = n_sim, ncol = 5)
  
  for (i in seq_along(results_list)) {
    res <- results_list[[i]]
    if (!is.null(res) && !inherits(res, "error") && 
        !is.null(res$mops) && !is.null(res$rmops)) {
      MOPS_results[i, ] <- res$mops
      RMOPS_results[i, ] <- res$rmops
    }
  }
  
  # Summarize
  n_valid <- sum(complete.cases(MOPS_results) & complete.cases(RMOPS_results))
  summary_MOPS <- colMeans(MOPS_results, na.rm = TRUE)
  summary_RMOPS <- colMeans(RMOPS_results, na.rm = TRUE)
  
  list(
    params = data.frame(alpha = alpha, sigma2 = sigma2, 
                        block_size = block_size, eta_n = eta_n),
    summary_table = rbind(MOPS = summary_MOPS, RMOPS = summary_RMOPS),
    n_valid = n_valid,
    attempts = n_sim
  )
}

# ---- run all combinations ----
# Build custom grid with sigma2-dependent block_size
grid <- do.call(rbind, lapply(sigma2_grid, function(s2) {
  block_sizes <- get_block_grid(s2)
  expand.grid(
    alpha = alpha_grid,
    sigma2 = s2,
    block_size = block_sizes,
    eta_n = eta_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}))

# Reset row names
rownames(grid) <- NULL

cat(sprintf("Total grid combinations: %d\n", nrow(grid)))
cat("Grid preview:\n")
print(head(grid, 20))

n_cores <- max(1, detectCores() - 2)
cat(sprintf("Detected %d cores, using %d for parallel execution\n", detectCores(), n_cores))

results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  cat(sprintf("Running %d/%d: alpha=%s sigma2=%s block=%s eta_n=%s\n",
              i, nrow(grid), grid$alpha[i], grid$sigma2[i], 
              grid$block_size[i], grid$eta_n[i]))
  
  start_time <- Sys.time()
  
  results[[i]] <- run_one_parallel(
    alpha = grid$alpha[i],
    sigma2 = grid$sigma2[i],
    block_size = grid$block_size[i],
    eta_n = grid$eta_n[i],
    n_cores = n_cores
  )
  
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("  Completed in %.2f minutes (n_valid = %d)\n", 
              elapsed, results[[i]]$n_valid))
}

# ---- tidy output ----
out <- do.call(rbind, lapply(seq_along(results), function(i) {
  res <- results[[i]]
  if (is.null(res) || is.null(res$summary_table) || 
      nrow(res$summary_table) == 0) return(NULL)
  
  tab <- res$summary_table
  
  data.frame(
    alpha = res$params$alpha,
    sigma2 = res$params$sigma2,
    block_size = res$params$block_size,
    eta_n = res$params$eta_n,
    method = rownames(tab),
    n_valid = res$n_valid,
    attempts = res$attempts,
    tab,
    row.names = NULL,
    check.names = FALSE
  )
}))

print(out)
saveRDS(out, "simulation_results_block_eta_n2000_21-2.rds")


# df must be a data.frame with columns:
# alpha, sigma2, block_size, eta_n, method, n_valid, attempts, and metrics in columns "1","2","3","4","5"
# where:
# 1=FDR, 2=TPR, 3=Pa, 4=Khat, 5=d
library(dplyr)
library(tidyr)

# out must have columns:
# alpha, sigma2, block_size, eta_n, method, n_valid, attempts, and metrics in columns "1","2","3","4","5"
# where:
# 1=FDR, 2=TPR, 3=Pa, 4=Khat, 5=d

# -----------------------
# 1) Keep only RMOPS + rename metrics
# -----------------------
df_rm <- out %>%
  mutate(method = toupper(method)) %>%
  filter(method == "RMOPS") %>%
  rename(
    FDR  = `1`,
    TPR  = `2`,
    Pa   = `3`,
    Khat = `4`,
    d    = `5`
  )

# -----------------------
# 2) Function: make LaTeX for one sigma2
#    (includes BOTH alphas as separate panels)
# -----------------------
make_latex_table_for_sigma <- function(df_rm, sigma2_value,
                                       caption = "Simulation results ($S=1000$, RMOPS)",
                                       label   = "tab:sim_rmops") {

  sub <- df_rm %>%
    filter(abs(sigma2 - sigma2_value) < 1e-12) %>%
    arrange(alpha, eta_n, block_size)

  blocks <- sort(unique(sub$block_size))
  metrics <- c("FDR","TPR","Pa","Khat","d")

  # Helper: format values
  fmt3 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", x))

  # Column format: eta + 5 * (#blocks)
  col_format <- paste0(
    "c",
    "@{\\hspace{15pt}}",
    paste(rep("c", 5 * length(blocks)), collapse = "")
  )

  # Header: block groups (shown as \ell = ...)
  group_header <- paste0(
    " & ",
    paste(
      vapply(blocks, function(b) sprintf("\\multicolumn{5}{c}{$\\ell = %d$}", b), character(1)),
      collapse = " & "
    ),
    " \\\\"
  )

  # cmidrules
  cmid <- paste(
    vapply(seq_along(blocks), function(i) {
      left  <- 2 + 5*(i-1)
      right <- 1 + 5*i
      sprintf("\\cmidrule(lr){%d-%d}", left, right)
    }, character(1)),
    collapse = "\n"
  )

  metric_header <- paste0(
    " & ",
    paste(rep("FDR & TPR & $P_a$ & $\\hat{K}$ & $d$", length(blocks)), collapse = " & "),
    " \\\\"
  )

  # Start LaTeX
  out_tex <- c(
    "\\begin{landscape}",
    "\\begin{table}[ht]",
    "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\setlength{\\tabcolsep}{4pt}",
    "\\begin{adjustbox}{max width=1.35\\textwidth}",
    sprintf("\\begin{tabular}{%s}", col_format),
    "\\toprule",
    paste0("$\\eta_n$", group_header),
    cmid,
    metric_header,
    "\\midrule",
    sprintf("\\multicolumn{%d}{c}{\\textbf{$\\sigma^2_{\\varepsilon} = %.2f$}} \\\\",
            1 + 5*length(blocks), sigma2_value),
    "\\midrule"
  )

  # Panels by alpha
  for (a in sort(unique(sub$alpha))) {

    out_tex <- c(out_tex,
      sprintf("\\multicolumn{%d}{l}{\\textbf{$\\alpha = %.1f$}} \\\\",
              1 + 5*length(blocks), a),
      "\\midrule"
    )

    wide <- sub %>%
      filter(abs(alpha - a) < 1e-12) %>%
      select(eta_n, block_size, all_of(metrics)) %>%
      pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") %>%
      mutate(metric = factor(metric, levels = metrics)) %>%
      pivot_wider(names_from = c(block_size, metric), values_from = value) %>%
      arrange(eta_n)

    # Ensure columns in order: block_size 1.. then metrics within each block
    ordered_cols <- unlist(lapply(blocks, function(b) paste0(b, "_", metrics)))
    names(wide) <- gsub("\\.", "_", names(wide))  # safety
    wide <- wide %>%
      select(eta_n, all_of(ordered_cols))

    # Write rows
    for (i in seq_len(nrow(wide))) {
      eta <- wide$eta_n[i]
      vals <- vapply(wide[i, -1], function(z) ifelse(is.na(z), "", sprintf("%.3f", z)), character(1))
      out_tex <- c(out_tex, paste0(eta, " & ", paste(vals, collapse = " & "), " \\\\"))
    }

    out_tex <- c(out_tex, "\\midrule")
  }

  # End LaTeX
  out_tex <- c(out_tex,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{adjustbox}",
    "\\end{table}",
    "\\end{landscape}"
  )

  paste(out_tex, collapse = "\n")
}

# -----------------------
# 3) Create one .tex per sigma2
# -----------------------
sigmas <- sort(unique(df_rm$sigma2))

for (s2 in sigmas) {
  tex <- make_latex_table_for_sigma(
    df_rm,
    sigma2_value = s2,
    caption = "Simulation results ($S=1000$, RMOPS)",
    label   = paste0("tab:sim_rmops_sigma", gsub("\\.", "", sprintf("%.2f", s2)))
  )
  outpath <- paste0("./sim_rmops_sigma", gsub("\\.", "", sprintf("%.2f", s2)), ".tex")
  writeLines(tex, outpath)
  message("✅ wrote ", outpath)
}