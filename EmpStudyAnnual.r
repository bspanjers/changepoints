library(strucchange)
library(zoo)
library(WeightedPortTest) # for the portmanteau test

source('./MethodCode/PELTtrendARp.R')
source('./MethodCode/MainCode.R')

K<-2 # doesnt matter but needed for refine() to run
tau_star <- c(1930, 1980) # dummy for refine() to run; not used in our evaluation

# load adj yearly data
file_adj_yearly = './data/temperature_anomalies_adj.RData'
load(file_adj_yearly)
df_adj_yearly <- subset(Tanom_annual_df, Year >= 1881)
names(df_adj_yearly)[1] <- "t"
dates_adj_yearly = df_adj_yearly$t
rownames(dates_adj_yearly) <- NULL

# For adj: NASA / NOAA blocksize=3, Berk/HadCRU blocksize=2

# load raw yearly data
file_raw_yearly = './data/temperature_anomalies.RData'
load(file_raw_yearly)
df_raw_yearly <- subset(Tanom_annual_df, year >= 1881 & year <= 2025)
names(df_raw_yearly)[1] <- "t"
dates_raw_yearly = df_raw_yearly$t
rownames(df_raw_yearly) <- NULL
#nasa 24: 1.28 25: 1.19
#noaa 24: 1.26 25: 1.14
#hadcrut 24: 1.17 25: 1.05
#Berkeley 24: 1.29 25: 1.19

# Add 2024 values
row_2024 <- which(df_raw_yearly$t == 2024)
if (length(row_2024) > 0) {
  df_raw_yearly$NASA[row_2024] <- 1.28
  df_raw_yearly$NOAA[row_2024] <- 1.26
  df_raw_yearly$HadCRUT[row_2024] <- 1.17
  df_raw_yearly$Berkeley[row_2024] <- 1.29
} else {
  # If 2024 row doesn't exist, add it
  new_row <- data.frame(
    t = 2024,
    NASA = 1.28,
    `Japan Met` = NA,
    NOAA = 1.26,
    HadCRUT = 1.17,
    Berkeley = 1.29,
    check.names = FALSE
  )
  new_row <- new_row[, names(df_raw_yearly)]
  df_raw_yearly <- rbind(df_raw_yearly, new_row)
}

# Add 2025 values
#row_2025 <- which(df_raw_yearly$t == 2025)
#if (length(row_2025) > 0) {
#  df_raw_yearly$NASA[row_2025] <- 1.19
#  df_raw_yearly$NOAA[row_2025] <- 1.14
#  df_raw_yearly$HadCRUT[row_2025] <- 1.05
#  df_raw_yearly$Berkeley[row_2025] <- 1.19
#} else {
#  # If 2025 row doesn't exist, add it
#  new_row <- data.frame(
#    t = 2025,
#    NASA = 1.19,
#    `Japan Met` = NA,
#    NOAA = 1.14,
#    HadCRUT = 1.05,
#    Berkeley = 1.19,
#    check.names = FALSE
#  )
#  new_row <- new_row[, names(df_raw_yearly)]
#  df_raw_yearly <- rbind(df_raw_yearly, new_row)
#}

# Reset row names after adding rows
rownames(df_raw_yearly) <- NULL

library(ggplot2)
library(gridExtra)

# -------------------------
# Model lists
# -------------------------
raw_models <- list(
  NASA     = df_raw_yearly$NASA,
  NOAA     = df_raw_yearly$NOAA,
  HadCRU   = df_raw_yearly$HadCRUT,   # raw column name in your data
  Berkeley = df_raw_yearly$Berkeley
)

adj_models <- list(
  NASA     = df_adj_yearly$NASA.adj,
  NOAA     = df_adj_yearly$NOAA.adj,
  HadCRU   = df_adj_yearly$HadCRU.adj,
  Berkeley = df_adj_yearly$Berk.adj
)

# -------------------------
# Helpers
# -------------------------
trim_even <- function(Y, dates) {
  m <- min(length(Y), length(dates))
  Y <- Y[1:m]
  dates <- dates[1:m]
  if (m %% 2 == 1) {
    Y <- Y[-m]
    dates <- dates[-m]
  }
  list(Y = Y, dates = dates)
}

sanitize_cps <- function(cps, n) {
  if (is.null(cps)) return(integer(0))
  if (is.list(cps)) cps <- unlist(cps, use.names = FALSE)
  cps <- as.integer(cps)
  cps <- cps[!is.na(cps)]
  cps <- cps[cps >= 1 & cps <= n]
  sort(unique(cps))
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  formatC(p, format = "f", digits = 3)
}

# block sizes you specified
get_block_size <- function(dataset, model) {
  if (dataset == "adj") {
    # For adj: NASA/NOAA blocksize=3/4, Berk/HadCRU blocksize=2
    if (model %in% c("NASA")) return(4L)
    if (model %in% c("NOAA")) return(4L)
    if (model %in% c("HadCRU")) return(4L)
    if (model %in% c("Berkeley")) return(4L)
  }
  if (dataset == "raw") {
    # For raw: HadCRU/NOAA/NASA blocksize=2, Berk blocksize=4
    if (model %in% c("NASA")) return(4L)
    if (model %in% c("NOAA")) return(4L)
    if (model %in% c("HadCRU")) return(4L)
    if (model %in% c("Berkeley")) return(4L)
  }
  stop("Unknown dataset/model: ", dataset, " / ", model)
}

# -------------------------
# Core runner for one series
# -------------------------
run_one_series <- function(Y, dates, block_size,
                           eta_n = 9, alpha = 0.5,
                           candidate_method = "PELT",
                           kernel_type = "moment",
                           p = 1,
                           LAG = 20,
                           plot_fit = FALSE) {

  te <- trim_even(Y, dates)
  Y <- te$Y
  dates <- te$dates

  n <- length(Y)
  X <- cbind(
    Intercept = rep(1, n),
    X_year = (1:n) / 12
  )

  # Moment kernel change-points
  result_moment <- robust_change_point_detection(
    Y, X,
    p=p,
    block_size = block_size,
    alpha = alpha,
    candidate_method = candidate_method,
    kernel_type = kernel_type,
    eta_n = eta_n
  )
  Location_e <- sanitize_cps(result_moment$change_points, n)

  # Candidate cps + refine
  candidate_cps_refined <- generate_candidate_change_points(
    Y, X,p=p,q=0,
    method = "PELT",
    p_n = floor(2 * n^(2/5)),
    eta_n = eta_n
  )

  R_MOPS <- refine(candidate_cps_refined, Location_e, K, tau_star)
  Location_p <- sanitize_cps(R_MOPS$Location_p, n)

  # Fit trend with refined cps
  fittrend_p <- fit.trendARp(
    Y, Location_p,
    p = p, plot = plot_fit, add.ar = FALSE, fit = TRUE, dates = dates
  )

  # Ljung-Box p-value on residuals (using your formula)
  resid <- Y - fittrend_p$fit_total
  lb_pval <- Weighted.Box.test(resid, lag = LAG, type = "Ljung", fitdf = 1)$p.value
  resid_finite <- resid[is.finite(resid)]
  normality_pval <- tryCatch({
    if (length(resid_finite) < 3 || length(resid_finite) > 5000) {
      NA_real_
    } else {
      shapiro.test(resid_finite)$p.value
    }
  }, error = function(e) NA_real_)

  list(
    Y = Y, dates = dates, X = X,
    Location_e = Location_e,
    Location_p = Location_p,
    fittrend_p = fittrend_p,
    lb_pval = lb_pval,
    normality_pval = normality_pval
  )
}
# ...existing code...

# -------------------------
# Plot builder (raw or adj) - Clean in-plot legend box
# -------------------------
make_cp_plots <- function(dataset = c("raw", "adj"),
                          df_raw_yearly, df_adj_yearly,
                          raw_models, adj_models,
                          eta_n = 9, alpha = 0.5, p = 1, LAG = 20,
                          return_results = FALSE) {

  dataset <- match.arg(dataset)

  if (dataset == "raw") {
    models <- raw_models
    dates_all <- df_raw_yearly$t
  } else {
    models <- adj_models
    dates_all <- df_adj_yearly$t
  }

  # Subtitle labels
  subtitle_labels <- c("(a) NASA", "(b) NOAA", "(c) HadCRU", "(d) Berkeley")
  names(subtitle_labels) <- c("NASA", "NOAA", "HadCRU", "Berkeley")

  # Store results for table
  results_list <- list()

  # First pass: compute shared axis limits
  all_Y <- c()
  all_fits <- c()
  all_dates <- c()
  
  for (model in names(models)) {
    Y <- models[[model]]
    bs <- get_block_size(dataset, model)
    res <- run_one_series(
      Y = Y,
      dates = dates_all,
      block_size = bs,
      eta_n = eta_n,
      alpha = alpha,
      p = p,
      LAG = LAG,
      plot_fit = FALSE
    )
    all_Y <- c(all_Y, res$Y)
    all_fits <- c(all_fits, res$fittrend_p$fit_trend)
    all_dates <- c(all_dates, res$dates)
  }
  
  shared_x_range <- range(all_dates, na.rm = TRUE)
  shared_y_range <- range(c(all_Y, all_fits), na.rm = TRUE)
  # Add small padding to y-axis
  y_pad <- 0.05 * diff(shared_y_range)
  shared_y_range <- c(shared_y_range[1] - y_pad, shared_y_range[2] + y_pad)

  plots <- lapply(names(models), function(model) {
    Y <- models[[model]]
    bs <- get_block_size(dataset, model)

    res <- run_one_series(
      Y = Y,
      dates = dates_all,
      block_size = bs,
      eta_n = eta_n,
      alpha = alpha,
      p = p,
      LAG = LAG,
      plot_fit = FALSE
    )

    # Store results
    results_list[[model]] <<- list(
      model = model,
      lb_pval = res$lb_pval,
      shapiro_pval = res$normality_pval
    )

    # Create data frames for ggplot
    df_data <- data.frame(dates = res$dates, Y = res$Y)
    df_fit <- data.frame(dates = res$dates, Y = res$fittrend_p$fit_trend)

    # Legend box position using SHARED ranges
    x_range <- shared_x_range
    y_range <- shared_y_range
    
    # Legend box position (top-left corner)
    box_x_min <- x_range[1] + 0.01 * diff(x_range)
    box_x_max <- x_range[1] + 0.22 * diff(x_range)
    box_y_max <- y_range[2] - 0.01 * diff(y_range)
    box_y_min <- y_range[2] - 0.18 * diff(y_range)
    
    # Line sample positions
    line_x_start <- box_x_min + 0.01 * diff(x_range)
    line_x_end <- box_x_min + 0.04 * diff(x_range)
    text_x <- line_x_end + 0.01 * diff(x_range)
    
    line_y1 <- box_y_max - 0.04 * diff(y_range)
    line_y2 <- box_y_max - 0.08 * diff(y_range)
    line_y3 <- box_y_max - 0.12 * diff(y_range)

    ggplot() +
      # Data line (grey, solid)
      geom_line(data = df_data, aes(x = dates, y = Y), 
                color = "grey50", size = 0.5) +
      # Fitted trend line (black, solid, thicker)
      geom_line(data = df_fit, aes(x = dates, y = Y), 
                color = "black", size = 0.8) +
      # Changepoint lines (black, dashed)
      { if (length(res$Location_p) > 0)
          geom_vline(xintercept = res$dates[res$Location_p],
                     linetype = "dashed", color = "black", size = 0.6)
        else NULL } +
      
      # ---- LEGEND BOX ----
      annotate("rect", 
               xmin = box_x_min, xmax = box_x_max, 
               ymin = box_y_min, ymax = box_y_max,
               fill = "white", color = "grey50", size = 0.15) +
      
      # Line 1: Data (grey solid)
      annotate("segment", 
               x = line_x_start, xend = line_x_end, 
               y = line_y1, yend = line_y1,
               color = "grey50", size = 0.5) +
      annotate("text", x = text_x, y = line_y1, 
               label = "Data", hjust = 0, size = 5, color = "grey30") +
      
      # Line 2: Fitted Trend (black solid)
      annotate("segment", 
               x = line_x_start, xend = line_x_end, 
               y = line_y2, yend = line_y2,
               color = "black", size = 0.8) +
      annotate("text", x = text_x, y = line_y2, 
               label = "Fitted Trend", hjust = 0, size = 5, color = "grey30") +
      
      # Line 3: Changepoints (black dashed)
      annotate("segment", 
               x = line_x_start, xend = line_x_end, 
               y = line_y3, yend = line_y3,
               color = "black", size = 0.6, linetype = "dashed") +
      annotate("text", x = text_x, y = line_y3, 
               label = "Changepoints", hjust = 0, size = 5, color = "grey30") +
      
      # Shared axis limits
      coord_cartesian(xlim = shared_x_range, ylim = shared_y_range) +
      
      # Labels - subtitle only, CENTERED and LARGER
      labs(
        title = subtitle_labels[model],
        x = "Year",
        y = "Temperature Anomaly"
      ) +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#e5e5e5"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "grey30", size = 9),
        axis.title = element_text(color = "grey30", size = 10),
        plot.title = element_text(color = "grey20", size = 14, face = "bold", hjust = 0.5)
      )
  })

  # Return both plots and results
  list(plots = plots, results = results_list)
}


# ...existing code...

# -------------------------
# Build + arrange grids
raw_output <- make_cp_plots(
  dataset = "raw",
  df_raw_yearly = df_raw_yearly,
  df_adj_yearly = df_adj_yearly,
  raw_models = raw_models,
  adj_models = adj_models,
  eta_n = 8, alpha = 0.5, p = 1, LAG = 20
)

adj_output <- make_cp_plots(
  dataset = "adj",
  df_raw_yearly = df_raw_yearly,
  df_adj_yearly = df_adj_yearly,
  raw_models = raw_models,
  adj_models = adj_models,
  eta_n = 8, alpha = 0.5, p = 1, LAG = 20
)

raw_plots <- raw_output$plots
adj_plots <- adj_output$plots

# For 2x2 grid:
# - Top row (plots 1,2): no x-axis label
# - Bottom row (plots 3,4): x-axis label "Year"
# - Left column (plots 1,3): y-axis label "Temperature Anomaly"
# - Right column (plots 2,4): no y-axis label

# Function to adjust labels based on position
adjust_plot_labels <- function(plots) {
  # Plot 1: top-left (y-label, no x-label)
  plots[[1]] <- plots[[1]] + labs(x = NULL, y = "Temperature Anomaly")
  # Plot 2: top-right (no labels)
  plots[[2]] <- plots[[2]] + labs(x = NULL, y = NULL)
  # Plot 3: bottom-left (both labels)
  plots[[3]] <- plots[[3]] + labs(x = "Year", y = "Temperature Anomaly")
  # Plot 4: bottom-right (x-label, no y-label)
  plots[[4]] <- plots[[4]] + labs(x = "Year", y = NULL)
  
  plots
}

raw_plots_clean <- adjust_plot_labels(raw_plots)
adj_plots_clean <- adjust_plot_labels(adj_plots)

# Arrange grids
gridExtra::grid.arrange(grobs = raw_plots_clean, ncol = 2)
gridExtra::grid.arrange(grobs = adj_plots_clean, ncol = 2)

# -------------------------
# P-value tables with Holm-Bonferroni adjustment
# -------------------------
# Raw data
raw_lb <- sapply(raw_output$results, function(x) x$lb_pval)
raw_shapiro <- sapply(raw_output$results, function(x) x$shapiro_pval)

raw_pval_table <- data.frame(
  Model = names(raw_output$results),
  `Ljung-Box (raw)` = sapply(raw_lb, format_p),
  `Ljung-Box (Holm)` = sapply(p.adjust(raw_lb, method = "holm"), format_p),
  `Shapiro (raw)` = sapply(raw_shapiro, format_p),
  `Shapiro (Holm)` = sapply(p.adjust(raw_shapiro, method = "holm"), format_p),
  check.names = FALSE
)

# Adjusted data
adj_lb <- sapply(adj_output$results, function(x) x$lb_pval)
adj_shapiro <- sapply(adj_output$results, function(x) x$shapiro_pval)

adj_pval_table <- data.frame(
  Model = names(adj_output$results),
  `Ljung-Box (raw)` = sapply(adj_lb, format_p),
  `Ljung-Box (Holm)` = sapply(p.adjust(adj_lb, method = "holm"), format_p),
  `Shapiro (raw)` = sapply(adj_shapiro, format_p),
  `Shapiro (Holm)` = sapply(p.adjust(adj_shapiro, method = "holm"), format_p),
  check.names = FALSE
)

# Print tables
cat("\n========================================\n")
cat("     RAW DATA - Residual Diagnostics\n")
cat("========================================\n")
print(raw_pval_table, row.names = FALSE)

cat("\n========================================\n")
cat("   ADJUSTED DATA - Residual Diagnostics\n")
cat("========================================\n")
print(adj_pval_table, row.names = FALSE)