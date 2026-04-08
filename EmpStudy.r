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
df_adj_yearly <- subset(Tanom_annual_df, Year >= 1880)
names(df_adj_yearly)[1] <- "t"
Y_adj_yearly <- df_adj_yearly$NOAA.adj
dates_adj_yearly = df_adj_yearly$t

# For adj: NASA / NOAA blocksize=3, Berk/HadCRU blocksize=2

#nasa 24: 1.28 25: 1.19
#noaa 24: 1.26 25: 1.14
#hadcrut 24: 1.17 25: 1.05
#Berkeley 24: 1.29 25: 1.19
# load raw yearly data
file_raw_yearly = './data/temperature_anomalies.RData'
load(file_raw_yearly)
df_raw_yearly <- subset(Tanom_annual_df, year >= 1880 & year < 2025)
names(df_raw_yearly)[1] <- "t"
Y_raw_yearly <- df_raw_yearly$HadCRUT
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
p <- 1
# Compare with moment kernel
result_moment <- robust_change_point_detection(
  Y,
  X, 
  p=p,
  block_size = 2,
  alpha = .5, 
  candidate_method = "PELT",
  kernel_type = "moment", eta_n=eta_n,
)
Location_e <- result_moment$change_points # detected change-points from moment kernel T
fittrend_e = fit.trendARp(Y,Location_e,p=p,plot=T,add.ar=F,fit=T,
                            dates=dates)#
K <- length(Location_e)
candidate_cps_refined <- generate_candidate_change_points(Y, X, p=p,method = "PELT", 
                                                   p_n = floor(2 * n^(2/5)), eta_n = eta_n) # A(Z)


R_MOPS <- refine(candidate_cps_refined, Location_e)
Location_p <- R_MOPS$Location_p #Refined MOPS change-point Location T(Z)

fittrend_p = fit.trendARp(Y, Location_p, p=p,plot=T,add.ar=F,fit=T,
                            dates=dates)#
LAG=20
resid <- Y - fittrend_p$fit_total
Weighted.Box.test(resid,lag=LAG,type="Ljung",fitdf=1)$p.value



# =========================
# Full code (single chunk)
# =========================

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
    if (model %in% c("NASA")) return(2L)
    if (model %in% c("NOAA")) return(2L)
    if (model %in% c("HadCRU")) return(2L)
    if (model %in% c("Berkeley")) return(2L)
  }
  if (dataset == "raw") {
    # For raw: HadCRU/NOAA/NASA blocksize=2, Berk blocksize=4
    if (model %in% c("NASA")) return(2L)
    if (model %in% c("NOAA")) return(2L)
    if (model %in% c("HadCRU")) return(2L)
    if (model %in% c("Berkeley")) return(2L)
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
    Y, X,p=p,
    method = "PELT",
    p_n = floor(2 * n^(2/5)),
    eta_n = eta_n
  )

  R_MOPS <- refine(candidate_cps_refined, Location_e)
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

# -------------------------
# Plot builder (raw or adj)
# -------------------------
make_cp_plots <- function(dataset = c("raw", "adj"),
                          df_raw_yearly, df_adj_yearly,
                          raw_models, adj_models,
                          eta_n = 9, alpha = 0.5, p = 1, LAG = 20) {

  dataset <- match.arg(dataset)

  if (dataset == "raw") {
    models <- raw_models
    dates_all <- df_raw_yearly$t
  } else {
    models <- adj_models
    dates_all <- df_adj_yearly$t
  }

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

    ggplot() +
      geom_line(aes(x = res$dates, y = res$Y), color = "blue") +
      geom_line(aes(x = res$dates, y = res$fittrend_p$fit_trend), color = "red") +
      { if (length(res$Location_p) > 0)
          geom_vline(xintercept = res$dates[res$Location_p], linetype = "dotted", color = "black")
        else NULL } +
      labs(
        title = paste0(
          toupper(dataset), " - ", model,
          " (block=", bs,
          ", Ljung Box p-value=", format_p(res$lb_pval),
          ", Shapiro p-value=", format_p(res$normality_pval), ")"
        ),
        x = "Year",
        y = "Temperature Anomaly"
      ) +
      theme_minimal()
  })

  plots
}

# -------------------------
# Build + arrange grids
# -------------------------
raw_plots <- make_cp_plots(
  dataset = "raw",
  df_raw_yearly = df_raw_yearly,
  df_adj_yearly = df_adj_yearly,
  raw_models = raw_models,
  adj_models = adj_models,
  eta_n = 9, alpha = 0.5, p = 1, LAG = 20
)

adj_plots <- make_cp_plots(
  dataset = "adj",
  df_raw_yearly = df_raw_yearly,
  df_adj_yearly = df_adj_yearly,
  raw_models = raw_models,
  adj_models = adj_models,
  eta_n = 9, alpha = 0.5, p = 1, LAG = 20
)

gridExtra::grid.arrange(grobs = raw_plots, ncol = 2, top = "Raw Data - Four Models")
gridExtra::grid.arrange(grobs = adj_plots, ncol = 2, top = "Adjusted Data - Four Models")
