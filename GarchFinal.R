## ------------- ONLY VaR and VaR Forecasting -----------


# Load required libraries
library(quantmod)
library(moments)
library(tseries)
library(FinTS)
library(rugarch)
library(rmgarch)
library(xts)
library(zoo)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

### -------------- GARCH(1,1) VaR estimation ----------------------
# Load packages
# Set date range and symbols
start_date <- as.Date("2019-04-30")
end_date <- Sys.Date()
symbols <- c("^GSPC", "^NDX", "^GDAXI", "^N225")
names(symbols) <- c("SP500", "NASDAQ100", "DAX", "NIKKEI")

# Download data
getSymbols(symbols, src = "yahoo", from = start_date, to = end_date)

# Prepare closing prices and log-returns
prices <- na.omit(merge(Cl(GSPC), Cl(NDX), Cl(GDAXI), Cl(N225)))
colnames(prices) <- names(symbols)
returns <- na.omit(diff(log(prices)))

# Train/test split (70/30)
n <- nrow(returns)
split_index <- floor(0.7 * n)
train <- returns[1:split_index, ]

# GARCH(1,1) specification
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "norm"
)

# Set up 2x2 plotting layout
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# Loop over each asset
for (i in 1:ncol(train)) {
  asset_name <- colnames(train)[i]
  ret_asset_xts <- train[, i]
  ret_asset <- as.numeric(ret_asset_xts)
  asset_dates <- index(ret_asset_xts)
  
  # Fit GARCH(1,1)
  fit <- ugarchfit(spec = garch_spec, data = ret_asset_xts, solver = "solnp")
  sigma_fit <- sigma(fit)
  
  # Compute 5% VaR
  VaR_5 <- qnorm(0.05) * sigma_fit
  violations <- ret_asset < VaR_5
  violation_pct <- mean(violations) * 100
  
  # Plot: returns, VaR, violations
  plot(asset_dates, ret_asset, type = "l", col = "gray",
       main = paste("GARCH(1,1) VaR 5% -", asset_name, "\nViolations:", round(violation_pct, 2), "%"),
       ylab = "Log-return", xlab = "")
  lines(asset_dates, VaR_5, col = "red")
  points(asset_dates[violations], ret_asset[violations], col = "blue", pch = 19)
  legend("bottomleft", legend = c("VaR 5%", "Violations"),
         col = c("red", "blue"), lty = c(1, NA), pch = c(NA, 19), bty = "n")
}
# ---- GARCH(1,1) portfolio VaR ----
# Define portfolio weights (equal weights or customize here)
w <- rep(1 / ncol(train), ncol(train))  # Equal-weighted

# Fit univariate GARCH models for each asset
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "norm"
)

garch_fits <- lapply(1:ncol(train), function(i) {
  ugarchfit(spec, data = train[, i], solver = "solnp")
})

# Extract conditional standard deviations
sigma_mat <- do.call(cbind, lapply(garch_fits, sigma))  # T x N matrix
colnames(sigma_mat) <- colnames(train)
rownames(sigma_mat) <- index(train)

# Approximate portfolio variance (ignoring correlations)
port_sd <- rowSums((sigma_mat * w)^2)^0.5
VaR_port <- qnorm(0.05) * port_sd

# Standardized residuals Q-Q plot for each asset
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
for (i in 1:ncol(train)) {
  fit <- ugarchfit(spec = garch_spec, data = train[, i], solver = "solnp")
  std_resid <- residuals(fit, standardize = TRUE)
  qqnorm(std_resid, main = paste("Q-Q Plot:", colnames(train)[i]),
         ylab = "Sample Quantiles", xlab = "Theoretical Quantiles")
  qqline(std_resid, col = "red")
}


# Compute actual portfolio returns
port_ret <- rowSums(train * matrix(rep(w, each = nrow(train)), ncol = ncol(train)))

# Create time series
VaR_xts <- xts(VaR_port, order.by = index(train))
port_ret_xts <- xts(port_ret, order.by = index(train))
violations <- port_ret_xts < VaR_xts
violation_pct <- mean(violations) * 100

# ---- Plot portfolio VaR ----
par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))
plot(port_ret_xts, type = "l", col = "gray",
     main = paste("GARCH(1,1) VaR 5% - Portfolio\nViolations:", round(violation_pct, 2), "%"),
     ylab = "Log-return", xlab = "")
lines(VaR_xts, col = "red")
points(index(port_ret_xts)[violations], port_ret_xts[violations], col = "blue", pch = 19)
legend("bottomleft", legend = c("VaR 5%", "Violations"),
       col = c("red", "blue"), lty = c(1, NA), pch = c(NA, 19), bty = "n")



### -------------- DCC-GARCH Var  ----------------------
# ---------- 1. VaR per ciascun indice con DCC-GARCH ----------
# GARCH specification per ogni serie
uspec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "norm"
)

# Specifica multivariata DCC
mspec <- dccspec(
  uspec = multispec(replicate(ncol(train), uspec)),
  dccOrder = c(1,1),
  distribution = "mvnorm"
)

# Fit DCC-GARCH multivariato
fit_dcc <- dccfit(mspec, data = train)

# Estrai le deviazioni standard condizionate
sigma_matrix <- sigma(fit_dcc)
returns_matrix <- train

# Plot 2x2: VaR singoli
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
for (i in 1:ncol(returns_matrix)) {
  asset_name <- colnames(returns_matrix)[i]
  ret_asset <- as.numeric(returns_matrix[, i])
  asset_dates <- index(returns_matrix)
  sigma_asset <- as.numeric(sigma_matrix[, i])
  
  VaR_5 <- qnorm(0.05) * sigma_asset
  violations <- ret_asset < VaR_5
  violation_pct <- mean(violations) * 100
  
  plot(asset_dates, ret_asset, type = "l", col = "gray",
       main = paste("DCC-GARCH VaR 5% -", asset_name, "\nViolations:", round(violation_pct, 2), "%"),
       ylab = "Log-return", xlab = "")
  lines(asset_dates, VaR_5, col = "red")
  points(asset_dates[violations], ret_asset[violations], col = "blue", pch = 19)
  legend("bottomleft", legend = c("VaR 5%", "Violations"),
         col = c("red", "blue"), lty = c(1, NA), pch = c(NA, 19), bty = "n")
}

# ---------- 2. VaR portafoglio pesato con DCC-GARCH ----------
# Estrai le matrici di covarianza
cov_matrices <- rcov(fit_dcc)
T <- dim(cov_matrices)[3]

# Pesi del portafoglio (equal weight)
w <- rep(1 / ncol(train), ncol(train))

# Inizializza vettori
ret_port <- VaR_port <- rep(NA, T)

for (t in 1:T) {
  Sigma_t <- cov_matrices[,,t]
  r_t <- as.numeric(train[t, ])
  ret_port[t] <- sum(w * r_t)
  VaR_port[t] <- qnorm(0.05) * sqrt(t(w) %*% Sigma_t %*% w)
}

# Serie temporali xts
dates <- index(train)
ret_port_xts <- xts(ret_port, order.by = dates)
VaR_port_xts <- xts(VaR_port, order.by = dates)
violations <- ret_port_xts < VaR_port_xts
violation_pct <- mean(violations) * 100

# Plot VaR portafoglio
par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))
plot(ret_port_xts, type = "l", col = "gray",
     main = paste("DCC-GARCH VaR 5% - Portfolio\nViolations:", round(violation_pct, 2), "%"),
     ylab = "Log-return", xlab = "")
lines(VaR_port_xts, col = "red")
points(index(ret_port_xts)[violations], ret_port_xts[violations], col = "blue", pch = 19)
legend("bottomleft", legend = c("VaR 5%", "Violations"),
       col = c("red", "blue"), lty = c(1, NA), pch = c(NA, 19), bty = "n")


## --------------- Portfolio rolling window --------------
# Parametri
window_size <- 300
test_start <- split_index + 1
n_forecast <- nrow(test)
w <- rep(1 / ncol(returns), ncol(returns))  # equal-weighted

# Inizializza
VaR_forecast <- rep(NA, n_forecast)
ret_port_forecast <- rep(NA, n_forecast)

# DCC-GARCH rolling 1-step-ahead forecast
for (i in 1:n_forecast) {
  start_idx <- test_start + i - window_size - 1
  end_idx <- start_idx + window_size - 1
  window_data <- returns[start_idx:end_idx, ]
  
  # Specifica
  uspec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
    distribution.model = "norm"
  )
  mspec <- dccspec(
    uspec = multispec(replicate(ncol(window_data), uspec)),
    dccOrder = c(1,1),
    distribution = "mvnorm"
  )
  
  # Fit e forecast
  fit_roll <- dccfit(mspec, data = window_data, fit.control = list(eval.se = FALSE), solver = "solnp")
  fcast <- dccforecast(fit_roll, n.ahead = 1)
  Sigma_1 <- rcov(fcast)[[1]][,,1]
  
  # VaR portafoglio 1-step ahead
  sd_port <- sqrt(t(w) %*% Sigma_1 %*% w)
  VaR_forecast[i] <- qnorm(0.05) * sd_port
  
  # Rendimento reale del portafoglio (test)
  ret_port_forecast[i] <- sum(w * as.numeric(test[i, ]))
}

# Serie xts
dates_test <- index(test)
ret_xts <- xts(ret_port_forecast, order.by = dates_test)
VaR_xts <- xts(VaR_forecast, order.by = dates_test)
violations <- ret_xts < VaR_xts
violation_pct <- mean(violations) * 100

# ---- Plot forecasted VaR ----
par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))
plot(ret_xts, type = "l", col = "gray",
     main = paste("DCC-GARCH Forecast VaR 5% - Portfolio\nRolling window: 300 | Violations:", round(violation_pct, 2), "%"),
     ylab = "Log-return", xlab = "")
lines(VaR_xts, col = "red")
points(index(ret_xts)[violations], ret_xts[violations], col = "blue", pch = 19)
legend("bottomleft", legend = c("VaR 5% Forecast", "Violations"),
       col = c("red", "blue"), lty = c(1, NA), pch = c(NA, 19), bty = "n")


## rolling forecast with GARCH(1,1) univariate
# Portfolio weights
w <- rep(1 / ncol(returns), ncol(returns))  # equally weighted

# Parameters
window_size <- 300
test_start <- split_index + 1
n_forecast <- nrow(test)

# Prepare storage
VaR_garch <- ret_port <- rep(NA, n_forecast)

# Define univariate GARCH(1,1) specification
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "norm"
)

# Rolling forecast
for (i in 1:n_forecast) {
  start_idx <- test_start + i - window_size - 1
  end_idx <- start_idx + window_size - 1
  window_data <- returns[start_idx:end_idx, ]
  
  sigma_next <- rep(NA, ncol(window_data))  # forecasted volatilities
  
  for (j in 1:ncol(window_data)) {
    fit <- tryCatch({
      ugarchfit(spec = garch_spec, data = window_data[, j],
                solver = "solnp", fit.control = list(scale = 1))
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      fcast <- ugarchforecast(fit, n.ahead = 1)
      sigma_next[j] <- sigma(fcast)[1]
    } else {
      sigma_next[j] <- NA
    }
  }
  
  if (all(!is.na(sigma_next))) {
    # Compute portfolio VaR assuming zero correlation
    port_sd <- sqrt(sum((w * sigma_next)^2))
    VaR_garch[i] <- qnorm(0.05) * port_sd
  } else {
    VaR_garch[i] <- NA
  }
  
  # Compute actual portfolio return
  ret_port[i] <- sum(w * as.numeric(test[i, ]))
}

# Convert to xts
dates_test <- index(test)
ret_xts <- xts(ret_port, order.by = dates_test)
VaR_xts <- xts(VaR_garch, order.by = dates_test)

# Violations
violations <- ret_xts < VaR_xts
violation_pct <- mean(violations, na.rm = TRUE) * 100

# Plot
par(mfrow = c(1,1), mar = c(4,4,3,2))
plot(ret_xts, type = "l", col = "gray",
     main = paste("GARCH(1,1) Forecast VaR 5% - Portfolio\nViolations:", round(violation_pct, 2), "%"),
     ylab = "Log-return", xlab = "")
lines(VaR_xts, col = "red")
points(index(ret_xts)[violations], ret_xts[violations], col = "blue", pch = 19)
legend("bottomleft", legend = c("VaR 5%", "Violations"),
       col = c("red", "blue"), lty = c(1, NA), pch = c(NA, 19), bty = "n")






