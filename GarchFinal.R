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

# Set date range
start_date <- as.Date("2019-04-30")
end_date <- Sys.Date()

# Define stock symbols
symbols <- c("^GSPC", "^NDX", "^GDAXI", "^N225")
names(symbols) <- c("SP500", "NASDAQ100", "DAX", "NIKKEI")

# Download data from Yahoo Finance
getSymbols(symbols, src = "yahoo", from = start_date, to = end_date)

# Create dataframes with open, close prices and log returns
SP500_df <- data.frame(
  date = index(GSPC),
  open = as.numeric(Op(GSPC)),
  close = as.numeric(Cl(GSPC)),
  log_return = as.numeric(dailyReturn(Cl(GSPC), type = "log"))
)

NASDAQ100_df <- data.frame(
  date = index(NDX),
  open = as.numeric(Op(NDX)),
  close = as.numeric(Cl(NDX)),
  log_return = as.numeric(dailyReturn(Cl(NDX), type = "log"))
)

DAX_df <- data.frame(
  date = index(GDAXI),
  open = as.numeric(Op(GDAXI)),
  close = as.numeric(Cl(GDAXI)),
  log_return = as.numeric(dailyReturn(Cl(GDAXI), type = "log"))
)

NIKKEI_df <- data.frame(
  date = index(N225),
  open = as.numeric(Op(N225)),
  close = as.numeric(Cl(N225)),
  log_return = as.numeric(dailyReturn(Cl(N225), type = "log"))
)

# Remove rows with NA values
SP500_df <- na.omit(SP500_df)
NASDAQ100_df <- na.omit(NASDAQ100_df)
DAX_df <- na.omit(DAX_df)
NIKKEI_df <- na.omit(NIKKEI_df)

# Summary statistics function
summary_stats <- function(df, name) {
  cat("\n---", name, "---\n")
  print(summary(df$log_return))
  cat("Standard Deviation:", sd(df$log_return), "\n")
  cat("Skewness:", skewness(df$log_return), "\n")
  cat("Kurtosis:", kurtosis(df$log_return), "\n")
}

# Stationarity test function
stationarity_tests <- function(x, name) {
  cat("\n---", name, "---\n")
  cat("ADF Test:\n")
  print(adf.test(x))
  cat("KPSS Test:\n")
  print(kpss.test(x))
}

# ARCH LM test
arch_lm_test <- function(x, name) {
  cat("\n---", name, "---\n")
  print(ArchTest(x, lags = 10))
}

# GARCH(1,1) specification
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "norm"
)

# Fit GARCH(1,1) models
fit_SP500 <- ugarchfit(spec, data = SP500_df$log_return)
fit_NASDAQ <- ugarchfit(spec, data = NASDAQ100_df$log_return)
fit_DAX <- ugarchfit(spec, data = DAX_df$log_return)
fit_NIKKEI <- ugarchfit(spec, data = NIKKEI_df$log_return)

# Function for outlier analysis
analyze_exceedances <- function(fit, real_data, name) {
  sigma <- sigma(fit)
  mu <- fitted(fit)
  real_data <- tail(real_data, length(mu))
  upper <- mu + 2 * sigma
  lower <- mu - 2 * sigma
  exceedances <- (real_data > upper) | (real_data < lower)
  cat("\n---", name, "---\n")
  cat("Total observations:", length(real_data), "\n")
  cat("Exceedances:", sum(exceedances, na.rm = TRUE), "\n")
  cat("Observed frequency:", round(mean(exceedances, na.rm = TRUE) * 100, 2), "%\n")
}

# Residual analysis
test_garch_adequacy <- function(fit, name) {
  cat("\n---", name, "---\n")
  resid <- residuals(fit, standardize = TRUE)
  print(Box.test(resid, lag = 20, type = "Ljung-Box"))
  print(Box.test(resid^2, lag = 20, type = "Ljung-Box"))
}

# DCC-GARCH
returns_xts <- na.omit(merge(
  xts(SP500_df$log_return, order.by = SP500_df$date),
  xts(NASDAQ100_df$log_return, order.by = NASDAQ100_df$date),
  xts(DAX_df$log_return, order.by = DAX_df$date),
  xts(NIKKEI_df$log_return, order.by = NIKKEI_df$date)
))
colnames(returns_xts) <- c("SP500", "NASDAQ100", "DAX", "NIKKEI")

uspec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                    mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
                    distribution.model = "norm")
mspec <- multispec(replicate(4, uspec))
dcc_spec <- dccspec(uspec = mspec, dccOrder = c(1,1), distribution = "mvnorm")
dcc_fit <- dccfit(dcc_spec, data = returns_xts)

# Extract rolling correlations
dates <- index(returns_xts)
cor_array <- rcor(dcc_fit)
asset_names <- colnames(returns_xts)
pair_combinations <- combn(asset_names, 2)
pair_names <- apply(pair_combinations, 2, function(x) paste(x[1], "vs", x[2]))

cor_df <- data.frame(Date = rep(dates, times = ncol(pair_combinations)))
for (i in 1:ncol(pair_combinations)) {
  s1 <- which(asset_names == pair_combinations[1, i])
  s2 <- which(asset_names == pair_combinations[2, i])
  cor_df[[pair_names[i]]] <- cor_array[s1, s2, ]
}

# Compute rolling statistics
window_size <- 90
rolling_list <- lapply(pair_names, function(pair) {
  series <- zoo(cor_df[[pair]], order.by = cor_df$Date)
  roll_mean <- rollapply(series, width = window_size, FUN = mean, align = "right", fill = NA)
  roll_sd <- rollapply(series, width = window_size, FUN = sd, align = "right", fill = NA)
  data.frame(Date = index(series), Pair = pair,
             RollingMean = coredata(roll_mean),
             RollingSD = coredata(roll_sd))
})

rolling_df <- bind_rows(rolling_list)

# Plot rolling correlations
plot <- ggplot(rolling_df, aes(x = Date, y = RollingMean)) +
  geom_line(color = "darkblue", size = 0.7) +
  geom_ribbon(aes(ymin = RollingMean - RollingSD, ymax = RollingMean + RollingSD),
              alpha = 0.2, fill = "red") +
  facet_wrap(~ Pair, scales = "free_y") +
  theme_minimal() +
  labs(title = "Rolling Mean and Variability of DCC-GARCH Correlations",
       subtitle = paste("Rolling window:", window_size, "days"),
       y = "Rolling Correlation", x = "Date")
print(plot)
