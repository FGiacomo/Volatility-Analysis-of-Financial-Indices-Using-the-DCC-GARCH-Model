# Install and load required packages
# install.packages("quantmod")  # only if necessary
library(quantmod)

# Set the date range
start_date <- as.Date("2019-04-30")
end_date <- Sys.Date()

# Define ticker symbols and corresponding names
symbols <- c("^GSPC", "^NDX", "^GDAXI", "^N225")
names(symbols) <- c("SP500", "NASDAQ100", "DAX", "NIKKEI")

# Download data from Yahoo Finance
getSymbols(symbols, src = "yahoo", from = start_date, to = end_date)

# Create dataframes with open, close prices, and log returns

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

# Preview data
head(SP500_df)
head(NASDAQ100_df)
head(DAX_df)
head(NIKKEI_df)

# Check for NA values in each dataframe
sapply(list(SP500 = SP500_df,
            NASDAQ100 = NASDAQ100_df,
            DAX = DAX_df,
            NIKKEI = NIKKEI_df), function(df) colSums(is.na(df)))

# Remove rows with NA values (just in case)
SP500_df <- na.omit(SP500_df)
NASDAQ100_df <- na.omit(NASDAQ100_df)
DAX_df <- na.omit(DAX_df)
NIKKEI_df <- na.omit(NIKKEI_df)

# Plot log returns (using base R)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(SP500_df$date, SP500_df$log_return, type = "l", main = "S&P 500 Log Returns", xlab = "Date", ylab = "Log Return")
plot(NASDAQ100_df$date, NASDAQ100_df$log_return, type = "l", main = "NASDAQ 100 Log Returns", xlab = "Date", ylab = "Log Return")
plot(DAX_df$date, DAX_df$log_return, type = "l", main = "DAX Log Returns", xlab = "Date", ylab = "Log Return")
plot(NIKKEI_df$date, NIKKEI_df$log_return, type = "l", main = "Nikkei 225 Log Returns", xlab = "Date", ylab = "Log Return")

# Descriptive statistics
summary(SP500_df$log_return)
summary(NASDAQ100_df$log_return)
summary(DAX_df$log_return)
summary(NIKKEI_df$log_return)


## summary statics
library(moments)

summary_stats <- function(df, name) {
  cat("\n---", name, "---\n")
  print(summary(df$log_return))
  cat("Standard Deviation:", sd(df$log_return), "\n")
  cat("Skewness:", skewness(df$log_return), "\n")
  cat("Kurtosis:", kurtosis(df$log_return), "\n")
}

summary_stats(SP500_df, "S&P 500")
summary_stats(NASDAQ100_df, "NASDAQ 100")
summary_stats(DAX_df, "DAX")
summary_stats(NIKKEI_df, "Nikkei 225")


par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

plot(SP500_df$date, SP500_df$close, type = "l", main = "S&P 500 Price", xlab = "Date", ylab = "Price")
plot(NASDAQ100_df$date, NASDAQ100_df$close, type = "l", main = "NASDAQ 100 Price", xlab = "Date", ylab = "Price")
plot(DAX_df$date, DAX_df$close, type = "l", main = "DAX Price", xlab = "Date", ylab = "Price")
plot(NIKKEI_df$date, NIKKEI_df$close, type = "l", main = "Nikkei 225 Price", xlab = "Date", ylab = "Price")


## ---- Stationarity tests (ADF and KPSS)
install.packages("tseries")  # if not already installed
library(tseries)

# ADF and KPSS tests for each series
stationarity_tests <- function(x, name) {
  cat("\n---", name, "---\n")
  cat("ADF Test:\n")
  print(adf.test(x))
  cat("KPSS Test:\n")
  print(kpss.test(x))
}

stationarity_tests(SP500_df$log_return, "S&P 500")
stationarity_tests(NASDAQ100_df$log_return, "NASDAQ 100")
stationarity_tests(DAX_df$log_return, "DAX")
stationarity_tests(NIKKEI_df$log_return, "Nikkei 225")



##  Check for volatility clustering and ARCH effects
# Plot ACF of squared returns
par(mfrow = c(2, 2))
acf(SP500_df$log_return^2, main = "S&P 500 Squared Returns ACF")
acf(NASDAQ100_df$log_return^2, main = "NASDAQ 100 Squared Returns ACF")
acf(DAX_df$log_return^2, main = "DAX Squared Returns ACF")
acf(NIKKEI_df$log_return^2, main = "Nikkei 225 Squared Returns ACF")

# Perform ARCH LM test
#install.packages("FinTS")
library(FinTS)

arch_lm_test <- function(x, name) {
  cat("\n---", name, "---\n")
  print(ArchTest(x, lags = 10))
}

arch_lm_test(SP500_df$log_return, "S&P 500")
arch_lm_test(NASDAQ100_df$log_return, "NASDAQ 100")
arch_lm_test(DAX_df$log_return, "DAX")
arch_lm_test(NIKKEI_df$log_return, "Nikkei 225")


## FIT GARCH(1,1) model
#install.packages("rugarch")  # Only if not already installed
library(rugarch)
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "norm"
)
# Fit GARCH model for each index
fit_SP500 <- ugarchfit(spec, data = SP500_df$log_return)
fit_NASDAQ <- ugarchfit(spec, data = NASDAQ100_df$log_return)
fit_DAX <- ugarchfit(spec, data = DAX_df$log_return)
fit_NIKKEI <- ugarchfit(spec, data = NIKKEI_df$log_return)
# Example for SP500
show(fit_SP500)
show(fit_NASDAQ)
show(fit_DAX)
show(fit_NIKKEI)

# Conditional standard deviation -----------------------
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))  # aumenta lo spazio sopra (mar[3])
plot(fit_SP500, which = 1)
mtext("S&P 500", side = 3, line = 0.1, font = 1, cex = 0.9)
plot(fit_NASDAQ, which = 1)
mtext("NASDAQ 100", side = 3, line = 0.1, font = 1, cex = 0.9)
plot(fit_DAX, which = 1)
mtext("DAX", side = 3, line = 0.1, font = 1, cex = 0.9)
plot(fit_NIKKEI, which = 1)
mtext("Nikkei 225", side = 3, line = 0.1, font = 1, cex = 0.9)
## 1.1 couting the exceed:
library(rugarch)
library(PerformanceAnalytics)

# Funzione per analizzare sforamenti ±2σ --> not EXECUTE (i'not sure about that)
analisi_sforamenti <- function(fit, real_data, nome_indice) {
  sigma <- sigma(fit)
  mu <- fitted(fit)
  
  # Allineamento temporale
  real_data <- tail(real_data, length(mu))
  
  upper <- mu + 2 * sigma
  lower <- mu - 2 * sigma
  
  sforamenti <- (real_data > upper) | (real_data < lower)
  n_total <- length(real_data)
  n_sforamenti <- sum(sforamenti, na.rm = TRUE)
  freq_sforamenti <- n_sforamenti / n_total
  
  cat("\n---", nome_indice, "---\n")
  cat("Totale osservazioni:", n_total, "\n")
  cat("Numero di sforamenti:", n_sforamenti, "\n")
  cat("Frequenza osservata:", round(freq_sforamenti * 100, 2), "%\n")
  cat("Frequenza attesa (±2σ): 4.6%\n")
}
# Applica la funzione ai 4 modelli stimati
analisi_sforamenti(fit_SP500, SP500_df$log_return, "S&P 500")
analisi_sforamenti(fit_NASDAQ, NASDAQ100_df$log_return, "NASDAQ 100")
analisi_sforamenti(fit_DAX, DAX_df$log_return, "DAX")
analisi_sforamenti(fit_NIKKEI, NIKKEI_df$log_return, "Nikkei 225")
## ------------------------------------------ 



# Standardized residuals
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(fit_SP500, which = 2)
mtext("S&P 500", side = 3, line = 0.1, font = 0.5, cex = 0.9)
plot(fit_NASDAQ, which = 2)
mtext("NASDAQ 100", side = 3, line = 0.1, font = 0.5, cex = 0.9)
plot(fit_DAX, which = 2)
mtext("DAX", side = 3, line = 0.1, font = 0.5, cex = 0.9)
plot(fit_NIKKEI, which = 2)
mtext("Nikkei 225", side = 3, line = 0.1, font = 0.5, cex = 0.9)


# ACF of squared residuals
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(fit_SP500, which = 3)
mtext("S&P 500", side = 3, line = 0.1, font = 1, cex = 0.9)
plot(fit_NASDAQ, which = 3)
mtext("NASDAQ 100", side = 3, line = 0.1, font = 1, cex = 0.9)
plot(fit_DAX, which = 3)
mtext("DAX", side = 3, line = 0.1, font = 1, cex = 0.9)
plot(fit_NIKKEI, which = 3)
mtext("Nikkei 225", side = 3, line = 0.1, font = 1, cex = 0.9)


# Q-Q plot
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(fit_SP500, which = 10)
mtext("S&P 500", side = 3, line = 0.1, font = 1, cex = 0.9)
plot(fit_NASDAQ, which = 10)
mtext("NASDAQ 100", side = 3, line = 0.1, font = 1, cex = 0.9)
plot(fit_DAX, which = 10)
mtext("DAX", side = 3, line = 0.1, font = 1, cex = 0.9)
plot(fit_NIKKEI, which = 10)
mtext("Nikkei 225", side = 3, line = 0.1, font = 1, cex = 0.9)


## fianl analysis over residuals --> se the specifications
test_garch_adequacy <- function(fit, name) {
  cat("\n---", name, "---\n")
  
  resid <- residuals(fit, standardize = TRUE)
  
  # Ljung-Box test on residuals
  test_resid <- Box.test(resid, lag = 20, type = "Ljung-Box")
  
  # Ljung-Box test on squared residuals
  test_resid_sq <- Box.test(resid^2, lag = 20, type = "Ljung-Box")
  
  cat("Ljung-Box on residuals (autocorrelation):\n")
  print(test_resid)
  
  cat("Ljung-Box on squared residuals (ARCH effect):\n")
  print(test_resid_sq)
  
  # Decision logic
  if (test_resid$p.value > 0.05 && test_resid_sq$p.value > 0.05) {
    cat("GARCH(1,1) appears adequate for", name, "\n")
  } else {
    cat("Consider extending the model for", name, ":\n")
    if (test_resid$p.value <= 0.05) cat("   → Residuals show autocorrelation.\n")
    if (test_resid_sq$p.value <= 0.05) cat("   → ARCH effects remain in squared residuals.\n")
    cat("   → Suggest trying GARCH(p,q), GJR-GARCH, or t-distribution.\n")
  }
}
test_garch_adequacy(fit_SP500, "S&P 500")
test_garch_adequacy(fit_NASDAQ, "NASDAQ 100")
test_garch_adequacy(fit_DAX, "DAX")
test_garch_adequacy(fit_NIKKEI, "Nikkei 225")

## ------------- DCC GARCH ----------------

library(rmgarch)
library(xts)
library(quantmod)

# 1. Combine log returns into one xts object
returns_xts <- na.omit(merge(
  xts(SP500_df$log_return, order.by = SP500_df$date),
  xts(NASDAQ100_df$log_return, order.by = NASDAQ100_df$date),
  xts(DAX_df$log_return, order.by = DAX_df$date),
  xts(NIKKEI_df$log_return, order.by = NIKKEI_df$date)
))
colnames(returns_xts) <- c("SP500", "NASDAQ100", "DAX", "NIKKEI")

# 2. Define univariate GARCH(1,1) specification for all series
uspec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                    mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
                    distribution.model = "norm")
mspec <- multispec(replicate(4, uspec))

# 3. Define and estimate DCC-GARCH model
dcc_spec <- dccspec(uspec = mspec, dccOrder = c(1,1), distribution = "mvnorm")
dcc_fit <- dccfit(dcc_spec, data = returns_xts)

# 4. Extract conditional correlations
correlations <- rcor(dcc_fit)
time_index <- index(returns_xts)

# 5. Plot pairwise dynamic conditional correlations
pairs <- combn(1:4, 2)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
for (i in 1:ncol(pairs)) {
  series1 <- pairs[1, i]
  series2 <- pairs[2, i]
  plot(time_index, correlations[series1, series2, ], type = "l",
       main = paste(colnames(returns_xts)[series1], "vs", colnames(returns_xts)[series2]),
       xlab = "Date", ylab = "Correlation")
}
# 6. Display the DCC-GARCH fit summary
summary(dcc_fit)


## ------------ Rolling Mean and Variability of DCC-GARCH Correlations ----------------
library(rmgarch)
library(zoo)
library(xts)
library(ggplot2)
library(reshape2)
library(dplyr)

# 1. Extract correlations from DCC fit
cor_array <- rcor(dcc_fit)  # 3D array: asset x asset x time
dates <- index(returns_xts)

# 2. Prepare pairwise labels
asset_names <- colnames(returns_xts)
pair_combinations <- combn(asset_names, 2)
pair_names <- apply(pair_combinations, 2, function(x) paste(x[1], "vs", x[2]))

# 3. Create correlation dataframe
cor_df <- data.frame(Date = rep(dates, times = ncol(pair_combinations)))
for (i in 1:ncol(pair_combinations)) {
  s1 <- which(asset_names == pair_combinations[1, i])
  s2 <- which(asset_names == pair_combinations[2, i])
  cor_df[[pair_names[i]]] <- cor_array[s1, s2, ]
}

# 4. Set rolling window size
window_size <- 90

# 5. Compute rolling statistics
rolling_list <- lapply(pair_names, function(pair) {
  series <- zoo(cor_df[[pair]], order.by = cor_df$Date)
  roll_mean <- rollapply(series, width = window_size, FUN = mean, align = "right", fill = NA)
  roll_sd <- rollapply(series, width = window_size, FUN = sd, align = "right", fill = NA)
  data.frame(Date = index(series),
             Pair = pair,
             RollingMean = coredata(roll_mean),
             RollingSD = coredata(roll_sd))
})

rolling_df <- bind_rows(rolling_list)

# 6. Plot
ggplot(rolling_df, aes(x = Date, y = RollingMean)) +
  geom_line(color = "darkblue", size = 0.7) +
  geom_ribbon(aes(ymin = RollingMean - RollingSD, ymax = RollingMean + RollingSD),
              alpha = 0.2, fill = "red") +
  facet_wrap(~ Pair, scales = "free_y") +
  theme_minimal() +
  labs(title = "Rolling Mean and Variability of DCC-GARCH Correlations",
       subtitle = paste("Rolling window:", window_size, "days"),
       y = "Rolling Correlation", x = "Date")

# VaR - 1 step ahed
# Funzione per calcolare il VaR a 1-step ahead da un oggetto ugarchfit
calcola_VaR <- function(fit, name, alpha = 0.05) {
  forecast <- ugarchforecast(fit, n.ahead = 1)
  mu <- fitted(forecast)        # previsione del rendimento atteso
  sigma <- sigma(forecast)      # deviazione standard condizionale
  VaR <- mu + qnorm(alpha) * sigma
  
  cat("\n---", name, "---\n")
  cat("Expected return (1-step ahead):", round(mu, 6), "\n")
  cat("Conditional SD:", round(sigma, 6), "\n")
  cat(paste0((1 - alpha) * 100, "% 1-step ahead VaR: "), round(VaR, 6), "\n")
}

# Calcolo VaR per ciascun indice (es. al 95%)
calcola_VaR(fit_SP500, "S&P 500", alpha = 0.05)
calcola_VaR(fit_NASDAQ, "NASDAQ 100", alpha = 0.05)
calcola_VaR(fit_DAX, "DAX", alpha = 0.05)
calcola_VaR(fit_NIKKEI, "Nikkei 225", alpha = 0.05)
