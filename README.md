# Financial Volatility Modeling and Forecasting with GARCH and DCC-GARCH

This repository provides a full R pipeline for modeling, analyzing, and forecasting financial time series using **GARCH(1,1)** and **DCC-GARCH** models. The analysis is applied to four major equity indices: **S&P 500**, **NASDAQ 100**, **DAX**, and **Nikkei 225**, using daily data over a 4-year period.

## 📈 Project Objectives

- Model volatility clustering in global equity indices
- Estimate GARCH(1,1) models and evaluate their adequacy
- Extend to multivariate DCC-GARCH models to capture dynamic correlations
- Forecast returns and compute Value-at-Risk (VaR)
- Evaluate forecast accuracy and correlation stability

## 🔧 Tools & Libraries

- `quantmod`: Data acquisition from Yahoo Finance
- `rugarch`: Univariate GARCH modeling and forecasting
- `rmgarch`: Multivariate GARCH (DCC) modeling
- `PerformanceAnalytics`: Risk and return measures
- `ggplot2`, `gridExtra`: Visualizations

## 🗃️ Dataset

Daily open and close prices for:

- S&P 500 (^GSPC)
- NASDAQ 100 (^NDX)
- DAX (^GDAXI)
- Nikkei 225 (^N225)

Collected using `getSymbols()` from Yahoo Finance between April 2019 and April 2023.

## 🧠 Methodology

1. **Preprocessing**:
   - Compute log returns from close prices
   - Remove missing values and check stationarity (ADF & KPSS tests)

2. **Univariate Volatility Modeling**:
   - Fit `GARCH(1,1)` models on each index
   - Analyze residuals, ACF plots, and standardized errors
   - Validate adequacy via Ljung-Box and ARCH tests

3. **Multivariate Volatility Modeling**:
   - Estimate `DCC-GARCH` on log return series
   - Extract dynamic correlations and assess rolling stability
   - Visualize time-varying correlations with confidence bands

4. **Forecasting and VaR**:
   - Generate 1-step ahead forecasts (last 10 days)
   - Plot forecasts with 95% confidence intervals
   - Compute Value-at-Risk (VaR) and assess predictive quality

5. **Forecast Evaluation**:
   - MAE, MSE, and RMSE calculated for each forecast series

## 📊 Outputs

- Conditional volatility plots
- Standardized residual diagnostics
- Dynamic correlation matrices (DCC)
- Rolling means and standard deviations of DCC correlations
- Forecast plots with confidence intervals
- Value-at-Risk estimates
- Forecast error metrics

## 📌 Conclusion

- **GARCH(1,1)** proved sufficient for modeling conditional volatility
- **DCC-GARCH** effectively captured evolving correlations across markets
- VaR predictions were within expected statistical bounds
- Model diagnostics and forecast accuracy suggest reliability for risk management

## 📝 Author

This project was developed as part of an academic and applied research initiative on financial econometrics and risk modeling. 

All rights reserved to Giacomo Fantato.


## 📄 License

This project is licensed under the MIT License.

---

Feel free to fork, star, or contribute!
