library(readxl)
library(lmtest)
library(forecast)
library(DIMORA)
library(fpp2)
library(dplyr)
library(car)
############################################################################################################
# Data presentation
############################################################################################################
df <- read.csv("data/global-co2-concentration.csv")
View(df)
attach(df)
head(df)
CO2 = Monthly.concentration.of.atmospheric.carbon.dioxide
Date = Day
tsdisplay(CO2, lag.max = 120)
CO2_ts <- ts(CO2, frequency = 12, start = c(1979, 1))
seasonplot(CO2_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide", xlab = "Month of Year",
           main = "Seasonal Plot: concentration of Atmospheric Carbon Dioxide",
           year.labels = TRUE, year.labels.left = TRUE,
           col = 1:40, pch = 19)

decomposition = stl(CO2_ts, s.window = 'periodic')
plot(decomposition)


################################################################################################################
# Linear Regression
################################################################################################################
tt <- 1:NROW(df)
model1 <- lm(CO2 ~ tt)
summary(model1)
anova(model1)
accuracy(model1)

##plot of the model
plot(tt, CO2, xlab = "Time(Month of Year)", ylab = "Monthly concentration of Atmospheric Carbon Dioxide")
abline(model1, col = 6, lwd = 2)

dwtest(model1)

resmodel1 <- residuals(model1)
plot(resmodel1, xlab = "Time(Month of Year)", ylab = "residuals")
Acf(resmodel1, lag.max = 22 * 12)
pacf(resmodel1, lag.max = 22 * 12)
###################################################################################################################
# TSLM model  
###################################################################################################################
ts.plot(CO2_ts, type = "o")

## we fit a linear model with the tslm function
TSLM_model <- tslm(CO2_ts ~ trend + season, data = df)

###obviously it gives the same results of the first model
summary(TSLM_model)
accuracy(TSLM_model)

plot(CO2_ts, xlab = "Time(Month of Year)", ylab = "Monthly concentration of Atmospheric Carbon Dioxide")
lines(fitted(TSLM_model), col = 2)

TSLM_res <- residuals(TSLM_model)
plot(TSLM_res, , xlab = "Time(Month of Year)", ylab = "residuals")
Acf(TSLM_res, lag.max = 22 * 12)
pacf(TSLM_res, lag.max = 22 * 12)
dwtest(TSLM_model)

TSLM_fore <- forecast(TSLM_model, h = 12)
plot(TSLM_fore)

################################################################################
##Exponential Smoothing methods
################################################################################
##1.Simple exponential smoothing
autoplot(CO2_ts) +
  ylab("Monthly concentration of Atmospheric Carbon Dioxide") +
  xlab("Time(Month of Year)")

fit1 <- ses(CO2_ts, alpha = 0.2, initial = "simple", h = 5)
fit2 <- ses(CO2_ts, alpha = 0.6, initial = "simple", h = 5)
fit3 <- ses(CO2_ts, h = 5)


par(mfrow = c(3, 1))

plot(CO2_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide",
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 1")
lines(fitted(fit1), col = "blue", type = "o", lwd = 2)

plot(CO2_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide",
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 2")
lines(fitted(fit2), col = "red", type = "o", lwd = 2)


plot(CO2_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide",
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 3")
lines(fitted(fit3), col = "green", type = "o", lwd = 2)
par(mfrow = c(1, 1))


fc <- ses(CO2_ts, h = 12)
accuracy(fc)

summary(fc)

autoplot(fc) +
  autolayer(fitted(fc), series = "Fitted") +
  ylab("Monthly concentration of Atmospheric Carbon Dioxide") +
  xlab("Time (Month of Year)")


##2.Trend methods (Holt method)
fc_Holt <- holt(CO2_ts, h = 12)
fc2_Holt <- holt(CO2_ts, damped = T, phi = 0.9, h = 12)

autoplot(CO2_ts) +
  autolayer(fc_Holt, series = "Holt's method", PI = F) +
  autolayer(fc2_Holt, series = "Damped Holt's method", PI = F)

###3.Trend and seasonality methods (Holt-Winters method)
autoplot(CO2_ts)

fit1 <- hw(CO2_ts, seasonal = "additive")
fit2 <- hw(CO2_ts, seasonal = "multiplicative")

autoplot(CO2_ts, series = "Original Time Series") +
  autolayer(fit1, series = "Additive Holt-Winters", PI = FALSE) +
  autolayer(fit2, series = "Multiplicative Holt-Winters", PI = FALSE) +
  ggtitle("Holt-Winters Forecasts with Different Seasonal Models") +
  xlab("Time (Month of Year)") +
  ylab("Monthly Concentration of Atmospheric Carbon Dioxide") +
  theme_minimal() +
  scale_color_manual(
    values = c("Original Time Series" = "black",
               "Additive Holt-Winters" = "blue",
               "Multiplicative Holt-Winters" = "red")
  )

# Evaluate performance of fit1 and fit2
accuracy_fit1 <- accuracy(fit1)
accuracy_fit2 <- accuracy(fit2)

par(mfrow = c(2, 2))
res_fit1 <- residuals(fit1)
Acf(res_fit1, main = "ACF of Residuals (Additive)", lag.max = 22 * 12)
pacf(res_fit1, main = "PACF of Residuals (Additive)", lag.max = 22 * 12)

res_fit2 <- residuals(fit2)
Acf(res_fit2, main = "ACF of Residuals (Multiplicative)", lag.max = 22 * 12)
pacf(res_fit2, main = "PACF of Residuals (Multiplicative)", lag.max = 22 * 12)
par(mfrow = c(1, 1))

########################################################################################################
# Arima Models
########################################################################################################
diff1 <- diff(CO2_ts)

# Seasonal differencing to remove the seasonal pattern
diff_seasonal <- diff(diff1, lag = 12)

# Display ACF and PACF after differencing
tsdisplay(diff_seasonal, lag.max = 22 * 12)


# first arima model
arima_model1 <- Arima(CO2_ts, order = c(1, 1, 1), seasonal = c(0, 1, 0))
fit1 <- fitted(arima_model1)

plot(CO2_ts)
lines(fit1, col = 2)

f1 <- forecast(arima_model1, h = 12)
plot(f1)

r1 <- residuals(arima_model1)
tsdisplay(r2, lag.max = 22 * 12)

accuracy_arima_fit1 <- accuracy(f1)


# second arima model (auto arima)
auto.a <- auto.arima(CO2_ts)
auto.a_fit1 <- fitted(auto.a)

plot(CO2_ts)
lines(auto.a_fit1, col = 2)

f2 <- forecast(auto.a, h = 12)
plot(f2)

r2 <- residuals(auto.a)
tsdisplay(r2, lag.max = 22 * 12, main = 'Residuals for Auto-Arima model')

accuracy_arima_fit2 <- accuracy(f2)


##################################################################################################################
# Methane
##################################################################################################################
df2 <- read.csv("D:\\University\\2024 - 2025 B\\Business,Economic, and Financial data\\global-methane-concentrations.csv")
View(df2)
attach(df2)
head(df2)
Methane = Monthly.concentration.of.atmospheric.methane
Date = Day
tsdisplay(Methane, lag.max = 120)
Methane_ts <- ts(Methane, frequency = 12, start = c(1983, 1))
seasonplot(Methane_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide", xlab = "Month of Year",
           main = "Seasonal Plot: concentration of Atmospheric Carbon Dioxide",
           year.labels = TRUE, year.labels.left = TRUE,
           col = 1:40, pch = 19)

decomposition = stl(Methane_ts, s.window = 'periodic')
plot(decomposition)


################################################################################################################
# Linear Regression
################################################################################################################
tt <- 1:NROW(df2)
model1 <- lm(Methane ~ tt)
summary(model1)
anova(model1)
accuracy(model1)

##plot of the model
plot(tt, Methane, xlab = "Time(Month of Year)", ylab = "Monthly concentration of Atmospheric Carbon Dioxide")
abline(model1, col = 6, lwd = 2)

dwtest(model1)

resmodel1 <- residuals(model1)
plot(resmodel1, xlab = "Time(Month of Year)", ylab = "residuals")
Acf(resmodel1, lag.max = 50 * 12)
pacf(resmodel1, lag.max = 50 * 12)
################################################################################
# TSLM model  
################################################################################
ts.plot(Methane_ts, type = "o")

## we fit a linear model with the tslm function
TSLM_model2 <- tslm(Methane_ts ~ trend + season, data = df2)

###obviously it gives the same results of the first model
summary(TSLM_model2)
accuracy(TSLM_model2)

plot(Methane_ts, xlab = "Time(Month of Year)", ylab = "Monthly concentration of Atmospheric Methane")
lines(fitted(TSLM_model2), col = 2)

TSLM_res <- residuals(TSLM_model2)
plot(TSLM_res, , xlab = "Time(Month of Year)", ylab = "residuals")
Acf(TSLM_res, lag.max = 50 * 12)
pacf(TSLM_res, lag.max = 50 * 12)
dwtest(TSLM_model)

TSLM_fore <- forecast(TSLM_model2, h = 12)
plot(TSLM_fore)

################################################################################
##Exponential Smoothing methods
################################################################################
##1.Simple exponential smoothing
autoplot(Methane_ts) +
  ylab("Monthly concentration of Atmospheric Methane") +
  xlab("Time(Month of Year)")

fit1 <- ses(Methane_ts, alpha = 0.2, initial = "simple", h = 5)
fit2 <- ses(Methane_ts, alpha = 0.6, initial = "simple", h = 5)
fit3 <- ses(Methane_ts, h = 5)


par(mfrow = c(3, 1))

plot(Methane_ts, ylab = "Monthly concentration of Atmospheric Methane",
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 1")
lines(fitted(fit1), col = "blue", type = "o", lwd = 2)

plot(Methane_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide",
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 2")
lines(fitted(fit2), col = "red", type = "o", lwd = 2)


plot(Methane_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide",
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 3")
lines(fitted(fit3), col = "green", type = "o", lwd = 2)
par(mfrow = c(1, 1))


fc <- ses(Methane_ts, h = 60)
accuracy(fc)

summary(fc)

autoplot(fc) +
  autolayer(fitted(fc), series = "Fitted") +
  ylab("Monthly concentration of Atmospheric Methane") +
  xlab("Time (Month of Year)")


##2.Trend methods (Holt method)
fc_Holt <- holt(Methane_ts, h = 60)
fc2_Holt <- holt(Methane_ts, damped = T, phi = 0.9, h = 60)

autoplot(Methane_ts) +
  autolayer(fc_Holt, series = "Holt's method", PI = F) +
  autolayer(fc2_Holt, series = "Damped Holt's method", PI = F)

###3.Trend and seasonality methods (Holt-Winters method)
autoplot(Methane_ts)

fit1 <- hw(Methane_ts, seasonal = "additive")
fit2 <- hw(Methane_ts, seasonal = "multiplicative")

autoplot(Methane_ts, series = "Original Time Series") +
  autolayer(fit1, series = "Additive Holt-Winters", PI = FALSE) +
  autolayer(fit2, series = "Multiplicative Holt-Winters", PI = FALSE) +
  ggtitle("Holt-Winters Forecasts with Different Seasonal Models") +
  xlab("Time (Month of Year)") +
  ylab("Monthly Concentration of Atmospheric Methane") +
  theme_minimal() +
  scale_color_manual(
    values = c("Original Time Series" = "black",
               "Additive Holt-Winters" = "blue",
               "Multiplicative Holt-Winters" = "red")
  )

# Evaluate performance of fit1 and fit2
accuracy_fit1 <- accuracy(fit1)
accuracy_fit2 <- accuracy(fit2)

par(mfrow = c(2, 2))
res_fit1 <- residuals(fit1)
Acf(res_fit1, main = "ACF of Residuals (Additive)", lag.max = 50 * 12)
pacf(res_fit1, main = "PACF of Residuals (Additive)", lag.max = 50 * 12)

res_fit2 <- residuals(fit2)
Acf(res_fit2, main = "ACF of Residuals (Multiplicative)", lag.max = 50 * 12)
pacf(res_fit2, main = "PACF of Residuals (Multiplicative)", lag.max = 50 * 12)
par(mfrow = c(1, 1))

########################################################################################################
# Arima Models
########################################################################################################
diff1 <- diff(Methane_ts)

# Seasonal differencing to remove the seasonal pattern
diff_seasonal <- diff(diff1, lag = 12)

# Display ACF and PACF after differencing
tsdisplay(diff_seasonal, lag.max = 50 * 12)


# first arima model
arima_model1 <- Arima(Methane_ts, order = c(1, 1, 1), seasonal = c(0, 1, 0))
fit1 <- fitted(arima_model1)

plot(Methane_ts)
lines(fit1, col = 2)

f1 <- forecast(arima_model1, h = 60)
plot(f1)

r1 <- residuals(arima_model1)
tsdisplay(r2, lag.max = 50 * 12)

accuracy_arima_fit1 <- accuracy(f1)


# second arima model (auto arima)
auto.a <- auto.arima(Methane_ts)
auto.a_fit1 <- fitted(auto.a)

plot(Methane_ts)
lines(auto.a_fit1, col = 2)

f2 <- forecast(auto.a, h = 60)
plot(f2)

r2 <- residuals(auto.a)
tsdisplay(r2, lag.max = 22 * 12, main = 'Residuals for Auto-Arima model')

accuracy_arima_fit2 <- accuracy(f2)
