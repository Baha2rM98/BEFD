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
df <- read.csv("data/sea-level.csv")
View(df)
attach(df)
head(df)
sea_level = Global.sea.level.according.to.Church.and.White..2011.[4:563]
Date = Day
tsdisplay(sea_level, lag.max = 600)
sea_level_ts <- ts(sea_level, frequency = 4, start = c(1881, 1))
seasonplot(sea_level_ts, ylab = "Quarterly Sea Level", xlab = "Quarterly Time Frame",
           main = "Seasonal Plot: Sea Level from 1881",
           year.labels = TRUE, year.labels.left = TRUE,
           col = 1:40, pch = 19)

sum(is.na(sea_level_ts))
sea_level_ts_new = na.omit(sea_level_ts)
decomposition = stl(sea_level_ts_new, s.window = 'periodic')
plot(decomposition)
################################################################################################################
# Linear Regression
################################################################################################################
tt <- 1:NROW(df[4:563,])
model1 <- lm(sea_level ~ tt)
summary(model1)
anova(model1)
accuracy(model1)

##plot of the model
plot(tt, sea_level, xlab = "Quarterly Time Frame", ylab = "Quarterly Sea Level")
abline(model1, col = 6, lwd = 2)

dwtest(model1)

resmodel1 <- residuals(model1)
plot(resmodel1, xlab = "Quarterly Time Frame", ylab = "Quarterly Sea Level")
Acf(resmodel1, lag.max = 600)
pacf(resmodel1, lag.max = 600)
###################################################################################################################
# TSLM model  
###################################################################################################################
ts.plot(sea_level_ts_new, type = "o")

## we fit a linear model with the tslm function
TSLM_model <- tslm(sea_level_ts_new ~ trend + season)

###obviously it gives the same results of the first model
summary(TSLM_model)
accuracy(TSLM_model)

plot(sea_level_ts_new, xlab = "Quarterly Time Frame", ylab = "Quarterly Sea Level")
lines(fitted(TSLM_model), col = 2)

TSLM_res <- residuals(TSLM_model)
plot(TSLM_res, xlab = "Time(Month of Year)", ylab = "residuals")
Acf(TSLM_res, lag.max = 600)
pacf(TSLM_res, lag.max = 600)
dwtest(TSLM_model)

TSLM_fore <- forecast(TSLM_model, h = 3)
plot(TSLM_fore)
####################################################################################################################
##Exponential Smoothing methods
####################################################################################################################
##1.Simple exponential smoothing
autoplot(sea_level_ts_new) +
  ylab("Quarterly Sea Level") +
  xlab("Quarterly Time Frame")

fit1 <- ses(sea_level_ts_new, alpha = 0.2, initial = "simple", h = 3)
fit2 <- ses(sea_level_ts_new, alpha = 0.6, initial = "simple", h = 3)
fit3 <- ses(sea_level_ts_new, h = 3)


par(mfrow = c(3, 1))

plot(sea_level_ts_new, ylab = "Quarterly Sea Level",
     xlab = "Quarterly Time Frame", main = "Original Time Series with Fit 1")
lines(fitted(fit1), col = "blue", type = "o", lwd = 2)

plot(sea_level_ts_new, ylab = "Quarterly Sea Level",
     xlab = "Quarterly Time Frame", main = "Original Time Series with Fit 2")
lines(fitted(fit2), col = "red", type = "o", lwd = 2)


plot(sea_level_ts_new, ylab = "Quarterly Sea Level",
     xlab = "Quarterly Time Frame", main = "Original Time Series with Fit 3")
lines(fitted(fit3), col = "green", type = "o", lwd = 2)
par(mfrow = c(1, 1))


fc <- ses(sea_level_ts_new, h = 3)
accuracy(fc)

summary(fc)

autoplot(fc) +
  autolayer(fitted(fc), series = "Fitted") +
  ylab("Quarterly Sea Level") +
  xlab("Quarterly Time Frame")

##2.Trend methods (Holt method)
fc_Holt <- holt(sea_level_ts_new, h = 10)
fc2_Holt <- holt(sea_level_ts_new, damped = T, phi = 0.9, h = 10)

autoplot(sea_level_ts_new) +
  autolayer(fc_Holt, series = "Holt's method", PI = F) +
  autolayer(fc2_Holt, series = "Damped Holt's method", PI = F)

###3.Trend and seasonality methods (Holt-Winters method)
autoplot(sea_level_ts_new)

fit1 <- hw(sea_level_ts_new, seasonal = "additive")

autoplot(sea_level_ts_new, series = "Original Time Series") +
  autolayer(fit1, series = "Additive Holt-Winters", PI = FALSE) +
  ggtitle("Holt-Winters Forecasts with Different Seasonal Models") +
  xlab("Quarterly Time Frame") +
  ylab("Quarterly Sea Level") +
  theme_minimal() +
  scale_color_manual(
    values = c("Original Time Series" = "black",
               "Additive Holt-Winters" = "blue")
  )

# Evaluate performance of fit1 and fit2
accuracy_fit1 <- accuracy(fit1)

par(mfrow = c(2, 1))
res_fit1 <- residuals(fit1)
Acf(res_fit1, main = "ACF of Residuals (Additive)", lag.max = 600)
pacf(res_fit1, main = "PACF of Residuals (Additive)", lag.max = 600)
par(mfrow = c(1, 1))

###################################################################################################################
# Arima Models
####################################################################################################################
diff1 <- diff(sea_level_ts_new)

# Seasonal differencing to remove the seasonal pattern
diff_seasonal <- diff(diff1, lag = 12)

# Display ACF and PACF after differencing
tsdisplay(diff_seasonal, lag.max = 600)


# SARIMA model
arima_model1 <- Arima(sea_level_ts_new, order = c(1, 1, 1), seasonal = list(order = c(1, 0, 0), period = 4))
fit1 <- fitted(arima_model1)

plot(sea_level_ts_new)
lines(fit1, col = 2)

f1 <- forecast(arima_model1, h = 4)
plot(f1)

r1 <- residuals(arima_model1)
tsdisplay(r2, lag.max = 500)

accuracy_arima_fit1 <- accuracy(f1)


# second arima model (auto arima)
auto.a <- auto.arima(sea_level_ts_new)
auto.a_fit1 <- fitted(auto.a)

plot(sea_level_ts_new)
lines(auto.a_fit1, col = 2)

f2 <- forecast(auto.a, h = 4)
plot(f2)

r2 <- residuals(auto.a)
tsdisplay(r2, lag.max = 500, main = 'Residuals for Auto-Arima model')

accuracy_arima_fit2 <- accuracy(f2)