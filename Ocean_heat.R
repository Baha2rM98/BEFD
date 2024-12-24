library(readxl)
library(lmtest) 
library(forecast)
library(DIMORA)
library(fpp2)
library(dplyr)
library(car)
#####################################################################################################################
# Data presentation
######################################################################################################################
df <- read.csv("data/monthly-ocean-heat-2000m.csv")
View(df)
attach(df)
head(df)
Ocean_heat = Monthly.average.ocean.heat.content.for.the.0.2000.meters.layer
Date = Day
tsdisplay(Ocean_heat, lag.max = 600)
Ocean_heat_ts <- ts(Ocean_heat, frequency = 12, start = c(2005, 1))
seasonplot(Ocean_heat_ts, ylab = "Monthly Averag Ocean Heat Content (2000 Meters Layer)", xlab = "Monthly Time Frame",
           main = "Seasonal Plot: Monthly Averag Ocean Heat Since 2005 ", 
           year.labels = TRUE, year.labels.left = TRUE, 
           col = 1:40, pch = 19)

sum(is.na(Ocean_heat_ts))
decomposition = stl(Ocean_heat_ts , s.window = 'periodic')
plot(decomposition)
################################################################################################################
# Linear Regression
################################################################################################################
tt<- 1:NROW(df)
model1 <- lm(Ocean_heat~ tt)
summary(model1)
anova(model1)
accuracy(model1)

##plot of the model
plot(tt, Ocean_heat, xlab="Monthly Time Frame", ylab="Monthly Averag Ocean Heat Content (2000 Meters Layer)")
abline(model1, col=6, lwd = 2)

dwtest(model1)

resmodel1<- residuals(model1)
plot(resmodel1,xlab="Monthly Time Frame", ylab="Monthly Averag Ocean Heat Content (2000 Meters Layer)" )
Acf(resmodel1, lag.max = 600)
pacf(resmodel1, lag.max = 600)
###################################################################################################################
# TSLM model  
###################################################################################################################
ts.plot(Ocean_heat_ts, type="o")

## we fit a linear model with the tslm function
TSLM_model<- tslm(Ocean_heat_ts~ trend + season)

###obviously it gives the same results of the first model
summary(TSLM_model)
accuracy(TSLM_model)

plot(Ocean_heat_ts, xlab="Monthly Time Frame", ylab="Monthly Averag Ocean Heat Content (2000 Meters Layer)")
lines(fitted(TSLM_model), col=2)

TSLM_res<- residuals(TSLM_model)
plot(TSLM_res,xlab="Monthly Time Frame", ylab="residuals")
Acf(TSLM_res, lag.max = 600)
pacf(TSLM_res, lag.max = 600)
dwtest(TSLM_model)

TSLM_fore <- forecast(TSLM_model, h = 12)
plot(TSLM_fore)
####################################################################################################################
##Exponential Smoothing methods
####################################################################################################################
##1.Simple exponential smoothing
autoplot(Ocean_heat_ts)+ylab("Monthly Averag Ocean Heat Content (2000 Meters Layer)")+xlab("Monthly Time Frame")

fit1<- ses(Ocean_heat_ts, alpha=0.2, initial="simple", h=12)
fit2<- ses(Ocean_heat_ts, alpha=0.6, initial="simple", h=12)
fit3<- ses(Ocean_heat_ts, h=12)


par(mfrow = c(3,1))

plot(Ocean_heat_ts, ylab = "Monthly Averag Ocean Heat Content (2000 Meters Layer)", 
     xlab = "Monthly Time Frame", main = "Original Time Series with Fit 1")
lines(fitted(fit1), col = "blue", type = "o", lwd = 2)

plot(Ocean_heat_ts, ylab = "Monthly Averag Ocean Heat Content (2000 Meters Layer)", 
     xlab = "Monthly Time Frame", main = "Original Time Series with Fit 2")
lines(fitted(fit2), col = "red", type = "o", lwd = 2)


plot(Ocean_heat_ts, ylab = "Monthly Averag Ocean Heat Content (2000 Meters Layer)", 
     xlab = "Monthly Time Frame", main = "Original Time Series with Fit 3")
lines(fitted(fit3), col = "green", type = "o", lwd = 2)
par(mfrow = c(1,1))


fc<- ses(Ocean_heat_ts, h=12)
accuracy(fc)

summary(fc)

autoplot(fc)+
  autolayer(fitted(fc), series="Fitted")+ylab("Quarterly Sea Level")+xlab("Quarterly Time Frame")

##2.Trend methods (Holt method)
fc_Holt<- holt(Ocean_heat_ts, h=12)
fc2_Holt<- holt(Ocean_heat_ts, damped=T, phi=0.9, h=12)

autoplot(Ocean_heat_ts)+
  autolayer(fc_Holt, series="Holt's method", PI=F)+
  autolayer(fc2_Holt, series="Damped Holt's method", PI=F)

###3.Trend and seasonality methods (Holt-Winters method)
autoplot(Ocean_heat_ts)

fit1<- hw(Ocean_heat_ts, seasonal="additive")
fit2<- hw(Ocean_heat_ts, seasonal="multiplicative")

autoplot(Ocean_heat_ts, series = "Original Time Series") +
  autolayer(fit1, series = "Additive Holt-Winters", PI = FALSE) +
  autolayer(fit2, series = "Multiplicative Holt-Winters", PI = FALSE) +
  ggtitle("Holt-Winters Forecasts with Different Seasonal Models") +
  xlab("Quarterly Time Frame") +
  ylab("Quarterly Sea Level") +
  theme_minimal() +
  scale_color_manual(
    values = c("Original Time Series" = "black",
               "Additive Holt-Winters" = "blue", 
               "Multiplicative Holt-Winters" = "red")
  )

# Evaluate performance of fit1 and fit2
accuracy_fit1 <- accuracy(fit1)
accuracy_fit2 <- accuracy(fit2)

par(mfrow = c(2,2))
res_fit1 <- residuals(fit1)
Acf(res_fit1, main = "ACF of Residuals (Additive)", lag.max = 22*12)
pacf(res_fit1, main = "PACF of Residuals (Additive)",lag.max = 22*12)

res_fit2 <- residuals(fit2)
Acf(res_fit2, main = "ACF of Residuals (Multiplicative)",lag.max = 22*12)
pacf(res_fit2, main = "PACF of Residuals (Multiplicative)",lag.max = 22*12)
par(mfrow = c(1,1))
###################################################################################################################
# Arima Models
####################################################################################################################
diff1 <- diff(Ocean_heat_ts)
tsdisplay(diff1, lag.max = 200)

# Seasonal differencing to remove the seasonal pattern
diff_seasonal <- diff(diff1, lag = 12)

# Display ACF and PACF after differencing
tsdisplay(diff_seasonal, lag.max = 200)


# SARIMA model
arima_model1 <- Arima(Ocean_heat_ts, order = c(1, 1, 1), seasonal = list(order = c(1, 0, 0), period = 12))
fit1<- fitted(arima_model1)

plot(Ocean_heat_ts)
lines(fit1, col=2)

f1<- forecast(arima_model1,h = 12)
plot(f1)

r1<- residuals(arima_model1)
tsdisplay(r2, lag.max = 300) 

accuracy_arima_fit1 <- accuracy(f1)


# second arima model (auto arima)
auto.a<- auto.arima(Ocean_heat_ts)
auto.a_fit1<- fitted(auto.a)

plot(Ocean_heat_ts)
lines(auto.a_fit1, col=2)

f2<- forecast(auto.a,h = 12)
plot(f2)

r2<- residuals(auto.a)
tsdisplay(r2, lag.max = 200, main = 'Residuals for Auto-Arima model')

accuracy_arima_fit2 <- accuracy(f2)