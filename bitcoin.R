library(TSA)
library(fUnitRoots)
library(AID)
library(lmtest)
library(forecast)
library(rugarch)
library(fGarch)
library(readxl)
library(FitAR)

btc.forecasts <- read_xlsx("~/RMIT/2018 Semester One/Time Series Analysis/Final Project/Bitcoin/Bitcoin_Prices_Forecasts.xlsx")
btc.forecasts <- btc.forecasts$`Closing price`


sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}

MASE = function(observed , fitted ){
  # observed: Observed series on the forecast period
  # fitted: Forecast values by your model
  Y.t = observed
  n = length(fitted)
  e.t = Y.t - fitted
  sum = 0 
  for (i in 2:n){
    sum = sum + abs(Y.t[i] - Y.t[i-1] )
  }
  q.t = e.t / (sum/(n-1))
  MASE = data.frame( MASE = mean(abs(q.t)))
  return(list(MASE = MASE))
}

residual.analysis <- function(model, std = TRUE,start = 2, class = c("ARIMA","GARCH","ARMA-GARCH")[1]){
  # If you have an output from arima() function use class = "ARIMA"
  # If you have an output from garch() function use class = "GARCH"
  # If you have an output from ugarchfit() function use class = "ARMA-GARCH"
  library(TSA)
  library(FitAR)
  if (class == "ARIMA"){
    if (std == TRUE){
      res.model = rstandard(model)
    }else{
      res.model = residuals(model)
    }
  }else if (class == "GARCH"){
    res.model = model$residuals[start:model$n.used]
  }else if (class == "ARMA-GARCH"){
    res.model = model@fit$residuals
  }else {
    stop("The argument 'class' must be either 'ARIMA' or 'GARCH' ")
  }
  par(mfrow=c(3,2))
  plot(res.model,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
  abline(h=0)
  hist(res.model,main="Histogram of standardised residuals")
  acf(res.model,main="ACF of standardised residuals")
  pacf(res.model,main="PACF of standardised residuals")
  qqnorm(res.model,main="QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
  par(mfrow=c(1,1))
}

# Import the data set.
bitcoin <- read.csv("~/RMIT/2018 Semester One/Time Series Analysis/Final Project/Bitcoin/Bitcoin_Historical_Price.csv")

# Create a daily date object.
date <- seq(as.Date("2013-04-27"), as.Date("2018-03-03"), by="day")

# Convert to a ts object.
bitcoin.ts <- ts(bitcoin$Close, start = c(2013, as.numeric(format(date[1], "%j"))), frequency = 365)

class(bitcoin.ts)

plot(bitcoin.ts, ylab='Closing Price ($)', main="Time Series Plot of Bitcoin Closing Price Data")


plot(bitcoin.ts, ylab='Closing Price ($)', main="Time Series Plot of Bitcoin Closing Price Data")
plot(bitcoin.ts, xlim=c(2015.5, 2016.0), ylim=c(200, 500), type='o', 
     main="Close-up Inspection of Bitcoin Time Seires Data between 2015 and 2016")
par(mfrow=c(1,1))
# The series shows multiple trends ands lots of periods where the volitility is high and low.
# This implies the existence of an ARCH effect.
# There's no sign of seasonality and and observations are bouncing around particular mean levels.

# Consider converting the series to a return series and apply GARCH models.

# Convert the series to a return series.
r.bitcoin <- diff(log(bitcoin.ts))

mean(r.bitcoin) # mean is close to zero.

plot(r.bitcoin, main="Daily Return Series for Bitcoin Price Data", ylab='Closing Price ($)')
# Obvious volatility, no sign of trend or seasonality.

McLeod.Li.test(y=r.bitcoin, main = "McLeod-Li Test Statistics for Daily Return Series.")
# Test is highly significant. Normality assumption is highly violated.


qqnorm(r.bitcoin, main = "QQ Plot of Daily Return Bitcoin Price Series.")
qqline(r.bitcoin, col = 2)


# ACF/PACF Plots
par(mfrow=c(1,2))
acf(r.bitcoin, main = "The Sample ACF Plot for the Daily Return Bitcoin Price Series.")
pacf(r.bitcoin, main = "The Sample PACF Plot for the Daily Return Bitcoin Price Series.")
par(mfrow=c(1,1))

eacf(r.bitcoin)

# Some significant correlations in ACF and PACF.
# EACF table supports an ARMA(0,0) model.


#---------------------------------GARCH MODEL SPECIFICATION-----------------------------------------------
# Use squared and abs value.

abs.r.bitcoin = abs(r.bitcoin)
sq.r.bitcoin = r.bitcoin^2

#After the square transformation, there are many significant lags in the both the ACF and PACF.
#The EACF does not suggest the ARMA(0,0) model. The EACF does not look good, but from it,
#Potential models are: ARMA(2,2), ARMA(2,3), ARMA(3,3).
par(mfrow=c(1,2))
acf(sq.r.bitcoin, ci.type="ma", main = "The Sample ACF Plot for the Squared Bitcoin Price Series.")
pacf(sq.r.bitcoin, main = "The Sample PACF Plot for the Squared Bitcoin Price Series.")
par(mfrow=c(1,1))

eacf(sq.r.bitcoin)

#After the absolute value transformation, there are many significant lags in the both the ACF and PACF.
#The EACF does not suggest the ARMA(0,0) model. The EACF looks better than the squared EACF.
#Potential models are: ARMA(1,1), ARMA(1,2), ARMA(2,2).
par(mfrow=c(1,2))
acf(abs.r.bitcoin, ci.type="ma", main = "The Sample ACF Plot for the Bitcoin Price Series.")
pacf(abs.r.bitcoin, main = "The Sample PACF Plot for the Bitcoin Price Series.")
par(mfrow=c(1,1))

eacf(abs.r.bitcoin)

#The models correspond to parameter settings of [max(1,1),1], [max(1,2),2], [max(2,2),2], [max(3,3),3]
# Set of possible GARCH models are: {GARCH(1,1), GARCH(2,2), GARCH(3,3)}.

# Fit the models.
model.11 = garch(r.bitcoin, order= c(1,1), trace = FALSE)

#Remove NaN values.
model.11$residuals <- na.omit(model.11$residuals)

summary(model.11)
#All significant.
model.11_2 = garchFit(formula ~garch(1,1), data = r.bitcoin)
summary(model.11_2)


model.22 = garch(r.bitcoin, order = c(2,2), trace = FALSE)
summary(model.22)
#a2 not significant.
model.22_2 = garchFit(formula ~garch(2,2), data = bitcoin.ts)
summary(model.22_2)


model.33 = garch(r.bitcoin, order = c(3,3), trace = FALSE)
summary(model.33)
#a2, b1, b3 not significant.
model.33_2 = garchFit(formula ~garch(3,3), data = bitcoin.ts)
summary(model.33_2)


sort.score(AIC(model.11, model.22, model.33), score = "aic")

# GARCH(1,1) has lowest AIC value.

#------------------GARCH MODEL DIAGNOSTICS-----------------
#GARCH(1,1)
model.11$residuals <- na.omit(model.11$residuals)
res.model.11 <- model.11$residuals

par(mfrow=c(3,2))
plot(res.model.11,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
abline(h=0)
hist(res.model.11,main="Histogram of standardised residuals")
acf(model.11$residuals, main="ACF of standardised residuals")
pacf(res.model.11,main="PACF of standardised residuals")
qqnorm(res.model.11,main="QQ plot of standardised residuals")
qqline(res.model.11, col = 2)
print(shapiro.test(res.model.11))
k=0
LBQPlot(res.model.11, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
par(mfrow=c(1,1))


#GARCH(2,2)
model.22$residuals <- na.remove(model.22$residuals)
res.model.22 <- model.22$residuals

par(mfrow=c(3,2))
plot(res.model.22,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
abline(h=0)
hist(res.model.22,main="Histogram of standardised residuals")
acf(model.22$residuals, main="ACF of standardised residuals")
pacf(res.model.22,main="PACF of standardised residuals")
qqnorm(res.model.22,main="QQ plot of standardised residuals")
qqline(res.model.22, col = 2)
print(shapiro.test(res.model.22))
k=0
LBQPlot(res.model.22, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
par(mfrow=c(1,1))


#GARCH(3,3)
model.33$residuals <- na.remove(model.33$residuals)
res.model.33 <- model.33$residuals

par(mfrow=c(3,2))
plot(res.model.33,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
abline(h=0)
hist(res.model.33,main="Histogram of standardised residuals")
acf(model.33$residuals, main="ACF of standardised residuals")
pacf(res.model.33,main="PACF of standardised residuals")
qqnorm(res.model.33,main="QQ plot of standardised residuals")
qqline(res.model.33, col = 2)
print(shapiro.test(res.model.33))
k=0
LBQPlot(res.model.33, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
par(mfrow=c(1,1))


#------------------ARMA + GARCH MODEL FITTING------------------
#ARMA(1,1) + GARCH(1,1)
model1 <- ugarchspec(variance.model = list(model="sGARCH",
                                           garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(1,1)),
                     distribution.model = "norm")                  

m.11.11 <- ugarchfit(spec=model1, data=bitcoin.ts)
m.11.11

modelTEST <- ugarchspec(variance.model = list(model="sGARCH",
                                           garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(5,4)),
                     distribution.model = "norm")  

test.model <- ugarchfit(spec=modelTEST, data=r.bitcoin)
test.model

forc.testmodel = ugarchforecast(test.model, data = r.bitcoin, n.ahead = 10)


# View the forecast values.
forc.testmodel@forecast$seriesFor

testBack <- diffinv(r.bitcoin, xi=log(bitcoin.ts)[1])
testBack <- exp(testBack)

testBack

MASE(as.vector(btc.forecasts), as.vector(forc.11.21@forecast$seriesFor))
# 3.613293


#MA(1) not significant.


#ARMA(1,2) + GARCH(1,1)
model2 <- ugarchspec(variance.model = list(model="sGARCH",
                                           garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(1,2)),
                     distribution.model = "norm")                  

m.11.12 <- ugarchfit(spec=model2, data=bitcoin.ts)
m.11.12
#MA(1),MA(2) not significant.


model3 <- ugarchspec(variance.model = list(model="sGARCH",
                                           garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(2,1)),
                     distribution.model = "norm")                  

m.11.21 <- ugarchfit(spec=model3, data=bitcoin.ts)
m.11.21
#All significant. Increase AR component.

model4 <- ugarchspec(variance.model = list(model="sGARCH",
                                           garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(3,1)),
                     distribution.model = "norm")                  

m.11.31 <- ugarchfit(spec=model4, data=bitcoin.ts)
m.11.31
#MA(1) not significant. Try increasing AR instead. 


model5 <- ugarchspec(variance.model = list(model="sGARCH",
                                           garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(3,2)),
                     distribution.model = "norm")                  

m.11.32 <- ugarchfit(spec=model5, data=bitcoin.ts)
m.11.32
#All significant. Increase AR again.

model6 <- ugarchspec(variance.model = list(model="sGARCH",
                                           garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(4,2)),
                     distribution.model = "norm")                  

m.11.42 <- ugarchfit(spec=model6, data=bitcoin.ts)
m.11.42

#All significant.

model7 <- ugarchspec(variance.model = list(model="sGARCH",
                                           garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(5,2)),
                     distribution.model = "norm")                  

m.11.52 <- ugarchfit(spec=model7, data=bitcoin.ts)
m.11.52
#AR(4),AR(5) insignifianct.

model8 <- ugarchspec(variance.model = list(model="sGARCH",
                                           garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(4,3)),
                     distribution.model = "norm")                  

m.11.43 <- ugarchfit(spec=model8, data=bitcoin.ts, solver.control = list(tol = 1e-5))
m.11.43
#Convergence problem

#INFORMATION CRITERIA
#               m.11.21   m.11.32*   m.11.42
# Akaike        9.2457    9.2353*    9.2377
# Bayes         9.2674    9.2632*    9.2686
# Shibata       9.2457    9.2353*    9.2376
# Hannan-Quinn  9.2537    9.2456*    9.2491
m.11.32
#ARMA(3,2) + GARCH(1,1) is the best model.
plot(m.11.32,which="all")

# FORECASTING
# Compute the forecast.
forc.11.21 = ugarchforecast(m.11.21, data = bitcoin.ts, n.ahead = 10)

# Plot the forecast.
plot(forc.11.21)

# View the forecast values.
forc.11.21


MASE(as.vector(btc.forecasts), as.vector(forc.11.21@forecast$seriesFor))
# 3.613293


forc.11.32 = ugarchforecast(m.11.32, data = bitcoin.ts, n.ahead = 10)
plot(forc.11.32)
forc.11.32

MASE(as.vector(btc.forecasts), as.vector(forc.11.32@forecast$seriesFor))
# 3.671382

forc.11.42 = ugarchforecast(m.11.42, data = bitcoin.ts,n.ahead = 10)
plot(forc.11.42, which=1)
forc.11.42

MASE(as.vector(btc.forecasts), forc.11.42@forecast$seriesFor)
#3.70586

#MASE for model fitting
#ARMA(2,1) + GARCH(1,1)
MASE(as.vector(bitcoin.ts),as.vector(m.11.21@fit$fitted.values))
#1.00

fitted.values


#ARMA(3,2) + GARCH(1,1)
MASE(as.vector(bitcoin.ts),as.vector(m.11.32@fit$fitted.values))
#1.01

#ARMA(4,2) + GARCH(1,1)
MASE(as.vector(bitcoin.ts),as.vector(m.11.42@fit$fitted.values))
#1.01





#-----------WTF
modelone <- ugarchspec(variance.model = list(model="sGARCH",
                                           garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(1,1)),
                     distribution.model = "norm")                  

m.11.11one <- ugarchfit(spec=modelone, data=r.bitcoin)
m.11.11one
#AR1 and MA1 insig

modeltwo <- ugarchspec(variance.model = list(model="sGARCH",
                                             garchOrder = c(1,1)),
                       mean.model = list(armaOrder = c(2,1)),
                       distribution.model = "norm")                  

m.11.11two <- ugarchfit(spec=modeltwo, data=r.bitcoin)
m.11.11two
#AR2 insig

modelthree <- ugarchspec(variance.model = list(model="sGARCH",
                                             garchOrder = c(1,1)),
                       mean.model = list(armaOrder = c(3,1)),
                       distribution.model = "norm")                  

m.11.11three <- ugarchfit(spec=modelthree, data=r.bitcoin)
m.11.11three
#AR2, AR3 insig

modelfour <- ugarchspec(variance.model = list(model="sGARCH",
                                               garchOrder = c(1,1)),
                         mean.model = list(armaOrder = c(2,2)),
                         distribution.model = "norm")                  

m.11.11four <- ugarchfit(spec=modelfour, data=r.bitcoin)
m.11.11four
#All good.

modelfive <- ugarchspec(variance.model = list(model="sGARCH",
                                              garchOrder = c(1,1)),
                        mean.model = list(armaOrder = c(2,3)),
                        distribution.model = "norm")                  

m.11.11five <- ugarchfit(spec=modelfive, data=r.bitcoin)
m.11.11five
#All good.

modelsix <- ugarchspec(variance.model = list(model="sGARCH",
                                              garchOrder = c(1,1)),
                        mean.model = list(armaOrder = c(3,3)),
                        distribution.model = "norm")                  

m.11.11six <- ugarchfit(spec=modelsix, data=r.bitcoin)
m.11.11six
#All good.

modelseven <- ugarchspec(variance.model = list(model="sGARCH",
                                             garchOrder = c(1,1)),
                       mean.model = list(armaOrder = c(4,3)),
                       distribution.model = "norm")                  

m.11.11seven <- ugarchfit(spec=modelseven, data=r.bitcoin)
m.11.11seven
#AR4 not good. overfit?

modeleight <- ugarchspec(variance.model = list(model="sGARCH",
                                               garchOrder = c(1,1)),
                         mean.model = list(armaOrder = c(3,4)),
                         distribution.model = "norm")                  

m.11.11eight <- ugarchfit(spec=modeleight, data=r.bitcoin)
m.11.11eight
#All good.


modelnine <- ugarchspec(variance.model = list(model="sGARCH",
                                               garchOrder = c(1,1)),
                         mean.model = list(armaOrder = c(3,5)),
                         distribution.model = "norm")                  

m.11.11nine <- ugarchfit(spec=modelnine, data=r.bitcoin)
m.11.11nine
#MA5 insig. overfit?

modelten <- ugarchspec(variance.model = list(model="sGARCH",
                                              garchOrder = c(1,1)),
                        mean.model = list(armaOrder = c(4,4)),
                        distribution.model = "norm")                  

m.11.11ten <- ugarchfit(spec=modelten, data=r.bitcoin)
m.11.11ten
#All good.

# Compute the forecast.
forc.11.21one = ugarchforecast(m.11.11ten, data = bitcoin.ts, n.ahead = 10)
# Plot the forecast.
plot(forc.11.21one)
# View the forecast values.
diffinv(forc.11.21one@forecast$seriesFor, xi=log(bitcoin.ts)[1])
MASE(as.vector(btc.forecasts), as.vector(exp(diffinv(forc.11.21one@forecast$seriesFor)))
# 3.613293