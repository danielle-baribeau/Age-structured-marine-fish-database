# Generating the F (fishing mortality) and M (age) time-series used in the following publication : 
#Charbonneau, J. A., Keith, D. M., MacNeil, M.A., and Hutchings, J. A. 2022. Effects of fishing mortality on the age structure of marine fishes. Submitted to Canadian Journal of Fisheries and Aquatic Sciences.

# General information  ----

#1: Fishing mortality - Specify the mean and variance of this process
## mn: The mean exploitation rate in the fishery (on instantaneous scale)
## var: The variance

# Next your age time series

#2: Age - Specify the mean and variance of this process
## mn: The mean age of the populations
## var: The variance 

#3: ARIMA - Specify parameters 
## proc: 'arima' runs full 'arima' model with p/d/q as specified.  'normal' is arima with p=q=d=0. 'rand_walk" is arima with p =1 , q=d=0
## p: ARIMA ar() term
## d: ARIMA degree of differencing term
## q: ARIMA moving average term
## cor: The correlation coefficient between fm and age time series -1 <= cor <= 1
## lag: The lag between the fm and age time series, negative indicates that fm leads age, positive is age leading fm.

#4: ts  Specify Time series length 
## years: Number of years in the ARIMA time series that you will be sampled (default = 1000)
## start.year: The first year from the  ARIMA time series to sample from (avoid going too close to the start)
## final.ts.len:  Length of final time series.
# Simulation function  ----

### 1. Generating fm and mn series.
sim.fun <- function(fm = list(mn = 0.5,var=0.1),
                    age = list(mn =2,var=0.1),
                    arima = list(proc = 'arima',p=0,d=1,q=0,cor = -0.5,lag = -2),
                    ts = list(years=1000,start.year =200, final.ts.len = 20))
  
{

require("mvtnorm")
  # Covariance structure for the 'error' (log scale)
  # Get the correlation matrix   
  corr.mat <- matrix(c(1,arima$cor,arima$cor,1),nrow = 2)
  # And convert to a covariance matrix.
  cov.mat <- corr.mat * sqrt(fm$var) * sqrt(age$var)
  n <- ts$years
  # This gets us two correlated 'error' time series
  eps <- mvtnorm::rmvnorm(n = n, mean = log(c(fm$mn,age$mn)), sigma = cov.mat)
  
### 2. Specifing parameters 
  if(arima$proc == 'normal') 
  {
    fm.ts <- arima.sim(model =list(0,0,0),n = n, innov = eps[,1])
    age.ts <- arima.sim(n = n, model = list(0,0,0), innov = eps[,2])  
  }
  # If you wanted it to be an ARIMA model with specfied p/q/d terms
  if(arima$proc == 'arima')  
  {
    fm.ts <- arima.sim(model =list(arima$p,arima$d,arima$q),n = n, innov = eps[,1])
    age.ts <- arima.sim(n = n, model = list(arima$p,arima$d,arima$q), innov = eps[,2]) 
  }
  # If you wanted it to be a random walk model (p=1)
  if(arima$proc == 'rand_walk')
  {
    fm.ts <- arima.sim(model =list(1,0,0),n = n, innov = eps[,1])
    age.ts <- arima.sim(n = n, model = list(1,0,0), innov = eps[,2]) 
  }
  
### 3. Save the time series (f input, age output) - (LL, fishing scenraio)
dat <- data.frame(year = ts$start.year:(ts$start.year+ts$final.ts.len-1),
                    fm.ts = exp(fm.ts[ts$start.year:(ts$start.year+ts$final.ts.len-1)]),
                    age.ts = exp(age.ts[(ts$start.year-arima$lag):(ts$start.year+ts$final.ts.len-1-arima$lag)]))
  return(dat)
} 

# Example with mod 1  ----

fm.mn<-1.56
age.mn<-5.57
p<-0
d<-1
q<-0
proc<-"arima"
sim.years <- 200 
syear <- 100 
ts.len <- 38 # chg for last 4 sim
n.sims <- 10000 
obs.cor <- NA
dat <- NULL

fm.var <- 0.01
age.var <- 0.03
lag <- -4
cor <- -0.9

for(i in 1:n.sims) 
{
  # Run the simulations
  dat[[i]] <- sim.fun(fm = list(mn = fm.mn,var=fm.var),
                      age = list(mn = age.mn,var=age.var),
                      arima= list(proc = proc,p=p,d=d,q=q,cor =cor,lag = lag),
                      ts = list(years=sim.years,start.year =syear, final.ts.len = ts.len))
  tmp <- ccf(dat[[i]]$age.ts,dat[[i]]$fm.ts,plot=F)
  obs.cor[i] <- tmp$acf[tmp$lag == lag]
}
# Checking results : obs.cor should be the cor you specified with some uncertainty  
hist(obs.cor)
mean(obs.cor)
require(ggplot2)
ggplot(dat[[1]]) + geom_line(aes(x=year,y = age.ts)) + geom_line(aes(x=year,y=fm.ts)) + scale_y_log10()


