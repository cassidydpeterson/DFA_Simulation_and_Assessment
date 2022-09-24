# Index rescaling minimization code from Nicholas Ducharme-Barth
# April 2022

optim_fn = function(x,data)
{
  c = exp(x)
  dat = data * c
  datL=log(dat)
  y.bar = apply(datL, 1, mean, na.rm=TRUE) # get mean of each index; ENSURE y.bar >0
  dat.dm = (datL-y.bar) / y.bar  # get demeaned index
  gsd = sd(c(dat.dm), na.rm=T) # get GSD or global standard deviation for all indices
  # ensure gsd<<1; if gsd>1, try new vector c
  dat.z = dat.dm/gsd
  
  criteria_1 = ifelse(gsd<0.5,0,1)
  criteria_2 = ifelse(y.bar>0,0,1)
  criteria_3 = ifelse(datL>0,0,1)
  criteria_4 = sum((apply(dat.z, 1, sd, na.rm=T) - 1)^2)
  
  return(c(criteria_1,criteria_2,criteria_3,criteria_4))
}

set.seed(123)
constrained_fit = rgenoud::genoud(fn=optim_fn,nvars=3,starting.values=c(5,5,5),lexical=TRUE,data.type.int=FALSE,data=dat.a)
c = exp(constrained_fit$par)


# In the above code I have 3 indices so nvars and starting.values would need to be changed for more/less time series.
