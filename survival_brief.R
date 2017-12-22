
library(rjags)
library(dclone)
library(doParallel)
library(foreach)
library(coda)
library(Hmisc)
library(lubridate)
library(xlsx)
library(gridExtra)


###
### set working directory
###

setwd('/home/aketz/Documents/Survival/171206_Survival')

###
### model statement
###


sink("model.cwd.R")
cat("
    model{
    
    # Priors
    for(j in 1:records){
    llambda[j] <- beta0 + beta1*x1[j]+beta2*x2[j]+beta3*x3[j]
    }
    
    beta0 ~ dnorm(0, 1.0E-6) # intercept log hazard rate
    beta1 ~ dnorm(0, 1.0E-6) # cwd effect log hazard rate
    beta2 ~ dnorm(0, 1.0E-6) # bow harvest effect log hazard rate
    beta3 ~ dnorm(0, 1.0E-6) # gun harvest effect log hazard rate
    
    # Likelihood for the total hazard
    for (j in 1:records) {
    for (k in left[j]:(right[j]-1)) {
    UCH[j,k] <- exp(llambda[j])
    }
    # total prob of surviving
    SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censor[j] ~ dbern(SLR[j])      
    }
    
    #Derived parameters
    for (t in 1:n.weeks){
    llambda.out[t,1]<-beta0 + beta2*x4[t] +beta3*x5[t]#CWD-
    llambda.out[t,2]<-beta0 + beta1 + beta2*x4[t] +beta3*x5[t] #CWD+
    for(j in 1:2){
    UCH0[t,j]<-exp(llambda.out[t,j])
    CH0[t,j]<-sum(UCH0[1:t,j])
    S0[t,j]<-exp(-CH0[t,j])
    }
    }
    }
    ",fill = TRUE)
sink()

#specify initial values
inits<-list(list("beta0"=1,"beta1"=-1.5,"beta2"=1,"beta3"=.1),list("beta0"=-.5,"beta1"=.25,"beta2"=2,"beta3"=-.3),list("beta0"=1,"beta1"=.5,"beta2"=.5,"beta3"=1.6))

#identify params to monitor
parameters<-c("beta0","beta1","beta2","beta3","S0") 

###
### x1 = cwd affect <- beta1 <- N.temp[,5].... N.data.fit[,4]
### x2 = bow effect <- beta2 <- N.temp[,9] ... N.data.fit[,8]
### x3 = gun effect <- beta3 <- N.temp[,7] ... N.data.fit[,6]
### x4 = for derived parameters, bow season survival
### x5 = for derived parameters, gun season survival
###




# Bundle data
# jags.data <- list(records=n.fit,
#                   left=N.data.fit[,1],
#                   right=N.data.fit[,2],
#                   censor=N.data.fit[,3],
#                   x1=N.data.fit[,4],
#                   x2=N.data.fit[,8],
#                   x3=N.data.fit[,6],
#                   n.weeks=study.end,
#                   x4=bow.harvest.haz,
#                   x5=gun.harvest.haz)

# save(jags.data,file="jags.data.R")
load("jags.data.Rdata")

# MCMC settings
nt = 1
ni = 1000
nb = 10000
nc = 3

# call parallel version of jags using dclone
cl=makeCluster(3)
out.cwd=jags.parfit(cl, data=jags.data, params=parameters, model="model.cwd.R",inits=inits,n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
stopCluster(cl)

n.weeks=max(N.data.fit,na.rm=TRUE)

fit.sum=summary(out.cwd)[[1]]
fit.quant=summary(out.cwd)[[2]]
fit.sum
fit.quant
