
#Certain-cause

## @knitr model.certain
model1 <- function(){
  
  ###Model - Cause of death assigned with certainty
  
  # Priors
  llambda  ~ dnorm(0, 1.0E-6) # log hazard rate
  
  # for(i in 1:4){
  #   alpha[i]<-1   
  # }
  # 
  # p ~ ddirich(alpha[])  #true probability of cause of death
  
  # Likelihood for the total hazard
  for (j in 1:records) {
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp(llambda)
    }
    # total prob of surviving
    SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censor[j] ~ dbern(SLR[j])      
  }
  
  # 
  # 
  # # Likelihood cause of death
  # for (i in 1:c.records) {
  #   cause[i] ~ dcat(p[])
  # }  
  # 
  # #Derived parameters
  # for (i in 1:4) {
  #   causehaz[i] <- exp(llambda)*p[i] #cause-specific hazards
  # }               
}


#Prior-cause

## @knitr model.uncertain
model2<-function(){
  # Priors
  llambda  ~ dnorm(0, 1.0E-6) # log hazard rate
  
  for(i in 1:4){
    alpha[i]<-1   
  }
  for(i in 1:c.records){
    cause[i] ~ dcat(p.obs[i,1:4])
  }
  
  #true probability of cause of death
  p ~ ddirich(alpha[])  
  
  # Likelihood for the total hazard
  for (j in 1:records) {
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp(llambda)
    }
    
    # total prob of surviving
    SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censor[j] ~ dbern(SLR[j])      
  }
  
  
  # Likelihood cause of death
  for (i in 1:c.records) {
    phi[i] <--log(p[cause[i]])
    
    #zeros trick since cause[i] given prior can't give ~ again
    z[i] ~ dpois(phi[i])  
  }  
  
  #Derived parameters
  for (i in 1:4) {
    causehaz[i] <- exp(llambda)*p[i] #cause-specific hazards
  }               
}  



#specify initial values
init1<-function(){list("llambda"=-3, p=rep(0.25,4))}

#identify parms to monitor
parms<-c("llambda","p","causehaz") 