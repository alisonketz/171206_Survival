
    model{
    
    # Priors
    for(j in 1:records){
        llambda[j] <- beta0 + beta1*x1[j]+beta2*x2[j]
    }
    
    beta0 ~ dnorm(0, 1.0E-6) # intercept log hazard rate
    beta1 ~ dnorm(0, 1.0E-6) # cwd effect log hazard rate
    beta2 ~ dnorm(0, 1.0E-6) # hunter harvest effect log hazard rate

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
        llambda.out[t,1]<-beta0
        llambda.out[t,2]<-beta0+beta1
        llambda.out[t,3]<-beta0+beta2
        llambda.out[t,4]<-beta0+beta1+beta2
        for(j in 1:4){
            UCH0[t,j]<-exp(llambda.out[t,j])
            CH0[t,j]<-sum(UCH0[1:t,j])
            S0[t,j]<-exp(-CH0[t,j])
        }
    }
}
    
