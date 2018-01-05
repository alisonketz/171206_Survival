
    model{
    beta0 ~ dnorm(0, 1.0E-6) # intercept log hazard rate
    beta1 ~ dnorm(0, 1.0E-6) # cwd effect log hazard rate
    beta2 ~ dnorm(0, 1.0E-6) # bow harvest effect log hazard rate
    beta3 ~ dnorm(0, 1.0E-6) # gun harvest effect log hazard
    beta4 ~ dnorm(0, 1.0E-6) # sex effect log hazard

    # Zero-inflated sensitivity of RAMALT
    for(j in 1:records){
        temp[j] <-p*z[j]
        x1[j] ~ dbin(temp[j],1)
        z[j] ~ dbin(psi,1)
    }

    p ~ dbeta(92.5,44.5)#based on sensitivity of Thomsen et al 2012
    psi~dbeta(1,1)

    # Likelihood for the total hazard
    for (j in 1:records) {
        for (k in left[j]:(right[j]-1)) {
             UCH[j,k] <- exp(beta0 + beta1*z[j]+beta2*x4[k]+beta3*x5[k]+beta4*x2[j])
        }
        # total prob of surviving
        SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
        censor[j] ~ dbern(SLR[j])      
    }

    #Derived parameters
    for (t in 1:n.weeks){
        llambda.out[t,1]<-beta0 + beta2*x4[t] +beta3*x5[t] #Female CWD-
        llambda.out[t,2]<-beta0 + beta2*x4[t] +beta3*x5[t] + beta4 #Male CWD -
        llambda.out[t,3]<-beta0 + beta1 + beta2*x4[t] + beta3*x5[t] #Female CWD+
        llambda.out[t,4]<-beta0 + beta1 + beta2*x4[t] + beta3*x5[t] + beta4 #Male CWD +
        for(j in 1:4){
            UCH0[t,j]<-exp(llambda.out[t,j])
            CH0[t,j]<-sum(UCH0[1:t,j])
            S0[t,j]<-exp(-CH0[t,j])
        }
    }
}
    
