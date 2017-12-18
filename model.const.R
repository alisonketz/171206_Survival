
    model{
    
    # Priors
    llambda  ~ dnorm(0, 1.0E-6) # log hazard rate
    
    # Likelihood for the total hazard
    for (j in 1:records) {
        for (k in left[j]:(right[j]-1)) {
            UCH[j,k] <- exp(llambda)
        }
        # total prob of surviving
        SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
        censor[j] ~ dbern(SLR[j])      
        }
    }

