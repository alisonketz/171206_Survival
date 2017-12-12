###
### Survival Analysis 
### A.Ketz 12/6/2017
###
### Note that this must be run in 32-bit R because RODBC can only work in 32bit.
###
###

###
### Preliminaries
###

rm(list=ls())

# setwd('F:/171206_Survival')

setwd('/home/aketz/Documents/Survival/171206_Survival')

library(rjags)
library(dclone)
library(doParallel)
library(foreach)
library(coda)
library(Hmisc)
library(lubridate)


###
### Load Data and fix column names
###

#load Access Database into R using Hmisc package with the function mdb.get (must use this on Linux)
# database =mdb.get('~/Documents/Data/SWDPPdeerDB.MDB')

d.cap=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Captures_adult_form_2017")
d.fawncap=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Fawn Capture")
d.mort=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Mortalities")
d.cens=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Censor")
d.cap.ramalt=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "CapRAMALT")


names(d.fawncap)=tolower(gsub("[[:punct:]]","",names(d.fawncap)))
names(d.cap)=c("InDB","DateEntered","CaptureStatus","LowTag","Date","Datalogger","Fieldcrew", "StudyArea","Location","Lat", "Long", "TrapType", "TrapNumber","Drughandler", "Timeofinitialstress", "Timeofinitialinjection", "DoseBAM","Bottle", "Timeofinduction","Collar", "Frequency" , "CollarID", "Magnetremoved" ,"Earpunch", "EartagLeft","EartagRight", "Ageclass", "Confidenceinage", "Sex","Hindlegcm", "Girthcm","Fecalsample", "RAMALT","RAMALTlocation", "Samplequality","RAMALTsamplerlastname", "Toothextraction","Bodycondition", "Bodyconditionvalue","Ultrasound", "UltrasoundPicID","VIT", "VITfrequency","Pregnant", "Xoffetuses","Weight", "ReversalAtipamezoleDose","Atipamezolebottle", "ReversalTime","TimeofRecovery", "Observations","EnteredCheck", "Mortality","Censor")
names(d.mort)=c("InDB","Lowtag","Technicians","DataRecorder","Lefteartag","Lefteartag","Righteartag","Righteartag","Collar","Collar","Sex","Captured","MORTAlertdate","CollarFound", "EstMortDate","Comment1","PropertyType","Landtype","County","Property", "Photos","Incisor","Lymphnodes","CWDSuspect","FieldNecropsy","LabNecropsy","SecondaryFateSheet","Cause1", "Cause1wt","Cause2","Cause2wt","Cause3","Cause3wt","Comments2","HunterName","Huntingseason","Weapon","Lat","Long","CustomerID","ConfidenceinCause","NumberAntlerPoints","Township","ApproximateLocationofHarvest","NumberDeerInGroup","Howfarawaywasdeerwhenshot", "Emailaddress","CWDsampleID","Legmeasurment","Girth","Dressedweight","Hunterphone","HunterAddress","Frequency")
names(d.cens)=c("Lowtag", "Lefteartag","Righteartag","Sex","AgeAtCapture","Latitude","Longitude","County","Censorcause","CensWt","Othercause1","Othercause1wt","Othercause2","Othercause2wt","Cencomments","CaptureDate","CollarRecovered","PropertyCollarRecovered","Landtype","Frequency","CensorDate","Technicians","CollarNumber")

for(i in 1:length(names(d.cap))){names(d.cap)[i]=paste(tolower(substr(names(d.cap)[i], 1, 1)),substr(names(d.cap)[i], 2, nchar(names(d.cap)[i])),sep="")}
for(i in 1:length(names(d.mort))){names(d.mort)[i]=paste(tolower(substr(names(d.mort)[i], 1, 1)),substr(names(d.mort)[i], 2, nchar(names(d.mort)[i])),sep="")}
for(i in 1:length(names(d.cens))){names(d.cens)[i]=paste(tolower(substr(names(d.cens)[i], 1, 1)),substr(names(d.cens)[i], 2, nchar(names(d.cens)[i])),sep="")}



###
### EDA of capture data
###


# d.cap$mortality
# names(d.cap)
# 
# sum(d.cap$censor)
# d.cap[which(as.logical(d.cap$censor)),1:5]
# table(d.cap$pregnant)
# names(d.mort)
# 
# head(d.mort[1:15])
# d.mort$mORTAlertdate
# 

### 
### Dates
### Initial date of captures, ie. start of the study
### 01/07/2017
### formatting date vectors
###

d.cap$date = as.character(d.cap$date)
d.cap$date=strtrim(d.cap$date,nchar(d.cap$date)-9)
d.cap$date = as.Date(d.cap$date,format=c('%m/%d/%y'))
d.cap$dateEntered = as.character(d.cap$dateEntered)
d.cap$dateEntered=strtrim(d.cap$dateEntered,nchar(d.cap$dateEntered)-9)
d.cap$dateEntered = as.Date(d.cap$dateEntered,format=c('%m/%d/%y'))

d.fawncap$date = as.character(d.fawncap$date)
d.fawncap$date=strtrim(d.fawncap$date,nchar(d.fawncap$date)-9)
d.fawncap$date = as.Date(d.fawncap$date,format=c('%m/%d/%y'))
d.fawncap$dateentered = as.character(d.fawncap$dateentered)
d.fawncap$dateentered=strtrim(d.fawncap$dateentered,nchar(d.fawncap$dateentered)-9)
d.fawncap$dateentered = as.Date(d.fawncap$dateentered,format=c('%m/%d/%y'))

d.mort$captured = as.character(d.mort$captured)
d.mort$captured=strtrim(d.mort$captured,nchar(d.mort$captured)-9)
d.mort$captured = as.Date(d.mort$captured,format=c('%m/%d/%y'))
d.mort$mORTAlertdate = as.character(d.mort$mORTAlertdate)
d.mort$mORTAlertdate[d.mort$mORTAlertdate==""]=NA
d.mort$mORTAlertdate[!is.na(d.mort$mORTAlertdate)]=strtrim(d.mort$mORTAlertdate[!is.na(d.mort$mORTAlertdate)],nchar(d.mort$mORTAlertdate[!is.na(d.mort$mORTAlertdate)])-9)
d.mort$mORTAlertdate = as.Date(d.mort$mORTAlertdate,format=c('%m/%d/%y'))
d.mort$estMortDate = as.character(d.mort$estMortDate)
d.mort$estMortDate[d.mort$estMortDate==""]=NA
d.mort$estMortDate[!is.na(d.mort$estMortDate)]=strtrim(d.mort$estMortDate[!is.na(d.mort$estMortDate)],nchar(d.mort$estMortDate[!is.na(d.mort$estMortDate)])-9)
d.mort$estMortDate = as.Date(d.mort$estMortDate,format=c('%m/%d/%y'))

d.cens$censorDate = as.character(d.cens$censorDate)
d.cens$censorDate=strtrim(d.cens$censorDate,nchar(d.cens$censorDate)-9)
d.cens$censorDate = as.Date(d.cens$censorDate,format=c('%m/%d/%y'))

d.cens$captureDate = as.character(d.cens$captureDate)
d.cens$captureDate=strtrim(d.cens$captureDate,nchar(d.cens$captureDate)-9)
d.cens$captureDate = as.Date(d.cens$captureDate,format=c('%m/%d/%y'))


library(lubridate)
start.week=week(ymd('2017-01-09'))

###
### Organizing data into processable form
###

n.cap=dim(d.cap)[1]
n.mort = dim(d.mort)[1]
n.cens = dim(d.cens)[1]
n.fawncap=dim(d.fawncap)[1]


# cap.ind = d.cap$lowTag
# mort.ind = d.mort$lowtag
# cens.ind = d.cens$lowtag
# fawncap.ind=d.fawncap$lowtag

###
### initial format of data - rows = individuals, cols = (e_i,r_i,s_i,censor)
###

#making dates continuous
d.cap$dateCts = (week(d.cap$date) - start.week + 1)+(wday(d.cap$date)%%7)/7
d.mort$estMortDateCts = (week(d.mort$estMortDate) - start.week + 1)+(wday(d.mort$estMortDate)%%7)/7
d.cens$censorDateCts = (week(d.cens$censorDate) - start.week + 1)+(wday(d.cens$censorDate)%%7)/7
d.fawncap$date = (week(d.fawncap$date) - start.week + 1)+(wday(d.fawncap$date)%%7)/7


###
### lumping adults and fawns together
###



#initialize matrix of captures/recaptures/morts
N.temp = matrix(NA,nr=n.cap+n.fawncap,nc = 4)

#initialize with all animals alive
N.temp[,4] = 1


for(i in 1:n.cap){
    N.temp[i,1] = d.cap$dateCts[i]
    for(j in 1:n.mort){
        if(d.cap$lowTag[i]==d.mort$lowtag[j]){
            N.temp[i,3] = d.mort$estMortDateCts[j]
            N.temp[i,2] = d.mort$estMortDateCts[j]-1/7
            N.temp[i,4] = 0
        }
    }
}


for(i in 1:n.fawncap){
    N.temp[n.cap+i,1] = d.cap$dateCts[i]
    for(j in 1:n.mort){
        if(d.fawncap$lowtag[i]==d.mort$lowtag[j]){
            N.temp[i+n.cap,3] = d.mort$estMortDateCts[j]
            N.temp[i+n.cap,2] = d.mort$estMortDateCts[j]-1/7
            N.temp[i+n.cap,4] = 0
        }
    }
}


# 
# for(k in 1:n.cens){
#     if(d.cap$lowTag[i]==d.cens$lowtag[k]){
#         N.temp[i,2] = Inf
#     }
# }



N.temp
for(i in 1:(n.cap+n.fawncap)){
    if(is.na(N.temp[i,2]))N.temp[i,2]=Inf
    if(is.na(N.temp[i,3]))N.temp[i,3]=Inf
}

N.temp

###
### Format data matrix to fit into jags
###

N.data.fit = matrix(NA,nr=n.cap+n.fawncap+n.mort,ncol=3)

indx = 1
for(i in 1:(n.cap+n.fawncap)){
    if(N.temp[i,4]==1){
        N.data.fit[indx,] = c(N.temp[i,1:2],1)
        indx=indx+1
    }
    if(N.temp[i,4]==0){
        N.data.fit[indx,] = c(N.temp[i,1:2],1)
        indx=indx+1
        N.data.fit[indx,] = c(N.temp[i,2:3],0)
        indx=indx+1
    }
}

N.data.fit

######################################################################################################################3
###
### Separating adults and fawns -- in the mortality data
###
######################################################################################################################3

#removing fawns from the mortality

mort.rm=c()
for(l in 1:n.fawncap){
    mort.rm=c(mort.rm,which(d.mort$lowtag==d.fawncap$lowtag[l]))
}

d.mort.ad=d.mort[-mort.rm,]

n.mort = dim(d.mort.ad)[1]

#initialize matrix of captures/recaptures/morts
N.temp = matrix(NA,nr=n.cap,nc = 4)

#initialize with all animals alive
N.temp[,4] = 1


for(i in 1:n.cap){
    N.temp[i,1] = floor(d.cap$dateCts[i])
    for(j in 1:n.mort){
        if(d.cap$lowTag[i]==d.mort.ad$lowtag[j]){
            N.temp[i,3] = floor(d.mort.ad$estMortDateCts[j])
            N.temp[i,2] = floor(d.mort.ad$estMortDateCts[j])-1
            N.temp[i,4] = 0
        }
    }
}
N.temp

non.harvest.survival.end=week(ymd('2017-11-18'))
harvest.start.week=week(ymd('2017-11-18'))+1
harvest.start.week

N.temp
for(i in 1:(n.cap+n.fawncap)){
    if(is.na(N.temp[i,2]))N.temp[i,2]=non.harvest.survival.end
    if(is.na(N.temp[i,3]))N.temp[i,3]=non.harvest.survival.end
}

N.temp


###
### Format data matrix to fit into jags
###

N.data.fit = matrix(NA,nr=n.cap+n.mort-1,ncol=3)

indx = 1
for(i in 1:n.cap){
    if(N.temp[i,4]==1){
        N.data.fit[indx,] = c(N.temp[i,1:2],1)
        indx=indx+1
    }
    if(N.temp[i,4]==0){
        N.data.fit[indx,] = c(N.temp[i,1:2],1)
        indx=indx+1
        N.data.fit[indx,] = c(N.temp[i,2:3],0)
        indx=indx+1
    }
}

N.data.fit


###
### define Jags model 
###

sink("model.const.R")
cat("
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
",fill = TRUE)
sink()

#specify initial values
inits<-list(list("llambda"=-1),list("llambda"=-2),list("llambda"=-3))

#identify params to monitor
parameters<-c("llambda") 

# Bundle data
jags.data <- list(records=dim(N.data.fit)[2],left=N.data.fit[,1],right=N.data.fit[,2],censor=N.data.fit[,3])



# MCMC settings
nt = 1
ni = 1000
nb = 10000
nc = 3

# call parallel version of jags using dclone
cl=makeCluster(3)
out.const=jags.parfit(cl, data=jags.data, params=parameters, model="model.const.R",inits=inits,n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
stopCluster(cl)

summary(out.const)



######################################################################################################################3
###
### Separating adults and fawns -- in the mortality data
### Adding covariates
###
######################################################################################################################3

#removing fawns from the mortality

mort.rm=c()
for(l in 1:n.fawncap){
    mort.rm=c(mort.rm,which(d.mort$lowtag==d.fawncap$lowtag[l]))
}

d.mort.ad=d.mort[-mort.rm,]

n.mort = dim(d.mort.ad)[1]

#initialize matrix of captures/recaptures/morts
N.temp = matrix(NA,nr=n.cap,nc = 6)



#initialize with all animals alive
N.temp[,4] = 1

for(i in 1:n.cap){
    N.temp[i,1] = floor(d.cap$dateCts[i])
    N.temp[i,4] = d.cap$bodyconditionvalue[i]
    N.temp[i,5] = d.cap$
    for(j in 1:n.mort){
        if(d.cap$lowTag[i]==d.mort.ad$lowtag[j]){
            N.temp[i,3] = floor(d.mort.ad$estMortDateCts[j])
            N.temp[i,2] = floor(d.mort.ad$estMortDateCts[j])-1
            N.temp[i,4] = 0
        }
    }
}





non.harvest.survival.end=week(ymd('2017-11-18'))
harvest.start.week=week(ymd('2017-11-18'))+1
harvest.start.week

N.temp
for(i in 1:(n.cap+n.fawncap)){
    if(is.na(N.temp[i,2]))N.temp[i,2]=non.harvest.survival.end
    if(is.na(N.temp[i,3]))N.temp[i,3]=non.harvest.survival.end
}

N.temp


###
### Format data matrix to fit into jags
###

N.data.fit = matrix(NA,nr=n.cap+n.mort-1,ncol=3)

indx = 1
for(i in 1:n.cap){
    if(N.temp[i,4]==1){
        N.data.fit[indx,] = c(N.temp[i,1:2],1)
        indx=indx+1
    }
    if(N.temp[i,4]==0){
        N.data.fit[indx,] = c(N.temp[i,1:2],1)
        indx=indx+1
        N.data.fit[indx,] = c(N.temp[i,2:3],0)
        indx=indx+1
    }
}

N.data.fit


###
### define Jags model 
###

sink("model.const.R")
cat("
    model{
    
    # Priors
    for(i in 1:records){
        cloglog(llambda[i]) <- beta0 + beta1*x1[i] + beta2*x2[i]
    }


    beta0 ~ dnorm(0, 1.0E-6) # intercept log hazard rate
    beta1 ~ dnorm(0, 1.0E-6) # cwd effect log hazard rate

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
    ",fill = TRUE)
sink()

#specify initial values
inits<-list(list("llambda"=-1),list("llambda"=-2),list("llambda"=-3))

#identify params to monitor
parameters<-c("llambda") 

# Bundle data
jags.data <- list(records=dim(N.data.fit)[2],left=N.data.fit[,1],right=N.data.fit[,2],censor=N.data.fit[,3])



# MCMC settings
nt = 1
ni = 1000
nb = 10000
nc = 3

# call parallel version of jags using dclone
cl=makeCluster(3)
out.const=jags.parfit(cl, data=jags.data, params=parameters, model="model.const.R",inits=inits,n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
stopCluster(cl)

summary(out.const)



