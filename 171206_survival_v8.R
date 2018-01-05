###
### Survival Analysis 
### A.Ketz 12/6/2017
###
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
library(xlsx)
library(gridExtra)
library(xtable)

###
### Load Previous Run
### 

load('survival_v8.Rdata')

###
### Load Data and fix column names
###

#load Access Database into R using Hmisc package with the function mdb.get (must use this on Linux)
# database =mdb.get('~/Documents/Data/SWDPPdeerDB.MDB')

d.cap=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Adult_Capture_2017_2018")
d.fawncap=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Fawn Capture")
d.mort=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Mortalities")
d.cens=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Censor")
d.cwd=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "RAMALT")
d.post.cwd=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Postmortem_CWD")

names(d.fawncap)=tolower(gsub("[[:punct:]]","",names(d.fawncap)))
names(d.cap)=c("InDB","DateEntered","CaptureStatus","lowtag","Date","Datalogger","Fieldcrew", "StudyArea","Location","Lat", "Long", "TrapType", "TrapNumber","Drughandler", "Timeofinitialstress", "Timeofinitialinjection", "DoseBAM","Bottle", "Timeofinduction","Collar", "Frequency" , "CollarID", "Magnetremoved" ,"Earpunch", "EartagLeft","EartagRight", "Ageclass", "Confidenceinage", "Sex","Hindlegcm", "Girthcm","Fecalsample", "RAMALT","RAMALTlocation", "Samplequality","RAMALTsamplerlastname", "Toothextraction","Bodycondition", "Bodyconditionvalue","Ultrasound", "UltrasoundPicID","VIT", "VITfrequency","Pregnant", "Xoffetuses","Weight", "ReversalAtipamezoleDose","Atipamezolebottle", "ReversalTime","TimeofRecovery", "Observations","EnteredCheck", "Mortality","Censor")
names(d.mort)=c("InDB","Lowtag","Technicians","DataRecorder","Lefteartag","Lefteartag","Righteartag","Righteartag","Collar","Collar","Sex","Captured","MORTAlertdate","CollarFound", "EstMortDate","Comment1","PropertyType","Landtype","County","Property", "Photos","Incisor","Lymphnodes","CWDSuspect","FieldNecropsy","LabNecropsy","SecondaryFateSheet","Cause1", "Cause1wt","Cause2","Cause2wt","Cause3","Cause3wt","Comments2","HunterName","Huntingseason","Weapon","Lat","Long","CustomerID","ConfidenceinCause","NumberAntlerPoints","Township","ApproximateLocationofHarvest","NumberDeerInGroup","Howfarawaywasdeerwhenshot", "Emailaddress","CWDsampleID","Legmeasurment","Girth","Dressedweight","Hunterphone","HunterAddress","Frequency")
names(d.cens)=c("Lowtag", "Lefteartag","Righteartag","Sex","AgeAtCapture","Latitude","Longitude","County","Censorcause","CensWt","Othercause1","Othercause1wt","Othercause2","Othercause2wt","Cencomments","CaptureDate","CollarRecovered","PropertyCollarRecovered","Landtype","Frequency","CensorDate","Technicians","CollarNumber")
names(d.cwd)=tolower(gsub('[[:punct:]]',"",names(d.cwd)))
names(d.post.cwd)=tolower(gsub('[[:punct:]]',"",names(d.post.cwd)))


for(i in 1:length(names(d.cap))){names(d.cap)[i]=paste(tolower(substr(names(d.cap)[i], 1, 1)),substr(names(d.cap)[i], 2, nchar(names(d.cap)[i])),sep="")}
for(i in 1:length(names(d.mort))){names(d.mort)[i]=paste(tolower(substr(names(d.mort)[i], 1, 1)),substr(names(d.mort)[i], 2, nchar(names(d.mort)[i])),sep="")}
for(i in 1:length(names(d.cens))){names(d.cens)[i]=paste(tolower(substr(names(d.cens)[i], 1, 1)),substr(names(d.cens)[i], 2, nchar(names(d.cens)[i])),sep="")}

### Double check loading correctly

# head(d.cap)
# head(d.mort)
# head(d.cens)
# head(d.fawncap)
# head(d.cwd)

###
### Number individuals within each dataframe
###

n.cap=dim(d.cap)[1]
n.cens = dim(d.cens)[1]
n.fawncap=dim(d.fawncap)[1]
n.mort = dim(d.mort)[1]
n.cwdtest=dim(d.cwd)[1]
n.postcwd=dim(d.post.cwd)[1]

#change class of CWD tests to integer for d.cwd
class(d.cwd$lowtag)="integer"

#add CWD status to capture dataframe
d.cap$cwdstatus=rep(NA,n.cap)
for(i in 1:n.cap){
    for(j in 1:n.cwdtest){
        if(d.cap$lowtag[i]==d.cwd$lowtag[j]){
            d.cap$cwdstatus[i]=as.character(d.cwd$result[j])
        }
    }
    if(d.cap$cwdstatus[i]=='ISF' | is.na(d.cap$cwdstatus[i]))d.cap$cwdstatus[i]='Negative'
}
d.cap$cwdstatus=as.factor(d.cap$cwdstatus)
d.cap$cwdstatus=as.numeric(d.cap$cwdstatus)-1

###
### Adding the december gun harvested individual
###

empty=matrix(NA,nr=1,nc=dim(d.mort)[2])
colnames(empty)=colnames(d.mort)
d.mort=rbind(d.mort,empty)
n.mort=dim(d.mort)[1]

d.mort$lowtag[n.mort]=6843
d.mort$cause1[n.mort]="Hunter harvest"
d.mort$weapon[n.mort]="Rifle"
d.mort$estMortDate=as.character(d.mort$estMortDate)
d.mort$estMortDate[n.mort]="12/24/17 00:00:00"

d.mort[dim(d.mort)[1],]

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

d.mort$collarFound= as.character(d.mort$collarFound)
d.mort$collarFound[d.mort$collarFound!=""]=strtrim(d.mort$collarFound[d.mort$collarFound!=""],nchar(d.mort$collarFound[d.mort$collarFound!=""])-9)
d.mort$collarFound= as.Date(d.mort$collarFound,format=c('%m/%d/%y'))

d.cens$censorDate = as.character(d.cens$censorDate)
d.cens$censorDate=strtrim(d.cens$censorDate,nchar(d.cens$censorDate)-9)
d.cens$censorDate = as.Date(d.cens$censorDate,format=c('%m/%d/%y'))

d.cens$captureDate = as.character(d.cens$captureDate)
d.cens$captureDate=strtrim(d.cens$captureDate,nchar(d.cens$captureDate)-9)
d.cens$captureDate = as.Date(d.cens$captureDate,format=c('%m/%d/%y'))


#start of the study
start.week=week(ymd('2017-01-09'))


###
### removing fawns from the mortality
###

mort.rm=c()
for(l in 1:n.fawncap){
    mort.rm=c(mort.rm,which(d.mort$lowtag==d.fawncap$lowtag[l]))
}
d.mort.ad=d.mort[-mort.rm,]

n.mort=dim(d.mort.ad)[1]



###
### initial format of data - rows = individuals, cols = (e_i,r_i,s_i,censor)
###

#making dates discrete for fitting model

#e
d.cap$dateDis = (week(d.cap$date) - start.week + 1)
d.fawncap$dateDis = (week(d.fawncap$date) - start.week + 1)

#r
d.mort.ad$mORTAlertdateDis=(week(d.mort.ad$mORTAlertdate) - start.week + 1)
d.mort.ad$estMortDateDis = (week(d.mort.ad$estMortDate) - start.week + 1)

#s
d.mort.ad$collarFoundDis=(week(d.mort.ad$collarFound) - start.week + 1)



#week cut-off for harvest 
non.harvest.survival.end=week(ymd('2017-09-16'))-start.week
bow.start.week=week(ymd('2017-09-16'))-start.week+1
bow.end.week = week(ymd('2018-01-07'))+52-start.week+1
gun.start.week=week(ymd('2017-11-18'))-start.week+1
gun.end.week = week(ymd('2017-11-26'))-start.week+1
gun.holiday.start.week = week(ymd('2017-12-24'))-start.week+1
gun.holiday.end.week = week(ymd('2018-01-01'))+52-start.week+1
study.end=week(ymd('2018-01-07'))+52-start.week+1

d.mort.ad$bow=0
d.mort.ad$bow[which(d.mort.ad$weapon!="Rifle" & d.mort.ad$cause1 == "Hunter harvest")[-1]]=1

# likely gun harvest
which(d.mort.ad$weapon!="Rifle" & d.mort.ad$cause1 == "Hunter harvest")

#setting harvest.season for derived parameters
n.weeks=study.end

# non.harvest.survival.end
# bow.start.week
# bow.end.week
# gun.start.week
# gun.end.week
# gun.holiday.start.week
# gun.holiday.end.week


bow.harvest.hazard=bow.harvest.haz = c(rep(0,non.harvest.survival.end),rep(1,bow.end.week-bow.start.week+1))
bow.harvest.hazard
length(bow.harvest.haz)

gun.harvest.haz=c(rep(0,gun.start.week-1),rep(1,gun.end.week-gun.start.week+1),rep(0,gun.holiday.start.week-gun.end.week),rep(1,gun.holiday.end.week-gun.holiday.start.week))
length(gun.harvest.haz)
gun.harvest.haz


#########################################################################################################################3
###
### Formatting/combining of data
###
###
#########################################################################################################################3

rm.second.yr=which(d.cap$dateDis>12)
d.cap=d.cap[-rm.second.yr,]
n.cap=dim(d.cap)[1]

###
### N.temp[,1] = e
### N.temp[,2] =r 
### N.temp[,3] = s
### N.temp[,4] = censor
### N.temp[,5] = CWD status
### N.temp[,6] = body condition
### N.temp[,7] = gun harvest
### N.temp[,8] = sex
### N.temp[,9] = bow harvest
###

#initialize matrix of captures/recaptures/morts
N.temp = matrix(NA,nr=n.cap,nc = 10)

#initialize with all animals alive
N.temp[,4] = 1

#initialize with all CWD neg
N.temp[,5] = 0

#initialize with gun harvest mortality=0
N.temp[,7] = 0

#initialize with sex=0 males sex=1 is females
N.temp[,8] = 0

#initialize for bow hunt
N.temp[,9] = 0


#Filling N.temp with all times
for(i in 1:n.cap){
    N.temp[i,1] = d.cap$dateDis[i]#e
    N.temp[i,5] = d.cap$cwdstatus[i]
    N.temp[i,6] = d.cap$bodyconditionvalue[i]
    N.temp[i,8] = as.numeric(d.cap$sex[i])-1
    N.temp[i,10] = d.cap$lowtag[i]
    for(j in 1:n.mort){
        if(d.cap$lowtag[i]==d.mort.ad$lowtag[j]){
            N.temp[i,2] = d.mort.ad$mORTAlertdateDis[j]#r
            N.temp[i,3] = d.mort.ad$estMortDateDis[j]#s
            N.temp[i,4] = 0
            if(d.mort.ad$cause1[j]=="Hunter harvest")N.temp[i,7]=1
            N.temp[i,9] = d.mort.ad$bow[j]
        }
        if(N.temp[i,4]==0){
            if(is.na(N.temp[i,2])){
                N.temp[i,2]=d.mort.ad$estMortDateDis[j]
            }
            if(is.na(N.temp[i,3])){N.temp[i,3]=study.end}#set to the max week, because that's the censor week of last known alive
            if(N.temp[i,2]>N.temp[i,3]){N.temp[i,2]=N.temp[i,3]}
        }
    }
}
N.temp

n.temp=dim(N.temp)[1]

###
### censoring individuals that were killed at capture
###

#censor morts caused by capture
mort.killcap=which(d.mort.ad$cause1=="Capture Related")
mort.killcap
d.mort.ad[mort.killcap,]
cens.id=d.mort.ad$lowtag[mort.killcap]
N.temp
for(i in 1:n.temp){
    for(j in cens.id){
        if(N.temp[i,10]==j){
            N.temp[i,2:4]=c(study.end,study.end,1)
        }
    }
}


###
### Fill in right censor
###

for(i in 1:(n.cap)){
    if(is.na(N.temp[i,2]))N.temp[i,2]=study.end
    if(is.na(N.temp[i,3]))N.temp[i,3]=study.end
}
tail(N.temp)

#set harvest coefficient

# which(N.temp[,3]>non.harvest.survival.end)

# N.temp[which(N.temp[,3]>non.harvest.survival.end),7]=1
# which(N.temp[,3]>non.harvest.survival.end & N.temp[,4]==0)
# sum(N.temp[,7])#number of gun harvested
# N.temp[,7]
# 
# head(N.temp)
# N.temp[6,1:3]=c(1,2,2)

#making last known alive be week before mort
for(i in 1:n.temp){
    if(N.temp[i,4]==0 & N.temp[i,2]==N.temp[i,3]){
        N.temp[i,2]=N.temp[i,2]-1
    }
}

N.temp[which(N.temp[,2]<N.temp[,1]),]


###
### Post-Mortem CWD test
###

d.post.cwd$cwdresult[d.post.cwd$cwdresult=='Pending']=NA
d.post.cwd$cwdresult=as.character(d.post.cwd$cwdresult)
d.post.cwd$cwdresult[d.post.cwd$cwdresult=='']=NA
d.post.cwd$cwdresult[d.post.cwd$cwdresult=='Negative']
d.post.cwd$cwdresult[d.post.cwd$cwdresult=='Postive']=1
d.post.cwd$cwdresult=as.numeric(as.factor(d.post.cwd$cwdresult))-1

###
### removing fawns from the post mort tests
###

fawn.cwd.rm=c()
for(l in 1:n.fawncap){
    fawn.cwd.rm=c(fawn.cwd.rm,which(d.post.cwd$lowtag==d.fawncap$lowtag[l]))
}
d.post.cwd=d.post.cwd[-fawn.cwd.rm,]

n.postcwd=dim(d.post.cwd)[1]

N.temp=cbind(N.temp,rep(NA,n.temp))

for(i in 1:n.postcwd){
    for(j in 1:n.temp){
        if(d.post.cwd$lowtag[i]==N.temp[j,10]){
            N.temp[j,11]=d.post.cwd$cwdresult[i]
        }
    }
}

N.temp[which(N.temp[,5]==1),11]=1

#marking this crossbowed deer as bowed instead of rifle
#7=gun,9=bow
N.temp[N.temp[,10]==5716,]
N.temp[N.temp[,10]==5716,9]=1

###
### Format data matrix to fit into jags
###

N.data.fit = matrix(NA,nr=n.temp+n.mort-7,ncol=9)

indx = 1
for(i in 1:n.temp){
    if(N.temp[i,4]==1){
        N.data.fit[indx,] = c(N.temp[i,1:2],1,N.temp[i,5:9],NA)
        if(!is.na(N.temp[i,11])){N.data.fit[indx,9]=N.temp[i,11]}
        indx=indx+1
    }
    if(N.temp[i,4]==0){
        N.data.fit[indx,] = c(N.temp[i,1:2],1,N.temp[i,5:9],NA)
        if(!is.na(N.temp[i,11])){N.data.fit[indx,9]=N.temp[i,11]}
        indx=indx+1
        N.data.fit[indx,] = c(N.temp[i,2:3],0,N.temp[i,5:9],NA)
        if(!is.na(N.temp[i,11])){N.data.fit[indx,9]=N.temp[i,11]}
        indx=indx+1
    }
}

#indexing records
n.fit=dim(N.data.fit)[1]

#filling in known CWD + status in column 9 of N.data.fit, which is equivalent to known z's in the analysis, remaining z's must be predicted
for(i in 1:n.fit){
    if(is.na(N.data.fit[i,9]) & N.data.fit[i,4]==1)N.data.fit[i,4]=1
}

#for fast morts, remove the "living"contribution to survival lines
mort.check=which(N.data.fit[,3]==0)-1
rm.indx=c()
for(i in mort.check){
    if(N.data.fit[i,1]>=N.data.fit[i,2])rm.indx=c(rm.indx,i)
}
N.data.fit[sort(c(rm.indx,rm.indx+1)),]
N.data.fit=N.data.fit[-rm.indx,]

# indexing records
n.fit=dim(N.data.fit)[1]
n.fit

#replace those bow harvested labeled with gun harvest, as 0...
N.data.fit[which(N.data.fit[,8]==1),6]=0

#set initial values for z's
z = N.data.fit[,9]
z.indx = which(is.na(z))
n.z = length(z.indx)
# z.init=N.data.fit[,9]
# z.init[is.na(z.init)]=0
sum(N.data.fit[,4])
sum(N.data.fit[,9],na.rm=TRUE)

######################################################################################################################3
###
### Adding covariates i.e. CWD status at capture and hunter harvest and sex
###
######################################################################################################################3

###
### define Jags model 
###

sink("model.cwd.sex.R")
cat("
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
    ",fill = TRUE)
sink()

#specify initial values
inits<-list(list("beta0"=1,"beta1"=-1.5,"beta2"=1,"beta3"=.1,"p"=.6,"psi"=.04),
            list("beta0"=-.5,"beta1"=.25,"beta2"=2,"beta3"=-.3,"p"=.5,"psi"=.1),
            list("beta0"=1,"beta1"=.5,"beta2"=.5,"beta3"=1.6,"p"=.7,"psi"=.14)
            )

#identify params to monitor
parameters<-c("beta0","beta1","beta2","beta3","beta4","S0","p","psi") 


###
### x1 = cwd affect <- beta1 <- N.temp[,5].... N.data.fit[,4]
### x2 = sex effect <- beta4 <- N.temp[,8] ... N.data.fit[,7]
### x4 = for derived parameters, bow season survival
### x5 = for derived parameters, gun season survival
###

# Bundle data
jags.data <- list(records = n.fit,
                  left = N.data.fit[,1],
                  right = N.data.fit[,2],
                  censor = N.data.fit[,3],
                  x1 = N.data.fit[,4],
                  x2=N.data.fit[,8],
                  n.weeks = study.end,
                  x4 = bow.harvest.haz,
                  x5 = gun.harvest.haz,
                  z = N.data.fit[,9])

# MCMC settings
nt = 1
ni = 10000
nb = 10000
nc = 3

# call parallel version of jags using dclone
cl=makeCluster(3)
out.cwd=jags.parfit(cl, data=jags.data, params=parameters, model="model.cwd.sex.R",inits=inits,n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
stopCluster(cl)

# gelman.diag(out.cwd,multivariate = FALSE)
# plot(out.cwd)

n.weeks=max(N.data.fit,na.rm=TRUE)

fit.sum=summary(out.cwd)[[1]]
fit.quant=summary(out.cwd)[[2]]
fit.sum
fit.quant


# The Colorblind palette with grey:
cbPalette <- c( "#0072B2", "#D55E00", "#CC79A7","#999999", "#E69F00", "#56B4E9")


Weeks=rep(1:n.weeks,4)
# Sex=c(rep(0,n.weeks),rep(1,n.weeks),rep(0,n.weeks),rep(1,n.weeks))
# CWD.status=c(rep(0,2*n.weeks),rep(1,2*n.weeks))
SexCWD=c(rep(1,n.weeks),rep(2,n.weeks),rep(3,n.weeks),rep(4,n.weeks))

Survival=fit.sum[1:(4*n.weeks),1]
Lower=fit.quant[1:(4*n.weeks),1]
Upper=fit.quant[1:(4*n.weeks),5]

out=data.frame(cbind(Weeks,Survival,Lower,Upper,SexCWD))
out$SexCWD=as.factor(out$SexCWD)



pdf("Survival_v8_sex.pdf", width=8, height =6)
ggplot(data =out,aes(x = Weeks,y=Survival,group=SexCWD,color=SexCWD))+
    annotate('rect',xmax=gun.holiday.end.week, xmin=gun.holiday.start.week, ymax=Inf, ymin=-Inf, linetype=2,alpha=.5,fill="grey90")+
    annotate('rect',xmax=gun.end.week, xmin=gun.start.week-1, ymax=Inf, ymin=-Inf, linetype=2,alpha=.5,fill="grey90")+
    geom_line()+theme_bw()+
    geom_ribbon(aes(ymin=Lower,ymax=Upper,fill=SexCWD),alpha=.1,show.legend=NA,linetype=0)+
    scale_fill_manual(values=cbPalette,name="Sex(CWD Status)",labels=c("Female(-)","Male(-)","Female(+)","Male(+)"))+scale_colour_manual(values=cbPalette,name="Sex(CWD Status)",labels=c("Female(-)","Male(-)","Female(+)","Male(+)"))+
    ggtitle("Survival")+xlab("Time(Weeks)")+ylab("Survival Probability")+
    geom_vline(aes(xintercept=bow.start.week-1),linetype=3,color="grey50")+
    geom_vline(aes(xintercept=gun.start.week-1),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.end.week),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.holiday.start.week),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.holiday.end.week),linetype=2,color="grey50")+
    geom_text(x=36.5,y=1.01,label="Bow",color="grey50",size=3.3)+
    geom_text(x=45.5,y=1.01,label="Gun",color="grey50",size=3.3)
dev.off()


 ###
### Hazard analysis coefficients
### 
fit.sum[(4*n.weeks+1):(4*n.weeks+5),]

output=cbind(fit.sum[(4*n.weeks+1):(4*n.weeks+5),1:2],fit.quant[(4*n.weeks+1):(4*n.weeks+5),c(1,5)])
row.names(output)=c("CWD(-), No Harvest","CWD(+), No Harvest","Bow Harvest","Gun Harvest","Sex")
output=round(output,2)
output
write.csv(output,file="Hazard_table_v8.csv",row.names = TRUE)


print.xtable(xtable(output,caption = 'Coefficients of the hazard rate analysis where the baseline hazard shows the hazard for chronic wasting disease positive adult females not exposed to harvest, with additive hazard effects in the following rows. The first column denotes the mean of the posterior distribution with the standard deviation (SD), and equal-tailed Bayesian credible intervals in the right two columns.', 
                    align = c('c','c','c','c','c'),label="tab:hazard_v8",digits=c(0,2,2,2,2)),
             sanitize.text.function = function(x) {x},
             hline.after=c(0,5),
             include.rownames=TRUE,
             rotate.colnames = FALSE,
             caption.placement = "top")


###
### Survival Summary Table
###

out.all.harvest=cbind(fit.sum[c(n.weeks,2*n.weeks,3*n.weeks,4*n.weeks),1],fit.sum[c(n.weeks,2*n.weeks,3*n.weeks,4*n.weeks),2],fit.quant[c(n.weeks,2*n.weeks,3*n.weeks,4*n.weeks),1],fit.quant[c(n.weeks,2*n.weeks,3*n.weeks,4*n.weeks),5])
out.all.harvest=round(out.all.harvest,3)
out.all.harvest
out.all.harvest=data.frame(out.all.harvest)
out.all.harvest=cbind(rep(NA,4),rep(NA,4),out.all.harvest)
out.all.harvest[,1]=c("Negative","","Positive","")
out.all.harvest[,2]=c("Female","Male","Female","Male")
names(out.all.harvest)=c("CWD Status","Sex","Mean","SD","0.025","0.975")

survival_summary=out.all.harvest
row.names(survival_summary)=NULL
survival_summary

write.csv(survival_summary,file="survival_summary_v8.csv",row.names = FALSE)


###
### Table with survival before/after bow and gun hunts
###

pre.end=non.harvest.survival.end
bow.end=gun.start.week-1
gun1.end=gun.end.week
gun2.end=gun.holiday.end.week

fit.sum[c(pre.end,n.weeks+pre.end,2*n.weeks+pre.end,3*n.weeks+pre.end),1]

pre.harvest=cbind(fit.sum[c(pre.end,n.weeks+pre.end,2*n.weeks+pre.end,3*n.weeks+pre.end),1],fit.sum[c(pre.end,n.weeks+pre.end,2*n.weeks+pre.end,3*n.weeks+pre.end),2],fit.quant[c(pre.end,n.weeks+pre.end,2*n.weeks+pre.end,3*n.weeks+pre.end),1],fit.quant[c(pre.end,n.weeks+pre.end,2*n.weeks+pre.end,3*n.weeks+pre.end),5])
pre.harvest=round(pre.harvest,3)
pre.harvest
pre.harvest=data.frame(pre.harvest)
pre.harvest=cbind(rep(NA,4),rep(NA,4),pre.harvest)
pre.harvest[,1]=c("Negative","","Positive","")
pre.harvest[,2]=c("Female","Male","Female","Male")
names(pre.harvest)=c("CWD Status","Sex","Mean","SD","0.025","0.975")

post.bow.harvest=cbind(fit.sum[c(bow.end,n.weeks+bow.end,2*n.weeks+bow.end,3*n.weeks+bow.end),1],fit.sum[c(bow.end,n.weeks+bow.end,2*n.weeks+bow.end,3*n.weeks+bow.end),2],fit.quant[c(bow.end,n.weeks+bow.end,2*n.weeks+bow.end,3*n.weeks+bow.end),1],fit.quant[c(bow.end,n.weeks+bow.end,2*n.weeks+bow.end,3*n.weeks+bow.end),5])
post.bow.harvest=round(post.bow.harvest,3)
post.bow.harvest
post.bow.harvest=data.frame(post.bow.harvest)
post.bow.harvest=cbind(rep(NA,4),rep(NA,4),post.bow.harvest)
post.bow.harvest[,1]=c("Negative","","Positive","")
post.bow.harvest[,2]=c("Female","Male","Female","Male")
names(post.bow.harvest)=c("CWD Status","Sex","Mean","SD","0.025","0.975")

post.gun1.harvest=cbind(fit.sum[c(gun1.end,n.weeks+gun1.end,2*n.weeks+gun1.end,3*n.weeks+gun1.end),1],fit.sum[c(gun1.end,n.weeks+gun1.end,2*n.weeks+gun1.end,3*n.weeks+gun1.end),2],fit.quant[c(gun1.end,n.weeks+gun1.end,2*n.weeks+gun1.end,3*n.weeks+gun1.end),1],fit.quant[c(gun1.end,n.weeks+gun1.end,2*n.weeks+gun1.end,3*n.weeks+gun1.end),5])
post.gun1.harvest=round(post.gun1.harvest,3)
post.gun1.harvest
post.gun1.harvest=data.frame(post.gun1.harvest)
post.gun1.harvest=cbind(rep(NA,4),rep(NA,4),post.gun1.harvest)
post.gun1.harvest[,1]=c("Negative","","Positive","")
post.gun1.harvest[,2]=c("Female","Male","Female","Male")
names(post.gun1.harvest)=c("CWD Status","Sex","Mean","SD","0.025","0.975")

# post.bow2.harvest=cbind(fit.sum[c(bow2.end,n.weeks+bow2.end,2*n.weeks+bow2.end,3*n.weeks+bow2.end),1],fit.sum[c(bow2.end,n.weeks+bow2.end,2*n.weeks+bow2.end,3*n.weeks+bow2.end),2],fit.quant[c(bow2.end,n.weeks+bow2.end,2*n.weeks+bow2.end,3*n.weeks+bow2.end),1],fit.quant[c(bow2.end,n.weeks+bow2.end,2*n.weeks+bow2.end,3*n.weeks+bow2.end),5])
# post.bow2.harvest=round(post.bow2.harvest,3)
# post.bow2.harvest
# post.bow2.harvest=data.frame(post.bow2.harvest)
# post.bow2.harvest=cbind(rep(NA,4),rep(NA,4),post.bow2.harvest)
# post.bow2.harvest[,1]=c("Negative","","Positive","")
# post.bow2.harvest[,2]=c("Female","Male","Female","Male")
# names(post.bow2.harvest)=c("CWD Status","Sex","Mean","SD","0.025","0.975")

post.gun2.harvest=cbind(fit.sum[c(n.weeks,2*n.weeks,3*n.weeks,4*n.weeks),1],fit.sum[c(n.weeks,2*n.weeks,3*n.weeks,4*n.weeks),2],fit.quant[c(n.weeks,2*n.weeks,3*n.weeks,4*n.weeks),1],fit.quant[c(n.weeks,2*n.weeks,3*n.weeks,4*n.weeks),5])
post.gun2.harvest=round(post.gun2.harvest,3)
post.gun2.harvest
post.gun2.harvest=data.frame(post.gun2.harvest)
post.gun2.harvest=cbind(rep(NA,4),rep(NA,4),post.gun2.harvest)
post.gun2.harvest[,1]=c("Negative","","Positive","")
post.gun2.harvest[,2]=c("Female","Male","Female","Male")
names(post.gun2.harvest)=c("CWD Status","Sex","Mean","SD","0.025","0.975")

Period=c("Pre-Harvest","","","","Post Bow/Archery ","","","","Post Gun 1","","","","Post Gun 2","","","")

survival_all_sum=data.frame(cbind(Period,rbind(pre.harvest,post.bow.harvest,post.gun1.harvest,post.gun2.harvest)))
names(survival_all_sum)=c("Period","CWD Status","Sex","Mean","SD","0.025","0.975")

write.csv(survival_all_sum,"survival_all_sum_v8.csv",row.names=FALSE)


print.xtable(xtable(survival_all_sum,caption = 'The summary statistics ', 
                    align = c('c','c','c','c','c','c','c','c'),label="tab:survival_v8",digits=c(0,0,2,2,2,2,2,2)),
             sanitize.text.function = function(x) {x},
             hline.after=c(-1,-1,0,4,8,12,16,16),
             add.to.row = list(pos=list(c(2,6,10,14)), command="\\hdashline \n"),
             include.rownames=FALSE,
             rotate.colnames = FALSE,
             caption.placement = "top")  


#saveworking directory
save.image("survival_v8.Rdata")
