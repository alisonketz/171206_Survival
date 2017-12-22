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
library(xlsx)


###
### Load Data and fix column names
###

#load Access Database into R using Hmisc package with the function mdb.get (must use this on Linux)
# database =mdb.get('~/Documents/Data/SWDPPdeerDB.MDB')

d.cap=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Captures_adult_form_2017")
d.fawncap=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Fawn Capture")
d.mort=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Mortalities")
d.cens=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Censor")
d.cwd=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "RAMALT")


names(d.fawncap)=tolower(gsub("[[:punct:]]","",names(d.fawncap)))
names(d.cap)=c("InDB","DateEntered","CaptureStatus","LowTag","Date","Datalogger","Fieldcrew", "StudyArea","Location","Lat", "Long", "TrapType", "TrapNumber","Drughandler", "Timeofinitialstress", "Timeofinitialinjection", "DoseBAM","Bottle", "Timeofinduction","Collar", "Frequency" , "CollarID", "Magnetremoved" ,"Earpunch", "EartagLeft","EartagRight", "Ageclass", "Confidenceinage", "Sex","Hindlegcm", "Girthcm","Fecalsample", "RAMALT","RAMALTlocation", "Samplequality","RAMALTsamplerlastname", "Toothextraction","Bodycondition", "Bodyconditionvalue","Ultrasound", "UltrasoundPicID","VIT", "VITfrequency","Pregnant", "Xoffetuses","Weight", "ReversalAtipamezoleDose","Atipamezolebottle", "ReversalTime","TimeofRecovery", "Observations","EnteredCheck", "Mortality","Censor")
names(d.mort)=c("InDB","Lowtag","Technicians","DataRecorder","Lefteartag","Lefteartag","Righteartag","Righteartag","Collar","Collar","Sex","Captured","MORTAlertdate","CollarFound", "EstMortDate","Comment1","PropertyType","Landtype","County","Property", "Photos","Incisor","Lymphnodes","CWDSuspect","FieldNecropsy","LabNecropsy","SecondaryFateSheet","Cause1", "Cause1wt","Cause2","Cause2wt","Cause3","Cause3wt","Comments2","HunterName","Huntingseason","Weapon","Lat","Long","CustomerID","ConfidenceinCause","NumberAntlerPoints","Township","ApproximateLocationofHarvest","NumberDeerInGroup","Howfarawaywasdeerwhenshot", "Emailaddress","CWDsampleID","Legmeasurment","Girth","Dressedweight","Hunterphone","HunterAddress","Frequency")
names(d.cens)=c("Lowtag", "Lefteartag","Righteartag","Sex","AgeAtCapture","Latitude","Longitude","County","Censorcause","CensWt","Othercause1","Othercause1wt","Othercause2","Othercause2wt","Cencomments","CaptureDate","CollarRecovered","PropertyCollarRecovered","Landtype","Frequency","CensorDate","Technicians","CollarNumber")
names(d.cwd)=tolower(gsub('[[:punct:]]',"",names(d.cwd)))


for(i in 1:length(names(d.cap))){names(d.cap)[i]=paste(tolower(substr(names(d.cap)[i], 1, 1)),substr(names(d.cap)[i], 2, nchar(names(d.cap)[i])),sep="")}
for(i in 1:length(names(d.mort))){names(d.mort)[i]=paste(tolower(substr(names(d.mort)[i], 1, 1)),substr(names(d.mort)[i], 2, nchar(names(d.mort)[i])),sep="")}
for(i in 1:length(names(d.cens))){names(d.cens)[i]=paste(tolower(substr(names(d.cens)[i], 1, 1)),substr(names(d.cens)[i], 2, nchar(names(d.cens)[i])),sep="")}


###
### Number individuals within each dataframe
###

n.cap=dim(d.cap)[1]
n.cens = dim(d.cens)[1]
n.fawncap=dim(d.fawncap)[1]
n.mort = dim(d.mort)[1]
n.cwdtest=dim(d.cwd)[1]

#change class of CWD tests to integer for d.cwd
class(d.cwd$lowtag)="integer"

#add CWD status to capture dataframe
d.cap$cwdstatus=rep(NA,n.cap)
for(i in 1:n.cap){
    for(j in 1:n.cwdtest){
        if(d.cap$lowTag[i]==d.cwd$lowtag[j]){
            d.cap$cwdstatus[i]=as.character(d.cwd$result[j])
        }
    }
    if(d.cap$cwdstatus[i]=='ISF' | is.na(d.cap$cwdstatus[i]))d.cap$cwdstatus[i]='Negative'
}
d.cap$cwdstatus=as.factor(d.cap$cwdstatus)
d.cap$cwdstatus=as.numeric(d.cap$cwdstatus)-1

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

#remove morts caused by capture

mort.rm.killcap=which(is.na(d.mort.ad$collarFound)&is.na(d.mort.ad$mORTAlertdate))
mort.rm.killcap
d.mort.ad=d.mort.ad[-mort.rm.killcap,]
n.mort = dim(d.mort.ad)[1]


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

max(d.mort.ad$collarFoundDis,na.rm=TRUE)


#week cut-off for harvest 
non.harvest.survival.end=week(ymd('2017-09-16'))-start.week
non.harvest.survival.end
harvest.start.week=week(ymd('2017-09-16'))-start.week+1
harvest.start.week
study.end= week(max(d.mort.ad$collarFound,na.rm=T))-start.week+1
study.end

#initialize matrix of captures/recaptures/morts
N.temp = matrix(NA,nr=n.cap,nc = 8)

#initialize with all animals alive
N.temp[,4] = 1

#initialize with all CWD neg
N.temp[,5] = 0

#initialize with harvest mortality=0
N.temp[,7] = 0
#initialize with harvest mortality=0
N.temp[,8] = 0

#Filling N.temp with all times
for(i in 1:n.cap){
    N.temp[i,1] = d.cap$dateDis[i]#e
    N.temp[i,5] = d.cap$cwdstatus[i]
    N.temp[i,6] = d.cap$bodyconditionvalue[i]
    N.temp[i,8] = as.numeric(d.cap$sex[i])-1
    for(j in 1:n.mort){
        if(d.cap$lowTag[i]==d.mort.ad$lowtag[j]){
            N.temp[i,2] = d.mort.ad$mORTAlertdateDis[j]#r
            N.temp[i,3] = d.mort.ad$estMortDateDis[j]#s
            N.temp[i,4] = 0
            if(d.mort.ad$cause1[j]=="Hunter harvest")N.temp[i,7]=1
        }
        if(N.temp[i,4]==0){
            if(is.na(N.temp[i,2])){
                N.temp[i,2]=d.mort.ad$estMortDateDis[j]
            }
            if(is.na(N.temp[i,3])){N.temp[i,3]=study.end}#set to the max week, because that's the following week of last known alive
            if(N.temp[i,2]>N.temp[i,3]){N.temp[i,2]=N.temp[i,3]}
        }
    }
}
N.temp


for(i in 1:(n.cap)){
    if(is.na(N.temp[i,2]))N.temp[i,2]=study.end
    if(is.na(N.temp[i,3]))N.temp[i,3]=study.end
}
N.temp

n.temp=dim(N.temp)[1]
n.temp

#set harvest coefficient

# which(N.temp[,3]>non.harvest.survival.end)

# N.temp[which(N.temp[,3]>non.harvest.survival.end),7]=1
# which(N.temp[,3]>non.harvest.survival.end & N.temp[,4]==0)
sum(N.temp[,7])

head(N.temp)
N.temp[6,1:3]=c(1,2,2)

for(i in 1:n.temp){
    if(N.temp[i,4]==0 & N.temp[i,2]==N.temp[i,3]){
        N.temp[i,2]=N.temp[i,2]-1
    }
}
head(N.temp)

N.temp[which(N.temp[,2]<N.temp[,1]),]

###
### Format data matrix to fit into jags
###

N.data.fit = matrix(NA,nr=n.temp+n.mort-1,ncol=7)

indx = 1
for(i in 1:n.temp){
    if(N.temp[i,4]==1){
        N.data.fit[indx,] = c(N.temp[i,1:2],1,N.temp[i,5:8])
        indx=indx+1
    }
    if(N.temp[i,4]==0){
        N.data.fit[indx,] = c(N.temp[i,1:2],1,N.temp[i,5:8])
        indx=indx+1
        N.data.fit[indx,] = c(N.temp[i,2:3],0,N.temp[i,5:8])
        indx=indx+1
    }
}

N.data.fit

#indexing records
n.fit=dim(N.data.fit)[1]



#remove the rows where morts are withiin 2 weeks of capture

mort.check=which(N.data.fit[,3]==0)-1
rm.indx=c()
for(i in mort.check){
    if(N.data.fit[i,1]>N.data.fit[i,2])rm.indx=c(rm.indx,i)
}
rm.indx


N.data.fit[sort(c(rm.indx,rm.indx+1)),]

N.data.fit=N.data.fit[-rm.indx,]


#indexing records
n.fit=dim(N.data.fit)[1]
n.fit


#ensure no indexes of 0(
for(i in 1:n.fit){
    if(N.data.fit[i,3]==0 & N.data.fit[i,1]==N.data.fit[i,2]) N.data.fit[i,1] = N.data.fit[i,1] -1
}

#indexing records
n.fit=dim(N.data.fit)[1]
n.fit


#fix index 6
head(N.data.fit)

min(N.data.fit[,2])

ch.indxx=which(N.data.fit[,2]==N.data.fit[,1])
ch.indxx
N.data.fit[ch.indxx,]
N.data.fit=N.data.fit[-ch.indxx,]
n.fit=dim(N.data.fit)[1]



#setting harvest.season for derived parameters
study.end
n.weeks=study.end
harvest.haz=c(rep(0,non.harvest.survival.end),rep(1,study.end-harvest.start.week+1))
length(harvest.haz)



which(d.cap$lowTag[43]==d.mort.ad$lowtag)
d.mort.ad[22,]


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
inits
#identify params to monitor
parameters<-c("llambda") 

# Bundle data
jags.data <- list(records=n.fit,left=N.data.fit[,1],right=N.data.fit[,2],censor=N.data.fit[,3])

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
### Adding covariates i.e. CWD status at capture and hunter harvest
###
######################################################################################################################3


###
### define Jags model 
###

sink("model.cwd.R")
cat("
    model{
    beta0 ~ dnorm(0, 1.0E-6) # intercept log hazard rate
    beta1 ~ dnorm(0, 1.0E-6) # cwd effect log hazard rate
    beta2 ~ dnorm(0, 1.0E-6) # bow harvest effect log hazard rate
    beta3 ~ dnorm(0, 1.0E-6) # gun harvest effect log hazard
    beta4 ~ dnorm(0, 1.0E-6) # sex effect log hazard

    # Likelihood for the total hazard
    for (j in 1:records) {
        for (k in left[j]:(right[j]-1)) {
             UCH[j,k] <- exp(beta0 + beta1*x1[j]+beta2*x4[k]+beta3*x5[k] +beta4*x2[j])
        }
        # total prob of surviving
        SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
        censor[j] ~ dbern(SLR[j])      
    }

    #Derived parameters
    for (t in 1:n.weeks){
        llambda.out[t,1]<-beta0 + beta2*x4[t] +beta3*x5[t] #Female CWD-
        llambda.out[t,2]<-beta0 + beta2*x4[t] +beta3*x5[t] + beta4 #Male CWD -
       llambda.out[t,3]<-beta0 + beta1 + beta2*x4[t] +beta3*x5[t] #Female CWD+
        llambda.out[t,4]<-beta0+beta1+beta2*x4[t]+beta3*x5[t]+beta4 #Male CWD +
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
inits<-list(list("beta0"=1,"beta1"=-1.5,"beta2"=1,"beta3"=-.35),list("beta0"=-.5,"beta1"=.25,"beta2"=2,"beta3"=.5),list("beta0"=1,"beta1"=.5,"beta2"=.5,"beta3"=-1.5))

#identify params to monitor
parameters<-c("beta0","beta1","beta2","S0") 

# Bundle data

jags.data <- list(records=n.fit,left=N.data.fit[,1],right=N.data.fit[,2],censor=N.data.fit[,3],x1=N.data.fit[,4],x2=N.data.fit[,6],x3=N.data.fit[,7],n.weeks=study.end,x4=harvest.haz)



# MCMC settings
nt = 1
ni = 1000
nb = 10000
nc = 3

# call parallel version of jags using dclone
cl=makeCluster(3)
out.cwd=jags.parfit(cl, data=jags.data, params=parameters, model="model.cwd.R",inits=inits,n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
stopCluster(cl)

# gelman.diag(out.cwd)
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

ggplot(data =out,aes(x = Weeks,y=Survival,group=SexCWD,color=SexCWD))+geom_line()+theme_bw()+
    geom_ribbon(aes(ymin=Lower,ymax=Upper,fill=SexCWD),alpha=.1,show.legend=NA,linetype=0)+
    scale_fill_manual(values=cbPalette,name="Sex(CWD Status)",labels=c("Female(-)","Male(-)","Female(+)","Male(+)"))+scale_colour_manual(values=cbPalette,name="Sex(CWD Status)",labels=c("Female(-)","Male(-)","Female(+)","Male(+)"))+
    ggtitle("Survival (Including Harvest)")+xlab("Time(Weeks)")


pdf("Survival_v4.pdf")
ggplot(data =out,aes(x = Weeks,y=Survival,group=SexCWD,color=SexCWD))+geom_line()+theme_bw()+
    geom_ribbon(aes(ymin=Lower,ymax=Upper,fill=SexCWD),alpha=.1,show.legend=NA,linetype=0)+
    scale_fill_manual(values=cbPalette,name="Sex(CWD Status)",labels=c("Male(-)","Female(-)","Male(+)","Female(+)"))+scale_colour_manual(values=cbPalette,name="Sex(CWD Status)",labels=c("Male(-)","Female(-)","Male(+)","Female(+)"))+
    ggtitle("Survival (Including Harvest)")+xlab("Time(Weeks)")
dev.off()

 





######################################################################################################################3
###
### Include harvested data, not in covariate
### Adding covariates i.e. CWD status at capture and hunter harvest
###
######################################################################################################################3

###
### define Jags model 
###

sink("model.cwd.R")
cat("
    model{
    
    # Priors
    for(j in 1:records){
    llambda[j] <- beta0 + beta1*x1[j]
    }
    
    beta0 ~ dnorm(0, 1.0E-6) # intercept log hazard rate
    beta1 ~ dnorm(0, 1.0E-6) # cwd effect log hazard rate
    
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
inits<-list(list("beta0"=1,"beta1"=-1.5),list("beta0"=-.5,"beta1"=.25),list("beta0"=1,"beta1"=.5))

#identify params to monitor
parameters<-c("beta0","beta1","S0") 

# Bundle data
jags.data <- list(records=n.fit,left=N.data.fit[,1],right=N.data.fit[,2],censor=N.data.fit[,3],x1=N.data.fit[,4],n.weeks=max(N.data.fit,na.rm=TRUE))



# MCMC settings
nt = 1
ni = 1000
nb = 10000
nc = 3

# call parallel version of jags using dclone
cl=makeCluster(3)
out.cwd=jags.parfit(cl, data=jags.data, params=parameters, model="model.cwd.R",inits=inits,n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
stopCluster(cl)


gelman.diag(out.cwd,multivariate = F)

n.weeks=max(N.data.fit,na.rm=TRUE)

fit.sum=summary(out.cwd)[[1]]
fit.quant=summary(out.cwd)[[2]]
fit.sum
fit.quant





# The Colorblind palette with grey:
cbPalette <- c( "#0072B2", "#D55E00", "#CC79A7","#999999", "#E69F00", "#56B4E9")






###
### Plot with only CWD Status
###


Weeks=rep(1:n.weeks,4)
 # Sex=c(rep(0,n.weeks),rep(1,n.weeks),rep(0,n.weeks),rep(1,n.weeks))
CWD.status=c(rep(0,n.weeks),rep(1,n.weeks))

Survival=fit.sum[1:(2*n.weeks),1]
Lower=fit.quant[1:(2*n.weeks),1]
Upper=fit.quant[1:(2*n.weeks),5]

out=data.frame(cbind(Weeks,Survival,Lower,Upper,CWD.status))
out$CWD.status=as.factor(out$CWD.status)


ggplot(data =out,aes(x = Weeks,y=Survival,group=CWD.status,color=CWD.status))+geom_line()+theme_bw()+
    geom_ribbon(aes(ymin=Lower,ymax=Upper,fill=CWD.status),alpha=.1,show.legend=NA,linetype=0)+
    scale_fill_manual(values=cbPalette,name="CWD Status",labels=c("Negative","Positive"))+scale_colour_manual(values=cbPalette,name="CWD Status",labels=c("Negative","Positive"))+
    ggtitle("Survival (Including Harvest)")+xlab("Time(Weeks)")


pdf("Survival_v4_b.pdf")
ggplot(data =out,aes(x = Weeks,y=Survival,group=CWD.status,color=CWD.status))+geom_line()+theme_bw()+
    geom_ribbon(aes(ymin=Lower,ymax=Upper,fill=CWD.status),alpha=.1,show.legend=NA,linetype=0)+
    scale_fill_manual(values=cbPalette,name="CWD Status",labels=c("Negative","Positive"))+scale_colour_manual(values=cbPalette,name="CWD Status",labels=c("Negative","Positive"))+
    ggtitle("Survival (Including Harvest)")+xlab("Time(Weeks)")

dev.off()


#saveworking directory
save.image("survival_v4.Rdata")


###
### Calculate summary statistics 
###

#raw stats

#(i.e., #alive/# collars deployed) for CWD+ and CWD-
non.cwd=N.temp[N.temp[,5]!=1,]
1-sum(non.cwd[,4]==0)/dim(non.cwd)[1]


cwd.pos=N.temp[N.temp[,5]==1,]
1-sum(cwd.pos[,4]==0)/dim(cwd.pos)[1]



fit.sum
cwd.id=d.cap$lowTag[d.cap$cwdstatus==1]

mort.cwd.sub=matrix(NA,nr=length(cwd.id),nc=dim(d.mort.ad)[2])
for(i in 1:length(cwd.id)){
    for(j in 1:dim(d.mort.ad)[1]){
        if(cwd.id[i]==d.mort.ad$lowtag[j]){
            mort.cwd.sub[i,]=d.mort.ad[j,]
        }
        else{
            mort.cwd.sub[i,]=rep(1,57)
        } 
    }
}
d.mort.ad$lowtag==

d.cap$cwdstatus
dim(d.cap)
dim(d.mort.ad)

names(d.mort.ad)
d.mort.ad.hunt=d.mort.ad[d.mort.ad$cause1=="Hunter Harvest",]

hist(d.mort.ad.hunt)
