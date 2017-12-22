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
library(gridExtra)


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
names(d.cap)=c("InDB","DateEntered","CaptureStatus","lowtag","Date","Datalogger","Fieldcrew", "StudyArea","Location","Lat", "Long", "TrapType", "TrapNumber","Drughandler", "Timeofinitialstress", "Timeofinitialinjection", "DoseBAM","Bottle", "Timeofinduction","Collar", "Frequency" , "CollarID", "Magnetremoved" ,"Earpunch", "EartagLeft","EartagRight", "Ageclass", "Confidenceinage", "Sex","Hindlegcm", "Girthcm","Fecalsample", "RAMALT","RAMALTlocation", "Samplequality","RAMALTsamplerlastname", "Toothextraction","Bodycondition", "Bodyconditionvalue","Ultrasound", "UltrasoundPicID","VIT", "VITfrequency","Pregnant", "Xoffetuses","Weight", "ReversalAtipamezoleDose","Atipamezolebottle", "ReversalTime","TimeofRecovery", "Observations","EnteredCheck", "Mortality","Censor")
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
        if(d.cap$lowtag[i]==d.cwd$lowtag[j]){
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

max(d.mort.ad$collarFoundDis,na.rm=TRUE)


#week cut-off for harvest 
non.harvest.survival.end=week(ymd('2017-09-16'))-start.week
non.harvest.survival.end
harvest.start.week=week(ymd('2017-09-16'))-start.week+1
harvest.start.week

gun.start.week=week(ymd('2017-11-18'))-start.week+1
gun.start.week

study.end=max(d.mort.ad$estMortDateDis,na.rm=TRUE)
study.end

d.mort.ad$bow=0
d.mort.ad$bow[which(d.mort.ad$cause1=="Hunter harvest" & d.mort.ad$estMortDateDis<gun.start.week)]=1
         

#########################################################################################################################3
###
### Formatting/combining of data
###
###
#########################################################################################################################3



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
N.temp



#set harvest coefficient

# which(N.temp[,3]>non.harvest.survival.end)

# N.temp[which(N.temp[,3]>non.harvest.survival.end),7]=1
# which(N.temp[,3]>non.harvest.survival.end & N.temp[,4]==0)
sum(N.temp[,7])#number of gun harvested
N.temp[,7]

head(N.temp)
# N.temp[6,1:3]=c(1,2,2)

#making last known alive be week before mort
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

N.data.fit = matrix(NA,nr=n.temp+n.mort-6,ncol=8)

indx = 1
for(i in 1:n.temp){
    if(N.temp[i,4]==1){
        N.data.fit[indx,] = c(N.temp[i,1:2],1,N.temp[i,5:9])
        indx=indx+1
    }
    if(N.temp[i,4]==0){
        N.data.fit[indx,] = c(N.temp[i,1:2],1,N.temp[i,5:9])
        indx=indx+1
        N.data.fit[indx,] = c(N.temp[i,2:3],0,N.temp[i,5:9])
        indx=indx+1
    }
}

#indexing records
n.fit=dim(N.data.fit)[1]



#for fast morts, remove the "living"contribution to survival lines
mort.check=which(N.data.fit[,3]==0)-1
rm.indx=c()
for(i in mort.check){
    if(N.data.fit[i,1]>N.data.fit[i,2])rm.indx=c(rm.indx,i)
}
rm.indx

N.data.fit[sort(c(rm.indx,rm.indx+1)),]
N.data.fit=N.data.fit[-rm.indx,]


# indexing records
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

N.data.fit[1:10,]
ch.indxx=which(N.data.fit[,2]==N.data.fit[,1])
ch.indxx
N.data.fit[ch.indxx,]
N.data.fit=N.data.fit[-ch.indxx,]
n.fit=dim(N.data.fit)[1]



#setting harvest.season for derived parameters
study.end
n.weeks=study.end
gun.harvest.haz=c(rep(0,gun.start.week-1),rep(1,study.end-gun.start.week+1))
length(gun.harvest.haz)

bow.harvest.haz = c(rep(0,non.harvest.survival.end),rep(1,gun.start.week-harvest.start.week),rep(0,study.end-gun.start.week+1))
bow.harvest.haz
# length(bow.harvest.haz)
sum(bow.harvest.haz)



#replace those bow harvested labeled with gun harvest, as 0...
N.data.fit[which(N.data.fit[,8]==1),6]=0

#remove gun and bow harvest from the living line prior to death line.
# N.data.fit[which(N.data.fit[,3]==0 & N.data.fit[,6]==1)-1,6]=0
# 
# N.data.fit[which(N.data.fit[,3]==0 & N.data.fit[,8]==1)-1,8]=0

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
    
    # Priors
    beta0 ~ dnorm(0, 1.0E-6) # intercept log hazard rate
    beta1 ~ dnorm(0, 1.0E-6) # cwd effect log hazard rate
    beta2 ~ dnorm(0, 1.0E-6) # bow harvest effect log hazard rate
    beta3 ~ dnorm(0, 1.0E-6) # gun harvest effect log hazard rate

    # Likelihood for the total hazard
    for (j in 1:records) {
        for (k in left[j]:(right[j]-1)) {
            UCH[j,k] <- exp(beta0 + beta1*x1[j]+beta2*x4[k]+beta3*x5[k])
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
### x4 = for derived parameters, bow season survival
### x5 = for derived parameters, gun season survival
###



# Bundle data
jags.data <- list(records=n.fit,
                  left=N.data.fit[,1],
                  right=N.data.fit[,2],
                  censor=N.data.fit[,3],
                  x1=N.data.fit[,4],
                  n.weeks=study.end,
                  x4=bow.harvest.haz,
                  x5=gun.harvest.haz)


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


Weeks=rep(1:n.weeks,2)
CWD.status=(c(rep(0,n.weeks),rep(1,n.weeks)))

Survival=fit.sum[1:(2*n.weeks),1]
Lower=fit.quant[1:(2*n.weeks),1]
Upper=fit.quant[1:(2*n.weeks),5]

out=data.frame(cbind(Weeks,CWD.status,Survival,Lower,Upper))
out$CWD.status=as.factor(CWD.status)


pdf("Survival_v5.pdf")
 ggplot(data = out, aes(x = Weeks,y=Survival,group=CWD.status,color=CWD.status))+geom_line()+theme_bw()+
    geom_ribbon(aes(ymin=Lower,ymax=Upper,fill=CWD.status),alpha=.1,show.legend=NA,linetype=0)+
    scale_fill_manual(values=cbPalette,name="CWD Status",labels=c("Negative","Positive"))+scale_colour_manual(values=cbPalette,name="CWD Status",labels=c("Negative","Positive"))+
    ggtitle("Survival")+xlab("Time(Weeks)")+ylab("Survival Probability")+
    geom_vline(aes(xintercept=harvest.start.week-1),linetype=2,color="grey50")+
    geom_vline(aes(xintercept=gun.start.week-1),linetype=2,color="grey50")+
    geom_text(x=37,y=1,label="Bow",color="grey50")+
    geom_text(x=46.5,y=1,label="Gun",color="grey50")
dev.off()

#saveworking directory
save.image("survival_v5.Rdata")


###
### Hazard analysis coefficients
### 

output=cbind(fit.sum[(2*n.weeks+1):(2*n.weeks+4),1:2],fit.quant[(2*n.weeks+1):(2*n.weeks+4),c(1,5)])
row.names(output)=c("CWD(-), No harvest","CWD(+), No harvest","Bow Harvest","Gun Harvest")
output=round(output,2)
output
write.csv(output,file="Hazard_table_v5.csv",row.names = TRUE)

######################################################################################################################3
###
### Summary Output
###
######################################################################################################################3


###
### N.temp[,1] = e
### N.temp[,2] = r 
### N.temp[,3] = s
### N.temp[,4] = censor
### N.temp[,5] = CWD status
### N.temp[,6] = body condition
### N.temp[,7] = gun harvest
### N.temp[,8] = sex
### N.temp[,9] = bow harvest
###



###
### Frequency Tables, EDA
###


eda.out=data.frame(N.temp)
names(eda.out)=c("Enter","AliveWeek","Exit","Censor","CWDstatus","BodyCondition","Gun","Sex","Bow")

pos.eda.out=eda.out[eda.out$CWDstatus==1,]
neg.eda.out=eda.out[eda.out$CWDstatus==0,]

apply(neg.eda.out,2,sum)[c(4,7,8,9)]
apply(pos.eda.out,2,sum)[c(4,7,8,9)]


freq.table=data.frame(matrix(NA,nr=2,ncol=5))
freq.table[,1]=c("Negative","Positive")
names(freq.table)=c("CWD.status","Alive","NumberFemale","GunHarvest","BowHarvest")
freq.table[1,2:5]=apply(neg.eda.out,2,sum)[c(4,8,7,9)]
freq.table[2,2:5]=apply(pos.eda.out,2,sum)[c(4,8,7,9)]



myfreq=data.frame(t(as.matrix(freq.table[,2:5])))
names(myfreq)=c("Negative","Positive")
pdf("CWD_freq_table.pdf")
barplot(as.matrix(myfreq), main="CWD Status Frequency Table", ylab="Total", beside=TRUE, col=cbPalette[3:6])
legend(6, 80, c("Number Collars Alive","Number Collared Female","Gun Harvested","Bow Harvested"), fill=cbPalette[3:6])
dev.off()

###
### No "Alive" Category
###

freq.table=data.frame(matrix(NA,nr=2,ncol=4))
freq.table[,1]=c("Negative","Positive")
names(freq.table)=c("CWD.status","NumberFemale","GunHarvest","BowHarvest")
freq.table[1,2:4]=apply(neg.eda.out,2,sum)[c(8,7,9)]
freq.table[2,2:4]=apply(pos.eda.out,2,sum)[c(8,7,9)]

myfreq=data.frame(t(as.matrix(freq.table[,2:4])))
names(myfreq)=c("Negative","Positive")
pdf("CWD_freq_table_v2.pdf")
barplot(as.matrix(myfreq), main="CWD Status Frequency Table", ylab="Total", beside=TRUE, col=cbPalette[4:6])
legend(5, 50, c("Number Collared Female","Gun Harvested","Bow Harvested"), fill=cbPalette[4:6])
dev.off()

###
### Males vs females, breakdown
###

male.eda.out=eda.out[eda.out$Sex==1,]
female.eda.out=eda.out[eda.out$Sex==0,]


freq.table=data.frame(matrix(NA,nr=2,ncol=5))
freq.table[,1]=c("Male","Female")
names(freq.table)=c("Sex","Alive","CWD.status","GunHarvest","BowHarvest")
freq.table[1,2:5]=apply(male.eda.out,2,sum)[c(4,5,7,9)]
freq.table[2,2:5]=apply(female.eda.out,2,sum)[c(4,5,7,9)]


myfreq=data.frame(t(as.matrix(freq.table[,2:5])))
myfreq=myfreq[,2:1]
names(myfreq)=c("Females","Males")
myfreq
pdf("Sex_freq_table.pdf")
barplot(as.matrix(myfreq), main="Sex Frequency Table", ylab="Total", beside=TRUE, col=cbPalette[3:6])
legend(7, 62, c("Alive","CWD(+) at Capture", "Gun Harvested","Bow Harvested"), fill=cbPalette[3:6])
dev.off()











###################################################################################################
###
### Harvest summary plots
###
###################################################################################################


N.harvest=N.temp[which(N.temp[,7]==1),]
N.harvest[which( N.harvest[,9]==1),7]=0

total.harvest=dim(N.harvest)[1]

male.harvest=N.harvest[N.harvest[,8]==1,]
female.harvest=N.harvest[N.harvest[,8]==0,]

dim(N.harvest)[1]-sum(N.harvest[,8])

male.harvest[,5]

harvest.tab=data.frame(matrix(NA,nr=2,ncol=5))

harvest.tab[,1]=c("Male","Female")
names(harvest.tab)=c("Sex","TotalHarvest","CWD.status","GunHarvest","BowHarvest")
harvest.tab[1,2:5]=apply(male.harvest,2,sum)[c(8,5,7,9)]
harvest.tab[2,2:5]=apply(female.harvest,2,sum)[c(8,5,7,9)]
harvest.tab[2,2]=dim(female.harvest)[1]

harvest.tab


myfreq=data.frame(t(as.matrix(harvest.tab[,2:5])))
myfreq=myfreq[,2:1]
names(myfreq)=c("Females","Males")
myfreq
pdf("Harvest_freq_table.pdf")
barplot(as.matrix(myfreq), main="Harvested Collars Frequency Table", ylab="Total", beside=TRUE, col=cbPalette[3:6])
legend(1.5, 10, c("Total Harvest","CWD(+) at Capture", "Gun Harvested","Bow Harvested"), fill=cbPalette[3:6])
dev.off()





###################################################################################################
###
### the "annual" survival estimates (i.e., the mean/median and 95%CI for surviving to the end of the analysis period) broken down by harvest types
###
###################################################################################################

###
### Table with survival before/after bow and gun hunts
###
pre.end=35
bow.end=44
gun.end=n.weeks

fit.sum[c(pre.end,n.weeks+pre.end),1]

pre.harvest=cbind(fit.sum[c(pre.end,n.weeks+pre.end),1],fit.sum[c(pre.end,n.weeks+pre.end),2],fit.quant[c(pre.end,n.weeks+pre.end),1],fit.quant[c(pre.end,n.weeks+pre.end),5])
pre.harvest=round(pre.harvest,3)
pre.harvest
pre.harvest=data.frame(pre.harvest)
pre.harvest=cbind(rep(NA,2),pre.harvest)
pre.harvest[,1]=c("Negative","Positive")
names(pre.harvest)=c("CWD Status","Mean","SD","0.025","0.975")



post.bow.harvest=cbind(fit.sum[c(bow.end,n.weeks+bow.end),1],fit.sum[c(bow.end,n.weeks+bow.end),2],fit.quant[c(bow.end,n.weeks+bow.end),1],fit.quant[c(bow.end,n.weeks+bow.end),5])
post.bow.harvest=round(post.bow.harvest,3)
post.bow.harvest
post.bow.harvest=data.frame(post.bow.harvest)
post.bow.harvest=cbind(rep(NA,2),post.bow.harvest)
post.bow.harvest[,1]=c("Negative","Positive")
names(post.bow.harvest)=c("CWD Status","Mean","SD","0.025","0.975")


post.gun.harvest=cbind(fit.sum[c(n.weeks,2*n.weeks),1],fit.sum[c(n.weeks,2*n.weeks),2],fit.quant[c(n.weeks,2*n.weeks),1],fit.quant[c(n.weeks,2*n.weeks),5])
post.gun.harvest=round(post.gun.harvest,3)
post.gun.harvest
post.gun.harvest=data.frame(post.gun.harvest)
post.gun.harvest=cbind(rep(NA,2),post.gun.harvest)
post.gun.harvest[,1]=c("Negative","Positive")
names(post.gun.harvest)=c("CWD Status","Mean","SD","0.025","0.975")

Period=c("Pre-Harvest","","Bow/Archery","","Gun","")

survival_all_sum=data.frame(cbind(Period,rbind(pre.harvest,post.bow.harvest,post.gun.harvest)))
survival_all_sum
names(survival_all_sum)=c("Period","CWD Status","Mean","SD","0.025","0.975")

write.csv(survival_all_sum,"survival_all_sum_v5.csv",row.names=FALSE)
