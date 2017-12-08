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

setwd('F:/171206_Survival')

library(rjags)
library(dclone)
library(doParallel)
library(foreach)
library(coda)
library(RODBC)


###
### Load Data and fix column names
###

database = odbcConnectAccess('C:/Users/aketz/Documents/Data/SWDPPdeerDB.MDB')

d.cap=sqlFetch(database, "Captures_adult_form_2017")
d.mort=sqlFetch(database, "Mortalities")
d.cap=sqlFetch(database, "")


names(d.cap)=gsub(" ","",names(d.cap)) 
names(d.cap)=gsub("[^A-Za-z0-9,;._-]","",names(d.cap))
names(d.cap)=gsub("-","",names(d.cap))
names(d.cap)=tolower(names(d.cap))
names(d.cap)

names(d.mort)=gsub(" ","",names(d.mort)) 
names(d.mort)=gsub("[^A-Za-z0-9,;._-]","",names(d.mort))
names(d.mort)=gsub("-","",names(d.mort))
names(d.mort)=tolower(names(d.mort))


names(d.cens)=gsub(" ","",names(d.cens)) 
names(d.cens)=gsub("[^A-Za-z0-9,;._-]","",names(d.cens))
names(d.cens)=gsub("-","",names(d.cens))
names(d.cens)=tolower(names(d.cens))





###
### EDA of capture data
###


d.cap$mortality

sum(d.cap$censor)
d.cap[which(as.logical(d.cap$censor)),1:5]
table(d.cap$pregnant)




