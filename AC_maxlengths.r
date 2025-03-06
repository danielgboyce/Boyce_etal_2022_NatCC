################################################################################################

## R CODE FOR: BOYCE, D. ET AL. 2021. A CLIMATE RISK INDEX FOR MARINE LIFE
## CODE AND ANALYSES DEVELOPED BY DANIEL BOYCE (DBOYCE@DAL.CA); PLEASE REFERENCE USE ACCORDINGLY
## IMPLEMENTED IN R VERSION 4.04 AND RUN ON COMPUTE CANADA COMPUTING CLUSTER
## LAST MODIFIED: AUGUST 2021

################################################################################################

## LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

## FUNCTION TO SET CODE DIRECTORY BASED ON WORKING DIRECTORY
cdir<-function(){
    codedir<-ifelse(grepl('sailfish',getwd(), fixed = TRUE),
                    'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                    'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
    source(paste(codedir,'GPARAMs_CC_Vuln.r',sep=''))
}
## cdir() ## Unused function, commented out

## LOAD MAXIMUM BODY LENGTH DATA FROM DIFFERENT SOURCES
## FishBase dataset
fbdat<-fread(paste(datadir,'traits/Reyes_MaxTL/FB_Fish_Maximum_Lengths.csv',sep=''),header=TRUE,skip=1)
names(fbdat)<-tolower(names(fbdat))
fbdat<-subset(fbdat,select=c('speciesid','maxtl','maxwd'))

## SeaLifeBase dataset for non-fish species
sldat<-fread(paste(datadir,'traits/Yap_MaxLength_Inverts/SLB_NonFish_MaxLengths_20200725.csv',sep=''),header=TRUE,skip=1)
names(sldat)<-tolower(names(sldat))
sldat<-subset(sldat,select=c('speciesid','maxtl','maxml','maxcl','maxshh'))

## COMBINE THE TWO DATASETS
ldat<-rbind.fill(fbdat,sldat)
## RETAIN ONLY RELEVANT COLUMNS (MAX TOTAL LENGTH & MAX MANTLE LENGTH)
ldat<-subset(ldat,select=c('speciesid','maxtl','maxml'))

## SAVE COMBINED DATASET AS RDATA FILE
save(ldat,file=paste(datadirout,'AC_maxlengths_1deg_unimputed.RData',sep=''))

##########################################################
## IMPUTE MISSING BODY LENGTHS
##########################################################

## LOAD ENVIRONMENTAL ENVELOPES (HABITAT RANGE)
hspen<-fread(paste(datadir,'Aquamaps/2020/raw_data/hspen.csv',sep=''),header=TRUE)
hspen<-data.frame(hspen)
names(hspen)<-tolower(names(hspen))

## CALCULATE LATITUDE AND DEPTH RANGE
hspen$latrange<-abs(hspen$nmostlat-hspen$smostlat)
hspen$depthrange<-hspen$depthmax-hspen$depthmin

## RETAIN RELEVANT VARIABLES
hspen<-subset(hspen,select=c('speciesid','depthmax','depthmin','tempmin','tempmax','salinitymin','salinitymax',
                             'primprodmin','primprodmax','iceconmax','oxymin','oxymax','landdistmax','landdistmin'))

## LOG TRANSFORM SOME VARIABLES TO IMPROVE IMPUTATION
hspen$depthmax<-log10(hspen$depthmax+1)
hspen$primprodmin<-log10(hspen$primprodmin+1)
hspen$primprodmax<-log10(hspen$primprodmax+1)
hspen$iceconmax<-log10(hspen$iceconmax+1)
hspen$landdistmax<-log10(hspen$landdistmax+1)
hspen$landdistmin<-log10(hspen$landdistmin+1)

########################################################
## LOAD GLOBAL AREA OF EACH SPECIES
########################################################
(load(paste(datadirout,'AC_hrange_1deg.RData',sep='')))
print(length(unique(hrange$speciesid)))
hrange<-subset(hrange,select=c('speciesid','areakm2'))
hrange$areakm2<-log10(hrange$areakm2+1)
print(summary(hrange))

########################################################
## LOAD GLOBAL HABITAT FRAGMENTATION DATA
########################################################
(load(paste(datadirout,'AC_hfrag_1deg.RData',sep='')))
print(length(unique(hfrag$speciesid)))
hfrag<-subset(hfrag,select=c('speciesid','hfrag'))
hfrag$hfrag<-log10(hfrag$hfrag+1)
print(summary(hfrag))

## LOAD TAXONOMIC INFORMATION
tdat<-read.csv(paste(datadir,'Aquamaps/2020/raw_data/speciesoccursum.csv',sep=''),header=TRUE)
names(tdat)<-tolower(names(tdat))
tdat<-subset(tdat,select=c('speciesid','kingdom','phylum','order','family','genus'))
print(length(unique(tdat$speciesid)))

## COMBINE ALL DATA LAYERS
## Merge environmental, range, fragmentation, length, and taxonomic data
dat<-left_join(hspen,hrange,by=c('speciesid'))
dat<-left_join(dat,hfrag,by=c('speciesid'))
dat<-left_join(dat,ldat,by=c('speciesid'))
dat<-left_join(dat,tdat,by=c('speciesid'))

## CONVERT TAXONOMIC CATEGORIES TO NUMERIC FACTORS
dat$kingdom<-as.numeric(as.factor(dat$kingdom))
dat$phylum<-as.numeric(as.factor(dat$phylum))
dat$order<-as.numeric(as.factor(dat$order))
dat$family<-as.numeric(as.factor(dat$family))
dat$genus<-as.numeric(as.factor(dat$genus))

## PRINT DATASET DIMENSIONS AND SUMMARY
print(dim(dat))
print(summary(dat))
print(length(unique(dat$speciesid)))

## IMPUTE MISSING BODY LENGTH DATA USING RANDOM FOREST METHOD
system.time(tempdata<-mice(dat[,-1],m=5,maxit=50,meth='rf',seed=500,print=FALSE))

## GET COMPLETE IMPUTED DATASET
impdat<-complete(tempdata)

## COMPUTE AVERAGE ACROSS ALL 5 MULTIPLE IMPUTATIONS
id1<-as.matrix(complete(tempdata,1))
id2<-as.matrix(complete(tempdata,2))
id3<-as.matrix(complete(tempdata,3))
id4<-as.matrix(complete(tempdata,4))
id5<-as.matrix(complete(tempdata,5))
impdatav<-(id1+id2+id3+id4+id5)/5
impdatav<-data.frame(impdatav)
impdatav$speciesid<-dat$speciesid

## RETAIN SPECIES ID AND MAX LENGTH FOR OUTPUT
ldat<-subset(impdatav,select=c('speciesid','maxtl'))
print(length(unique(hrange$speciesid)))

## SAVE FINAL IMPUTED DATASET AS RDATA FILE
save(ldat,file=paste(datadirout,'AC_maxlengths_1deg.RData',sep=''))

## CLEAR WORKSPACE AND RUN GARBAGE COLLECTION
rm(list=ls(all=TRUE))
gc()

## PRINT COMPLETION MESSAGE
print('AC_maxlengths_done')



