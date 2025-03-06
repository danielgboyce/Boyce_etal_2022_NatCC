################################################################################################

##R CODE FOR: BOYCE, D. ET AL. 2021. A CLIMATE RISK INDEX FOR MARINE LIFE
##CODE AND ANALYSES DEVELOPED BY DANIEL BOYCE (DBOYCE@DAL.CA); PLEASE REFERENCE USE ACCORDINGLY
##IMPLEMENTED IN R VERSION 4.04 AND RUN ON COMPUTE CANADA COMPUTING CLUSTER
##LAST MODIFIED: AUGUST 2021

################################################################################################

##LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

cdir<-function(){
  codedir<-ifelse(grepl('sailfish',getwd(), fixed = TRUE),
                  'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                  'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
  source(paste(codedir,'GPARAMs_CC_Vuln.r',sep=''))
}
##cdir()

##FUNCTION TO ESIMATE TOE
##############################################################
##############################################################
findToE<-function(focmaxT,
                      futTs,
                      yearsFut,
                      runLength=5){

  ##IF SPECIES DOES NOT EMERGE FROM ITS NICHE BY END OF PROJECTION, SET TO THIS
  ToEFirst<-max(yearsFut)+1

  emYears<-which(futTs>focmaxT)

  if(length(emYears)>0){

    ##BINARY VECTOR SHOWING EMERGENT (1) OR NON-EMERGENT (0) YEARS
    x<-rep(0,length(yearsFut))
    x[emYears]<-1

    ##YEAR WHEN SST EMERGES FROM NICHE FOR AT LEAST NRUN YEARS
    ##WHERE VALUE CHANGES
    diffs <- x[-1L] != x[-length(x)]
    ##EXTRACT INDICES, THEN DIFFERENCE IN INDICES
    idx <- c(which(diffs), length(x))
    ##RUN LENGTH
    runs<-diff(c(0, idx))
    ##DETERMINE IF EACH RUN IS 0S OR 1S
    x2<-x[-which(c(NA,diff(x))==0)]
    ##START OF EACH RUN
    yearsRunStart<-yearsFut[-which(c(NA,diff(x))==0)]
    ##STORE
    runsTab<-data.frame(year=yearsRunStart,val=x2,run=runs)
    ##REMOVE 0S
    runsTab<-runsTab[runsTab$val==1,]
    runStart<-which(runsTab$run>=runLength)
    if(length(runStart)>0){
      ##FIRST YEAR WHERE RUN IS >= RUNLENGTH
      ToEFirst<-min(runsTab$year[runStart])
    }

  }
  return(ToEFirst)
}



##TEMPERATURE TOLERANCES
hspen<-fread(paste(datadir,'Aquamaps/2020/raw_data/hspen.csv',sep=''),header=TRUE)
names(hspen)<-tolower(names(hspen))
hspen<-subset(hspen,speciesid %in% spinclude$speciesid)
hspen<-subset(hspen,select=c('speciesid','tempmax'))

##GET CELLS IN RANGE OF EACH SPECIES
load(paste(datadirout,'hcaf_species_native_1deg.RData',sep=''))
nr<-subset(nr, speciesid %in% spinclude$speciesid)
nr<-left_join(nr,hspen,by=c('speciesid'))

##SET GLOBAL DOMAIN
domain<-unique(subset(nr,select=c('cell','lon','lat')))

nr<-subset(nr,select=c('speciesid','cell','tempmax'))
f<-function(d){
    return(list('speciesid'=d$speciesid[1],'cell'=d$cell,'tmax'=d$tempmax[1]))
}
rangeCells<-dlply(nr,.(speciesid),.fun=f,.progress='text')

rm(nr)
gc()
##rangeCells[[4]]$speciesid
##rangeCells[[4]]$tmax
##names(rangeCells[[4]])




###################################################

##1 FUNCTION AND 2 LOOPS TO CALCULATE TOES FOR EACH SPECIES ACROSS ITS DOMAIN FOR EACH CLIMATE MODEL USING PARALLEL PROCESSING
##PARALLEL PROCESSING SPEEDS UP TO ~3HRS PER CLIMATE PROJ

###################################################
##rangeCells2<-rangeCells[1:500]
##SET UP COMPUTE CLUSTER; MINE USES 7 OF MAX 8 CORES
library(doSNOW)
library(foreach)
library(doParallel)
UseCores<-32
print(UseCores)
cl <- makeCluster(UseCores, type = "SOCK")
registerDoSNOW(cl)

##DIRECTORY WHERE CLIMATE PROJECTIONS ARE
projdir<-paste(datadir,'/analysis_outputs_local/one_deg/rfiles/cmodels/SSP585/annual_max/',sep='')

##LOOP THROUGH CLIMATE PROJECTIONS
n<-1
fls<-(paste0(projdir,list.files(projdir)))
nms<-list.files(projdir)
setwd(projdir)
##fls<-fls[1]
l<-list()
    for(x in seq(n,length(fls),1)){
    print(fls[x])
    ##LOAD SST PROJECTION FOR MODEL 'X'; NAMED 'DB'
    (load(fls[x]))
    db<-subset(db,select=c(paste('X',c(seq(2015,2100,1)),sep='')))

    ##FUNSPEC OPERATES ON EACH SPECIES IN RANGECELLS (Y); ON EACH LOOP, THIS IS PARALLELLIZED
    ##y<-rangeCells2[[1]]
    funspec<-function(y){
        ##PROJECTED SSTS FOR ALL CELLS ACROSS SPECIES RANGE
        sstts<-subset(db,rownames(db) %in% y$cell)

        ##L IS HOLDER FOR OUTPUT WITHIN LOOP
        l<-list()
        ##LOOP EACH ROW(CELL) OF SST PROJECTIONS FOR SPECIES Y
        for(i in seq(1:dim(sstts)[1])){
            ##SST TS
            z<-sstts[i,]
            ##OUTPUT CELL AND TOE FOR SPECIES USING FINDTOE FUNCTION
            l[[i]]<-data.frame(cell=rownames(z),
                                       ToE=findToE(focmaxT=y$tmax,futTs=z,yearsFut=2015:2100,runLength=2))
        }
        ##RETURN LIST WITH SPECIES ID, CELL ID, AND TOES
        return(list(data.frame(rbindlist(l))))
    }
    ##IN PARALLEL; NEED TO PASS REQUIRED PACKAGES AND DATA SOURCES TO CLUSTER
    print(system.time(
    toedat<-llply(rangeCells,
                  .fun=funspec,
                  .parallel=TRUE,
                  .progress='text',
                  .paropts=list(.packages=c('data.table'),
                                .export=c('db','rangeCells','findToE')))
    ))
    l[[x]]<-toedat
rm(toedat)
gc()
}
names(l)<-nms
save(l,file=paste(datadirout,'E_ToEdata_1deg_list_85_TANmax.RData',sep=''))
stopCluster(cl)
gc()

##LISTS TO DATAFRAME
f<-function(d){
    f2<-function(d2){ return(data.frame(d2))  }
    return(ldply(d,.fun=f2,.id='speciesid'))
}
toedat<-ldply(l,.fun=f,.progress='text',.id='model')
toedat$model<-gsub('\\.rda','',toedat$model)
toedat<-spread(toedat,key=model,value=ToE)
gc()

##CONVERT TO YEARS FROM PRESENT
for(i in seq(1,length(nms),1)){
    toedat[,i+2]<-ifelse(toedat[,i+2]>=2100,2100,toedat[,i+2])
    toedat[,i+2]<-ifelse(toedat[,i+2]<2020,2020,toedat[,i+2])
    toedat[,i+2]<-toedat[,i+2]-2020
}


################################################

## EFFECT OF CC ON FUTURE RANGE OF EACH SPECIES

################################################
##% OF NATIVE RANGE LOST BY 2100 FOR EACH MODEL
##a<-subset(rcdata,speciesid==unique(rcdata$speciesid)[1])
rcdata0 <- toedat[,-2]
f<-function(a){

    f2<-function(b){
        b<-b[b<200]
        em<-b[b<=79]
        return(data.frame(nrchng=length(em)/length(b)))
    }
    return(adply(a[,-1],2,.fun=f2,.expand=TRUE,.id='model'))
}
rcdat<-ddply(rcdata0,.(speciesid),.fun=f,.progress='text')
rcdat<-spread(rcdat,key='model',value='nrchng')

##AVERAGE MM NRCHNG
rcdat$nrchng<-rowMeans(rcdat[,2:dim(rcdat)[2]],na.rm=TRUE)
rcdata<-subset(rcdat,select=c('speciesid','nrchng'))
save(rcdata,file=paste(datadirout,'E_ToEdata_NR_1deg_85_TANmax.RData',sep=''))
rm(rcdata)
gc()


################################################

## EFFECT OF CC ON ECOSYSTEM INDEX

################################################
##% OF ALL SPECIES LOST IN EACH CELL
toedataeco0 <- toedat
##a<-subset(toedataeco0,cell %in% unique(toedataeco0$cell)[1])
f<-function(a){

    ##d<-a[,3]
    f2<-function(d){
        d<-d[d<=200]
        dlost<-d[d<=79]
        return(data.frame(plost=length(dlost)/length(d)))
    }
    ##dd<-adply(a[,3:dim(a)[2]],2,.fun=f2,.id='model')
    d2<-adply(a[,-c(1,2)],2,.fun=f2,.id='model')
    return(spread(d2,key='model',value='plost'))
}
toedateco<-ddply(toedataeco0,.(cell),.fun=f,.progress='text')

##AVERAGE MM NRCHNG
toedateco$plost<-rowMeans(toedateco[,2:dim(toedateco)[2]],na.rm=TRUE)
##hist(toedateco$plost,breaks=100)
toedataeco<-subset(toedateco,select=c('cell','plost'))
save(toedataeco,file=paste(datadirout,'E_ToEdataeco_1deg_85_TANmax.RData',sep=''))
rm(toedataeco)
gc()



################################################

## ToE

################################################
toedata<-toedat

##AVERAGE MM TOE
toedata$ToE<-rowMeans(toedata[3:length(nms)+2],na.rm=TRUE)
toedata<-subset(toedata,select=c('speciesid','cell','ToE'))
save(toedata,file=paste(datadirout,'E_ToEdata_1deg_85_TANmax.RData',sep=''))
rm(toedat)
gc()

print('E_ToEs_species_done 85')






##########################################################
##########################################################
##########################################################

##SAME BUT FOR SSP2.6

##########################################################

##SET UP COMPUTE CLUSTER; MINE USES 7 OF MAX 8 CORES
library(doSNOW)
library(foreach)
library(doParallel)
UseCores<-32
print(UseCores)
cl <- makeCluster(UseCores, type = "SOCK")
registerDoSNOW(cl)

##DIRECTORY WHERE CLIMATE PROJECTIONS ARE
projdir<-paste(datadir,'/analysis_outputs_local/one_deg/rfiles/cmodels/SSP126/annual_max/',sep='')

##LOOP THROUGH CLIMATE PROJECTIONS
n<-1
fls<-(paste0(projdir,list.files(projdir)))
nms<-list.files(projdir)
setwd(projdir)
##fls<-fls[1]
l<-list()
system.time(
    for(x in seq(n,length(fls),1)){
    print(fls[x])
    ##LOAD SST PROJECTION FOR MODEL 'X'; NAMED 'DB'
    (load(fls[x]))
    db<-subset(db,select=c(paste('X',c(seq(2015,2100,1)),sep='')))

    ##FUNSPEC OPERATES ON EACH SPECIES IN RANGECELLS (Y); ON EACH LOOP, THIS IS PARALLELLIZED
    funspec<-function(y){
        ##PROJECTED SSTS FOR ALL CELLS ACROSS SPECIES RANGE
        sstts<-subset(db,rownames(db) %in% y$cell)

        ##L IS HOLDER FOR OUTPUT WITHIN LOOP
        l<-list()
        ##LOOP EACH ROW(CELL) OF SST PROJECTIONS FOR SPECIES Y
        for(i in seq(1:dim(sstts)[1])){
            ##SST TS
            z<-sstts[i,]
            ##OUTPUT CELL AND TOE FOR SPECIES USING FINDTOE FUNCTION
            l[[i]]<-data.frame(cell=rownames(z),
                                       ToE=findToE(focmaxT=y$tmax,futTs=z,yearsFut=2015:2100,runLength=2))
        }
        ##RETURN LIST WITH SPECIES ID, CELL ID, AND TOES
        return(list(data.frame(rbindlist(l))))
    }
    ##IN PARALLEL; NEED TO PASS REQUIRED PACKAGES AND DATA SOURCES TO CLUSTER
    toedat<-llply(rangeCells,
                  .fun=funspec,
                  .parallel=TRUE,
                  .progress='text',
                  .paropts=list(.packages=c('data.table'),
                                .export=c('db','rangeCells','findToE')))
    l[[x]]<-toedat
rm(toedat)
gc()
    }
)
names(l)<-nms
save(l,file=paste(datadirout,'E_ToEdata_1deg_list_26_TANmax.RData',sep=''))
stopCluster(cl)
gc()


##LISTS TO DATAFRAME
f<-function(d){
    f2<-function(d2){ return(data.frame(d2))  }
    return(ldply(d,.fun=f2,.id='speciesid'))
}
toedat<-ldply(l,.fun=f,.progress='text',.id='model')
toedat$model<-gsub('\\.rda','',toedat$model)
toedat<-spread(toedat,key=model,value=ToE)
rm(l)
gc()

##CONVERT TO YEARS FROM PRESENT
for(i in seq(1,length(nms),1)){
    toedat[,i+2]<-ifelse(toedat[,i+2]>=2100,2100,toedat[,i+2])
    toedat[,i+2]<-ifelse(toedat[,i+2]<2020,2020,toedat[,i+2])
    toedat[,i+2]<-toedat[,i+2]-2020
}

################################################

## EFFECT OF CC ON FUTURE RANGE OF EACH SPECIES

################################################
##% OF NATIVE RANGE LOST BY 2100 FOR EACH MODEL
rcdata0 <- toedat[,-2]
f<-function(a){

    f2<-function(b){
        b<-b[b<200]
        em<-b[b<=79]
        return(data.frame(nrchng=length(em)/length(b)))
    }
    return(adply(a[,-1],2,.fun=f2,.expand=TRUE,.id='model'))
}
rcdat<-ddply(rcdata0,.(speciesid),.fun=f,.progress='text')
rcdat<-spread(rcdat,key='model',value='nrchng')

##AVERAGE MM NRCHNG
rcdat$nrchng<-rowMeans(rcdat[,2:dim(rcdat)[2]],na.rm=TRUE)
rcdata<-subset(rcdat,select=c('speciesid','nrchng'))
save(rcdata,file=paste(datadirout,'E_ToEdata_NR_1deg_26_TANmax.RData',sep=''))
rm(rcdata)
gc()


################################################

## EFFECT OF CC ON ECOSYSTEM INDEX

################################################
##% OF ALL SPECIES LOST IN EACH CELL
toedataeco0 <- toedat
f<-function(a){

    ##d<-a[,3]
    f2<-function(d){
        d<-d[d<=200]
        dlost<-d[d<=79]
        return(data.frame(plost=length(dlost)/length(d)))
    }
    d2<-adply(a[,-c(1,2)],2,.fun=f2,.id='model')
    return(spread(d2,key='model',value='plost'))
}
toedateco<-ddply(toedataeco0,.(cell),.fun=f,.progress='text')

##AVERAGE MM NRCHNG
toedateco$plost<-rowMeans(toedateco[,2:dim(toedateco)[2]],na.rm=TRUE)
toedataeco<-subset(toedateco,select=c('cell','plost'))
save(toedataeco,file=paste(datadirout,'E_ToEdataeco_1deg_26_TANmax.RData',sep=''))
rm(toedataeco)
gc()



################################################

## ToE

################################################
toedata<-toedat

##AVERAGE MM TOE
toedata$ToE<-rowMeans(toedata[3:length(nms)+2],na.rm=TRUE)
toedata<-subset(toedata,select=c('speciesid','cell','ToE'))
save(toedata,file=paste(datadirout,'E_ToEdata_1deg_26_TANmax.RData',sep=''))

rm(list=ls(all=TRUE))
gc()
print('E_ToEs_species_done 26')








###########################################################


##########################################################

##########################################################
efun<-function(){
a<-toedataeco
par(mfrow=c(3,4))
plot(a$plost,a$e.plost)
km<-0.1##VELOCITY AT WHICH 0.5 IS REACHED
vmax<-1##ASYMPTOTIC VALUE
a$e.plost<-vmax*(a$plost/(km+a$plost))
plot(a$plost,a$e.plost,ylim=c(0,1))

toedataeco$e.plost<-(1-exp(-toedataeco$plost/.25))

a<- toedata[sample(nrow(toedata), 20000), ]
km<-5##VELOCITY AT WHICH 0.5 IS REACHED
vmax<-1##ASYMPTOTIC VALUE
par(mfrow=c(3,4))
a$e.ToE<-1-vmax*(a$ToE/(km+a$ToE))
plot(a$ToE,a$e.ToE,ylim=c(0,1))
abline(h=0)
hist(a$ToE,breaks=150)
hist(a$e.ToE,breaks=150)

##toedata1<-data.frame(toedat)
##toedata<-toedata1[,grepl('ToE',names(toedata1))==TRUE]
##names(toedata)<-gsub('.rda.ToE','',names(toedata))
##toedata$cell<-toedata1[,grepl('rda.cell',names(toedata1))==TRUE][,1]
##toedata$speciesid<-toedata1[,grepl('rda.speciesid',names(toedata1))==TRUE][,1]
##rm(toedata1)
gc()

}
