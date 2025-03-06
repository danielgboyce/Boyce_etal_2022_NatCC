################################################################################################

## R CODE FOR: BOYCE, D. ET AL. 2021. A CLIMATE RISK INDEX FOR MARINE LIFE
## CODE AND ANALYSES DEVELOPED BY DANIEL BOYCE (DBOYCE@DAL.CA); PLEASE REFERENCE USE ACCORDINGLY
## IMPLEMENTED IN R VERSION 4.04 AND RUN ON COMPUTE CANADA COMPUTING CLUSTER
## LAST MODIFIED: AUGUST 2021

################################################################################################

## LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')
rasterOptions(progress="text", timer=TRUE)

## FUNCTION TO SET CODE DIRECTORY BASED ON SYSTEM PATH
cdir<-function(){
    codedir<-ifelse(grepl('sailfish',getwd(), fixed = TRUE),
                    'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                    'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
    source(paste(codedir,'GPARAMs_CC_Vuln.r',sep=''))
}
## cdir()  # Uncomment if needed

###########################################################

## LOAD HISTORICAL SEA SURFACE TEMPERATURE (SST) VARIABILITY DATA

###########################################################

## LOAD REYNOLDS SST DATA
(load(paste(datadir,'analysis_outputs_local/one_deg/rfiles/SST/Reynolds_OISST_daily_1deg_1degCfreq_1981_2020.RData',sep='')))
## REFORMAT TO LONG FORMAT
rsst<-tidyr::gather(data=fdatar,sst,n,-c(lon,lat,cell))
rsst$sst<-gsub('X\.','-',rsst$sst)
rsst$sst<-gsub('X','',rsst$sst)
rsst$sst<-as.numeric(rsst$sst)
names(rsst)<-c('lon','lat','cell','sstr','nr')
gc()

## LOAD PATHFINDER SST DATA
(load(paste(datadir,'analysis_outputs_local/one_deg/rfiles/SST/Pathfinder_SST_daily_1deg_1degCfreq_1981_2020.RData',sep='')))
## REFORMAT TO LONG FORMAT
psst<-tidyr::gather(data=fdatap,sst,n,-c(lon,lat,cell))
psst$sst<-gsub('X\.','-',psst$sst)
psst$sst<-gsub('X','',psst$sst)
psst$sst<-as.numeric(psst$sst)
names(psst)<-c('lon','lat','cell','sstp','np')

## SORT DATA AND MERGE REYNOLDS AND PATHFINDER DATASETS
rsst<-rsst[order(rsst$lon,rsst$lat),]
psst<-psst[order(psst$lon,psst$lat),]
rsst<-data.frame(cell=rsst$cell,lon=rsst$lon,lat=rsst$lat,sstr=rsst$sstr,nr=rsst$nr,sstp=psst$sstp,np=psst$np)

## FILTER OUT NON-OCEANIC CELLS WHERE SST RANGE > 0
a<-subset(rsst,nr>0 | np>0)
dt<-data.frame(cell=sort(unique(a$cell)),
               rngr=tapply(a$sstr,a$cell,function(x) max(x)-min(x)),
               rngp=tapply(a$sstp,a$cell,function(x) max(x)-min(x)))
dt<-subset(dt,rngr>0 | rngp>0)
rsst<-subset(rsst,cell %in% dt$cell)
gc()
print('SST extract done')

####################################################

## LOAD AQUAMAPS NATIVE RANGES

####################################################

(load(paste(datadirout,'hcaf_species_native_1deg.RData',sep='')))
nr<-subset(nr,speciesid %in% spinclude$speciesid)

## LOAD SPECIES TEMPERATURE ENVELOPE INFORMATION
hspen<-fread(paste(datadir,'Aquamaps/2020/raw_data/hspen.csv',sep=''),header=TRUE)
names(hspen)<-tolower(names(hspen))
hspen<-subset(hspen,speciesid %in% spinclude$speciesid)
hspen<-subset(hspen,select=c('speciesid','tempmin','tempprefmin','tempprefmax','tempmax'))

#######################################################

## MERGE SPECIES RANGES WITH TEMPERATURE TOLERANCES
smdata<-left_join(nr,hspen,by=c('speciesid'))
smdata<-data.frame(smdata)
rm(list=c('hspen','dat','nr'))
gc()

########################################################

## CALCULATE AUC FOR SST DATA

########################################################

f0<-function(d){
    f1<-function(x,mint=d$tempmin[1],maxt=d$tempmax[1],minp=d$tempprefmin[1],maxp=d$tempprefmax[1]){
        ## CALCULATE SST RANGE FOR REYNOLDS
        x2<-subset(x,nr>0)
        trng.r<-ifelse(dim(x2)[1]>0, max(x2$sstr)-min(x2$sstr), 0)
        
        ## CALCULATE SST RANGE FOR PATHFINDER
        x2<-subset(x,np>0)
        trng.p<-ifelse(dim(x2)[1]>0, max(x2$sstp)-min(x2$sstp), 0)
        
        ## CALCULATE AREA UNDER CURVE (AUC) METRICS
        aucp.r <- tryCatch(AUC(x=x$sstr,y=x$nr,from=minp,to=maxp,method='spline'), error=function(e) NA)
        auct.r <- tryCatch(AUC(x=x$sstr,y=x$nr,from=mint,to=maxt,method='spline'), error=function(e) NA)
        
        aucp.p <- tryCatch(AUC(x=x$sstp,y=x$np,from=minp,to=maxp,method='spline'), error=function(e) NA)
        auct.p <- tryCatch(AUC(x=x$sstp,y=x$np,from=mint,to=maxt,method='spline'), error=function(e) NA)

        return(data.frame(auct.r, aucp.r, trng.r, auct.p, aucp.p, trng.p))
    }
    return(ddply(subset(rsst,cell %in% d$cell), .(cell), .fun=f1))
}

## PARALLEL PROCESSING
library(doSNOW)
library(foreach)
library(doParallel)
UseCores <- 32
cl <- makeCluster(UseCores, type = "SOCK")
registerDoSNOW(cl)

dout <- ddply(smdata, .(speciesid), .fun=f0, .progress='text', .parallel=TRUE,
              .paropts=list(.packages='DescTools', .export=c('rsst','smdata')))
stopCluster(cl)
print(summary(dout))

## STANDARDIZE AUC METRICS
## Normalize values between 0 and 1 for comparisons
dout$auctpct.r <- pmin(pmax(dout$auct.r / dout$auctot.r, 0), 1)
dout$aucppct.r <- pmin(pmax(dout$aucp.r / dout$auctot.r, 0), 1)
dout$trngs.r <- pmin(dout$trng.r, 30) / 30

## SAVE RESULTS
save(dout,file=paste(datadirout,'AC_SST_AUC_1deg.RData',sep=''))
gc()
print('AC_thermal_variab_done')

