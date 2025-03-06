## CODE FORMATS ALL SST DATA (REYNOLDS AND PATHFINDER)
## REINTERPOLATES TO 1 DEGREE
## CALCULATES MAXIMUM VALUES OVER 2000-2020
## GETS COUNTS OF VALUES IN BINS OVER 1981-2020

## LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

## Function to dynamically set the working directory based on system environment
cdir<-function(){
  codedir<-ifelse(grepl('sailfish',getwd(), fixed = TRUE),
                  'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                  'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
  source(paste(codedir,'GPARAMs_CC_Vuln.r',sep=''))
}
##cdir()  ## Uncomment if directory change is needed
rasterOptions(progress="text", timer=TRUE)  ## Set raster processing options

#############################
## REYNOLDS OISST PROCESSING
#############################

## Define global 1-degree raster grid
mct<- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crds<-expand.grid(lon=seq(0.5,359.5,1),
                  lat=seq(-89.5,89.5,1),
                  sst=NA)
r1deg<-rasterFromXYZ(crds)
crs(r1deg) <- mct  ## Assign coordinate reference system

#############################################################
## REINTERPOLATE REYNOLDS OISST DATA TO 1 DEGREE
#############################################################

fls<-list.files(paste(datadir,'reynolds_OISST/',sep=''))
yrs<-seq(1980,2020,1)  ## Define years of interest
pt<-paste(yrs,collapse = '|')
fls<-str_subset(fls,pattern=pt)  ## Filter files matching year patterns

## Function to read and interpolate SST data
download_and_interpolate <- function(d){
    print(d)
    ncfname<-paste(paste(datadir,'reynolds_OISST/',sep=''), d, sep="")
    sstr<-brick(ncfname, varname='sst')  ## Load SST data as raster brick
    gc()

    ## Perform bilinear interpolation to 1-degree resolution
    return(resample(sstr, r1deg, method='bilinear'))
    gc()
}

## Parallel processing setup
library(doSNOW)
library(foreach)
library(doParallel)
UseCores <- detectCores() - 8  ## Use available cores
UseCores <- 32  ## Manually setting to 32 cores
print(UseCores)
cl <- makeCluster(UseCores, type = "SOCK")
registerDoSNOW(cl)

## Process files in parallel
l<-llply(fls,
         .fun=download_and_interpolate,
         .parallel = TRUE,
         .paropts=list(.packages=c('raster'),.export=c('fls','r1deg','datadir')))
##stopCluster(cl)  ## Uncomment to stop cluster
sstr<-stack(l)  ## Combine processed rasters into a single stack

##############################################
## CALCULATE MAXIMUM SST VALUES FOR 2010-2020
##############################################

sstr1<-subset(sstr,names(sstr)[10350:dim(sstr)[3]])  ## Extract data from 2010-2020
system.time(dat<-calc(sstr1, function(x) max(x,na.rm=TRUE)))  ## Compute max SST per cell

## Extract coordinates and save results
xy<-data.frame(coordinates(dat))
pdatar <- data.frame(lon=xy$x, lat=xy$y, sst=values(dat))
pdatar$sst<-ifelse(pdatar$sst %in% c(Inf, -Inf),NA,pdatar$sst)
pdatar$lon<-ifelse(pdatar$lon > 180, -360 + pdatar$lon, pdatar$lon)
pdatar$cell<-gsub(' ','',paste(pdatar$lon,'_',pdatar$lat))
save(pdatar,file=paste(datadir,'analysis_outputs_local/one_deg/rfiles/SST/Reynolds_OISST_daily_1deg_MaxSST_2010_2020.RData',sep=''))
print('Reynolds max SST done')

##############################################
## COUNT VALUES IN EACH 1-DEGREE SST BIN
##############################################

tbreaks<-seq(-4,45,0.25)  ## Define bin breaks for histogram
bin_count_function <- function(d){ return(hist(d,breaks=tbreaks,plot=FALSE)$counts) }

r1<-calc(sstr, bin_count_function)  ## Compute histogram counts per bin
xy<-data.frame(coordinates(r1))
dt<-data.frame(y=values(r1))
names(dt)<-seq(min(tbreaks)+0.125,max(tbreaks)-0.125,0.25)
dt$lon<-xy$x
dt$lat<-xy$y
fdatar<-dt
fdatar$lon<-ifelse(fdatar$lon > 180, -360 + fdatar$lon, fdatar$lon)
fdatar$cell<-gsub(' ','',paste(fdatar$lon,'_',fdatar$lat))
save(fdatar,file=paste(datadir,'analysis_outputs_local/one_deg/rfiles/SST/Reynolds_OISST_daily_1deg_1degCfreq_1981_2020.RData',sep=''))
print('Reynolds SST frequency done')

################################################################
## PATHFINDER SST ARRAY PROCESSING
################################################################

## Load preprocessed Pathfinder SST data
(load(paste(datadir,'analysis_outputs_local/one_deg/netcdf/SST/Pathfinder_SST_bilinear_1deg_3D_array.RData',sep='')))
sstp<-subset(sstp,c(1:14255))  ## Subset relevant years

## Parallel processing setup
library(parallel)
ncores<-detectCores()-12
ncores<-32
beginCluster(ncores,type='SOCK')

## Compute max SST for 2010-2020
sstp1<-subset(sstp,names(sstp)[10339:dim(sstp)[3]])
system.time(dat<-clusterR(sstp1, function(x) max(x,na.rm=TRUE)))
endCluster()

## Extract coordinates and save results
xy<-data.frame(coordinates(dat))
pdatap <- data.frame(lon=xy$y, lat=xy$x, sst=values(dat))
pdatap$sst<-ifelse(pdatap$sst %in% c(Inf, -Inf),NA,pdatap$sst)
pdatap$lon<-ifelse(pdatap$lon > 180, -360 + pdatap$lon, pdatap$lon)
pdatap$lat<-pdatap$lat*-1
pdatap$cell<-gsub(' ','',paste(pdatap$lon,'_',pdatap$lat))
save(pdatap,file=paste(datadir,'analysis_outputs_local/one_deg/rfiles/SST/Pathfinder_SST_daily_1deg_MaxSST_2010_2020.RData',sep=''))
print('Pathfinder max SST done')

## Final cleanup
print('Processing complete')
gc()


################################################
##COUNT OF VALUES IN EACH 1DEGREE SST sstr<-brick(l)
tbreaks<-seq(-4,45,0.25)
f<-function(d){ return(hist(d,breaks=tbreaks,plot=FALSE)$counts) }
r1<-calc(sstr, f)##41 mins
xy<-data.frame(coordinates(r1))
dt<-data.frame(y=values(r1))
names(dt)<-seq(min(tbreaks)+0.125,max(tbreaks)-0.125,0.25)
dt$lon<-xy$x
dt$lat<-xy$y
fdatar<-dt
fdatar$lon<-ifelse(fdatar$lon > 180, -360 + fdatar$lon, fdatar$lon)
fdatar$cell<-gsub(' ','',paste(fdatar$lon,'_',fdatar$lat))
save(fdatar,file=paste(datadir,'analysis_outputs_local/one_deg/rfiles/SST/Reynolds_OISST_daily_1deg_1degCfreq_1981_2020.RData',sep=''))
##print(summary(fdatar))
gc()
print('reynolds SST freq done')




