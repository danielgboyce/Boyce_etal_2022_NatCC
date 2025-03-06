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

#######################################################

##GET AREA OCCUPIED BY ALL SPECIES COMBINED - BASICALLY ENTIRE HABITAT
(load(paste(datadirout,'hcaf_species_native_1deg.RData',sep='')))
print(dim(nr))
a<-subset(nr,probability>0,select=c('lon','lat'))
a$probability<-1
a<-unique(a)
ag<-rasterFromXYZ(a)
crs(ag) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
globarea<-tapply(area(ag), ag[], sum)
globarea<-379831425 #FROM TAPPLY
print('globarea done')

#TOTAL AREA OF SPECIES LARGER THAN TOTAL OCEAN AREA LISTED ON WIKIPEDIA, BUT IS VERY CLOSE
##361900000#FROM WIKIPEDIA
##379831425 #FROM TAPPLY
##a<-area(ag)
##a<-data.frame(rasterToPoints(a))
##(18919364/379831425)*100

###############################################
##RASTER STOCK CONTAINING ALL NATIVE RANGES
rdatab<-stack(paste(rasterdir,'hcaf_species_native_1deg_rstack_BINARY.grd',sep=''))
print(dim(rdatab))
print(names(subset(rdatab,1:10)))

##MERCATOR PROJECTION
crs(rdatab) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


#GET AREA OF EACH RASTER LAYER - MAY TAKE AWHILE
l<-matrix(NA,nrow=nlayers(rdatab),ncol=2)
for (i in 1:nlayers(rdatab)){
    r<-rdatab[[i]]
    ##r[r <= 0] <- NA
    ar<-round(sum(tapply(area(r), r[], sum)),digits=0)
    #compute area [km2] of all cells in geo_raster
    l[i,]<-c(speciesid=names(r),
                       rarea=ar)
}
hrange<-data.frame(l)
print('loop done')

names(hrange)<-c('speciesid','areakm2')
hrange$areakm2<-as.numeric(as.character(hrange$areakm2))
hrange$gareakm2<-361900000#FROM WIKIPEDIA
hrange$areapct<-(hrange$areakm2/hrange$gareakm2)*100
hrange$rarea<-hrange$areapct/100
hrange$speciesid<-gsub('\\.','\\-',hrange$speciesid)


###################################################

##LATITUDE RANGE FOR EACH SPECIES

###################################################
rdata<-stack(paste(rasterdir,'hcaf_species_native_1deg_rstack.grd',sep=''))

a<-rdata
l<-list()
for(i in seq(1,dim(a)[3],1)){
  print(i)
  d<-subset(a,i)
  n<- Which(is.na(d)==FALSE, cells=TRUE)
  if(length(n)>0){
    ext<-extent(trim(d,values=NA))
    l[[i]]<-data.frame(speciesid=names(d),
                       lat1=ext[3],
                       lat2=ext[4])
  } else NULL
}
latdat<-data.frame(rbindlist(l))
latdat$speciesid<-gsub('\\.','\\-',latdat$speciesid)
latdat$lrange<-ifelse(abs(latdat$lat1)>abs(latdat$lat2),abs(latdat$lat1),abs(latdat$lat2))
latdat$lrangetot<-abs(latdat$lat1-latdat$lat2)
##save(latdat,file=paste(datadirout,'AC_latrange_1deg.RData',sep=''))
##(load(paste(datadirout,'AC_latrange_1deg.RData',sep='')))
latdat<-subset(latdat,select=c('speciesid','lrange'))
hrange<-left_join(hrange,latdat,by=c('speciesid'))
hrange$lrange<-ifelse(is.na(hrange$lrange)==TRUE,0,hrange$lrange)


##############################################################
##SAVE AS RDATA
save(hrange,file=paste(datadirout,'AC_hrange_1deg.RData',sep=''))

rm(list=ls(all=TRUE))
gc()
print('AC_hrange_done')

