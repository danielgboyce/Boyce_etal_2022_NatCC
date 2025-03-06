##############################################################################################

##R CODE FOR: BOYCE, D. ET AL. 2021. A CLIMATE RISK INDEX FOR MARINE LIFE
##CODE AND ANALYSES DEVELOPED BY DANIEL BOYCE (DBOYCE@DAL.CA); PLEASE REFERENCE USE ACCORDINGLY
##IMPLEMENTED IN R VERSION 4.04 AND RUN ON COMPUTE CANADA COMPUTING CLUSTER
##LAST MODIFIED: AUGUST 2021

##############################################################################################

##LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

cdir<-function(){
    codedir<-ifelse(grepl('sailfish',getwd(), fixed = TRUE),
                    'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                    'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
    source(paste(codedir,'GPARAMs_CC_Vuln.r',sep=''))
}

print('S_HII_start')

#############################################################
##HALPERN 2013 HUMAN IMPACTS INDEX
hiidat<-raster(paste(datadir,'halpern_HII/2013/cumulative_impact_2013_degrees.tif',sep=''))
mrct<-'+proj=longlat +datum=WGS84 +no_defs'
s <- raster(nrow=180, ncol=360,crs=mrct)
hiidat.bl <- resample(hiidat, s, method='bilinear')

xy<-coordinates(hiidat.bl)
hiidata<-data.frame(lon=xy[,1],
               lat=xy[,2],
               hii=values(hiidat.bl))

##SAVE AS RDATA
save(hiidata,file=paste(datadirout,'S_HII_1deg.RData',sep=''))

rm(list=ls(all=TRUE))
gc()
print('S_HII_done')

