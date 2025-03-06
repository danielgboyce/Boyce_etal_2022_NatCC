## LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
## Source external R script that sets global parameters
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

## Function to determine code directory based on the working directory
cdir <- function(){
  codedir <- ifelse(grepl('sailfish', getwd(), fixed = TRUE),
                    'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                    'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
  
  ## Source the parameter script from the determined directory
  source(paste(codedir, 'GPARAMs_CC_Vuln.r', sep=''))
}
## Uncomment the following line to use the function
## cdir()

##############################################################
## COMBINE STANDARDIZED VULNERABILITY LAYERS
##############################################################

## SENSITIVITY LAYERS; LARGER VALUES DENOTE HIGHER SENSITIVITY
## ORDER OF MERGE IMPORTANT

## LOAD SPECIES NATIVE RANGES TO MERGE
load(paste(datadirout, 'hcaf_species_native_1deg.RData', sep=''))
nr <- subset(nr, speciesid %in% spinclude$speciesid, select = c('speciesid', 'cell', 'lon', 'lat', 'probability'))
names(nr) <- c('speciesid.am', 'cell', 'lon', 'lat', 'prob')

## THERMAL SAFETY MARGINS (SPECIES & CELL); SDATA + TSM = SDATA
load(paste(datadirout, 'S_TSM_1deg.RData', sep=''))
smdata <- subset(smdata, select = c('speciesid', 'cell', 'WEFFr', 'TSMr'))
names(smdata) <- c('speciesid.am', 'cell', 'WEFFr', 'S.TSMr')
sdata <- left_join(nr, smdata, by = c('speciesid.am', 'cell'))
rm(list = c('smdata'))
gc()

## ADD REDLIST STATUS; RLSTATUS + WORMS = SDATA
load(paste(datadirout, 'S_redlist_1deg.RData', sep=''))
rl <- subset(rl, select = c('speciesid', 's.rlstatus'))
names(rl) <- c('speciesid.am', 'S.rlstatus')
rl$S.rlstatus <- as.numeric(as.character(rl$S.rlstatus)) ## Convert factor to numeric if necessary
sdata <- left_join(sdata, rl, by = c('speciesid.am'))
rm(rl)
gc()

## ADD VERTICAL RANGE DATA
load(paste(datadirout, 'S_vrange_1deg.RData', sep=''))
vdata <- subset(vdata, select = c('speciesid', 'lon', 'lat', 'vertrange', 'depthmax'))
names(vdata) <- c('speciesid.am', 'lon', 'lat', 'S.vertrange', 'S.depthmax')
sdata <- left_join(sdata, vdata, by = c('speciesid.am', 'lon', 'lat'))
rm(list = c('vdata'))
gc()

## ADD HUMAN IMPACT INDEX (HII)
load(paste(datadirout, 'S_HII_1deg.RData', sep=''))
hiidata$cell <- gsub(' ', '', paste(hiidata$lon, '_', hiidata$lat))
hiidata <- subset(hiidata, select = c('cell', 'hii'))
names(hiidata) <- c('cell', 'S.HII')
sdata <- left_join(sdata, hiidata, by = c('cell'))
rm(list = c('hiidata'))
gc()

##############################################################
## ADAPTIVE CAPACITY; HIGHER VALUES INDICATE HIGHER ADAPTABILITY
##############################################################

## LOAD ADAPTIVE CAPACITY DATA
load(paste(datadirout, 'AC_SST_AUC_1deg.RData', sep=''))
dout <- subset(dout, select = c('speciesid', 'cell', 'auctpct.r', 'trng.r'))
names(dout) <- c('speciesid.am', 'cell', 'AC.tauc', 'AC.trng')
acdata <- left_join(nr, dout, by = c('speciesid.am', 'cell'))
rm(list = c('dout'))
gc()

## HABITAT FRAGMENTATION
load(paste(datadirout, 'AC_hfrag_1deg.RData', sep=''))
hfrag <- subset(hfrag, select = c('speciesid', 'hfrag'))
names(hfrag) <- c('speciesid.am', 'AC.hfrag')
acdata <- left_join(acdata, hfrag, by = c('speciesid.am'))
rm(list = c('hfrag'))
gc()

## ADD HABITAT RANGE
load(paste(datadirout, 'AC_hrange_1deg.RData', sep=''))
names(hrange) <- ifelse(names(hrange) %in% c('hrange'), 'rarea', names(hrange))
hrange <- subset(hrange, select = c('speciesid', 'rarea', 'lrange'))
names(hrange) <- c('speciesid.am', 'AC.rarea', 'AC.lrange')
acdata <- left_join(acdata, hrange, by = c('speciesid.am'))
rm(list = c('hrange'))
gc()

## ADD MAX SPECIES LENGTH
load(paste(datadirout, 'AC_maxlengths_1deg.RData', sep=''))
ldat <- subset(ldat, select = c('speciesid', 'maxtl'))
names(ldat) <- c('speciesid.am', 'AC.lmax')
acdata <- left_join(acdata, ldat, by = c('speciesid.am'))
rm(list = c('ldat'))
gc()

##############################################################
## COMBINE RAW INDICES
##############################################################

## CONVERT TO DATA TABLES FOR FASTER PROCESSING
sdata <- as.data.table(sdata)
edata8 <- as.data.table(edata8)
acdata <- as.data.table(acdata)

## REMOVE EXTRANEOUS VARIABLES
sdata <- sdata[, c('assessmentid.rl', 'kingdom', 'phylum', 'order', 'family', 'cname', 'acceptedname') := NULL]
acdata <- acdata[, c('lon', 'lat', 'prob') := NULL]
edata8 <- edata8[, c('lon', 'lat', 'prob') := NULL]

## VALIDATE DIMENSIONS OF DATASETS
print(dim(sdata))
print(dim(acdata))
print(dim(edata8))

## JOIN AND SAVE RAW INDICES FOR DIFFERENT CLIMATE SCENARIOS
vdatar <- join(sdata, acdata, by = c('speciesid.am', 'cell'))
vdatar <- join(vdatar, edata8, by = c('speciesid.am', 'cell'))
save(vdatar, file = paste(datadirout, 'VULNERABILITY_RAW_1deg_RCP8.5.RData', sep=''))
print(summary(vdatar))
gc()

## PROCESS RCP 2.6 DATA
edata2 <- as.data.table(edata2)
edata2 <- edata2[, c('lon', 'lat', 'prob') := NULL]
vdatar2 <- join(sdata, acdata, by = c('speciesid.am', 'cell'))
vdatar2 <- join(vdatar2, edata2, by = c('speciesid.am', 'cell'))
save(vdatar2, file = paste(datadirout, 'VULNERABILITY_RAW_1deg_RCP2.6.RData', sep=''))
print(summary(vdatar2))
gc()






stdfun<-function(d){
  
  # Function to compute and visualize summary statistics
  sfun<-function(a){
    print(mean(a,na.rm=TRUE))  # Print mean of the vector, ignoring NAs
    hist(a,xlim=c(0,1))        # Plot histogram with limits [0,1]
  }

  ## ADAPTIVE CAPACITY
  # Transform habitat fragmentation data
  d$AC.hfrag<-tf.AC.hfrag(d$AC.hfrag)
  print(mean(d$AC.hfrag,na.rm=TRUE))
  
  # Transform and summarize habitat range data
  d$AC.rarea<-tf.AC.rarea(d$AC.rarea)
  d$AC.lrange<-tf.AC.lrange(d$AC.lrange)
  d$AC.hrange<-tf.AC.hrange(d$AC.rarea,d$AC.lrange)
  print(summary(d$AC.hrange))

  # Transform and summarize body length data
  d$AC.lmax<-tf.AC.lmax(d$AC.lmax)
  print(summary(d$AC.lmax))

  # Transform and summarize temperature adaptation data
  d$AC.tauc<-tf.AC.tauc(d$AC.tauc)
  d$AC.trng<-tf.AC.trng(d$AC.trng)
  d$AC.tvar<-tf.AC.tvar(d$AC.tauc,d$AC.trng)
  print(summary(d$AC.tvar))

  ## SENSITIVITY
  # Transform and summarize Ocean Health Index
  d$S.HII<-tf.S.HII(d$S.HII)
  print(summary(d$S.HII))

  # Transform and summarize vertical habitat data
  d$S.depthmax<-tf.S.depthmax(d$S.depthmax)
  d$S.vertrange<-tf.S.vertrange(d$S.vertrange)
  d$S.vind<-tf.S.vind(d$S.vertrange,d$S.depthmax)
  print(summary(d$S.vind))

  # Transform and summarize Thermal Safety Margin
  d$S.TSMr<-tf.S.TSMr(d$S.TSMr)
  print(summary(d$S.TSMr))

  ## CLIMATE VELOCITY
  d$E.vel<-tf.E.vel(d$E.vel)
  print(summary(d$E.vel))

  ## PROPORTION OF SPECIES LOST DUE TO CLIMATE CHANGE
  d$E.plost<-tf.E.plost(d$E.plost)
  print(summary(d$E.plost))

  ## PROJECTED TIME OF EMERGENCE FOR SPECIES
  d$E.toe<-tf.E.toe(d$E.toe)
  print(summary(d$E.toe))

  ## PROJECTED CHANGE IN NATIVE RANGE
  d$E.nrchng<-tf.E.nrchng(d$E.nrchng)
  print(summary(d$E.nrchng))

  print('indices done')

  ######################################
  ## VULNERABILITY DIMENSIONS
  RowSD <- function(x) {
    sqrt(rowSums((x-rowMeans(x,na.rm=TRUE))^2)/(dim(x)[2]-1))
  }

  ## Compute Sensitivity Score
  d$sens<-rowMeans(subset(d,select=c('AC.hrange','S.rlstatus','S.TSMr','S.vind')),na.rm=TRUE)
  d$sens.sd<-RowSD(subset(d,select=c('AC.hrange','S.rlstatus','S.TSMr','S.vind')))

  ## Compute Missing Values in Sensitivity
  d$n.sens<-rowSums(is.na(subset(d,select=c('S.HII','S.rlstatus','S.TSMr','S.vind'))),na.rm=TRUE)
  
  ## Compute Adaptive Capacity Score
  d$adcap<-rowMeans(subset(d,select=c('AC.tvar','AC.hfrag','S.HII','AC.lmax')),na.rm=TRUE)
  d$adcap.sd<-RowSD(subset(d,select=c('AC.tvar','AC.hfrag','S.HII','AC.lmax')))

  ## Compute Missing Values in Adaptive Capacity
  d$n.adcap<-rowSums(is.na(subset(d,select=c('AC.tvar','AC.hfrag','AC.hrange','AC.lmax'))),na.rm=TRUE)

  ## Compute Exposure Score
  d$expo<-rowMeans(subset(d,select=c('E.toe','E.vel','E.plost','E.nrchng')),na.rm=TRUE)
  d$expo.sd<-RowSD(subset(d,select=c('E.toe','E.vel','E.plost','E.nrchng')))

  ## Compute Missing Values in Exposure
  d$n.expo<-rowSums(is.na(subset(d,select=c('E.toe','E.vel','E.plost','E.nrchng'))),na.rm=TRUE)

  print('dimensions done')

  #######################################################
  ## AVERAGE VULNERABILITY
  dr<-discrate(80,12,40,0.025)  # Compute discount rate
  print(dr)

  ## Adjust Adaptive Capacity so that higher values indicate higher vulnerability
  d$adcap<-1-d$adcap

  ## Apply discount rate
  d$expo.dr<-(1-dr)*d$expo
  d$sens.dr<-(1+dr)*d$sens
  d$adcap.dr<-1*d$adcap

  ## Compute Unweighted Vulnerability
  d$vuln.avg<-(d$expo + d$adcap + d$sens)/3
  print('average vuln done')

  ## Compute Discounted Vulnerability
  dsd<-setDT(d)[,list(vuln.uw=mean(c(expo.dr,sens.dr,adcap.dr),na.rm=TRUE),
                      vulnsd.uw=sd(c(expo.dr,sens.dr,adcap.dr),na.rm=TRUE)),by=.(speciesid.am,cell)]
  d<-left_join(d,dsd,by=c('speciesid.am','cell'))

  ## Convert Data to Tall Format for Weighted Vulnerability Calculation
  d2<-subset(d,select=c('speciesid.am','cell','expo.dr','sens.dr','adcap.dr'))
  dd2<-gather(d2,key='dimension',value='x',-c('speciesid.am','cell'))
  dd2$dimension<-gsub('.dr','',dd2$dimension)
  print('uw vuln done')

  ## Compute Weighted Vulnerability
  d3<-setDT(d2)[,list(vuln.w=wtd.mean(x,xcvi,na.rm=TRUE),
                      vulnsd.w=sqrt(wtd.var(x,xcvi,na.rm=TRUE,normwt=TRUE))),by=.(speciesid.am,cell)]
  print('w vuln done')

  ## Join Weighted to Unweighted Vulnerability
  d<-left_join(d,d3,by=c('speciesid.am','cell'))
  print('return done')
  return(d)
}

## Run Analysis for RCP8.5 and Save Output
vdata<-stdfun(vdatar)
save(vdata,file=paste(datadirout,'VULNERABILITY_1deg_RCP8.5.RData',sep=''))

## Run Analysis for RCP2.6 and Save Output
vdata2<-stdfun(vdatar2)
save(vdata2,file=paste(datadirout,'VULNERABILITY_1deg_RCP2.6.RData',sep=''))


