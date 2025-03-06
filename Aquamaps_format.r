## LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
# Load parameter file containing necessary settings and configurations
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

# Function to set the directory path based on the current working directory
cdir<-function(){
    codedir<-ifelse(grepl('sailfish',getwd(), fixed = TRUE),
                    'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                    'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
    source(paste(codedir,'GPARAMs_CC_Vuln.r',sep=''))
}

# Function for processing Aquamaps data
efun<-function(){
    print('Aquamaps_format_start')
    
    ##############################################################
    ## RASTER STACK OF ALL NATIVE RANGES
    ##############################################################
    
    # Load Aquamaps thermal envelope information
    hspen<-fread(paste(datadir,'Aquamaps/2020/raw_data/hspen.csv',sep=''),header=TRUE)
    names(hspen)<-tolower(names(hspen))
    hspen<-subset(hspen,select=c('speciesid','speccode','rank','depthmin'))
    
    # Read 1/2-degree native range data and re-process to 1-degree resolution
    nr<-fread(paste(datadir,'Aquamaps/2020/raw_data/hcaf_species_native.csv',sep=''),header=TRUE)
    names(nr)<-tolower(names(nr))
    names(nr)[3:4]<-c('lat','lon')
    nr$cell<-gsub(' ','',paste(nr$lon,'_',nr$lat))
    nr<-subset(nr,speciesid %in% spinclude$speciesid)
    
    gc()  # Run garbage collection to free up memory
    
    # Merge with species thermal envelope data
    nr<-left_join(nr,hspen,by=c('speciesid'))
    
    # Define coordinate reference system (CRS)
    mct<- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    
    # Create a global 1-degree raster template
    crds<-expand.grid(lon=seq(-179.5,179.5,1), lat=seq(-89.5,89.5,1), probability=NA)
    r1deg<-rasterFromXYZ(crds)
    crs(r1deg) <- mct
    
    # Adjust resolution for processing
    crds<-expand.grid(lon=seq(-179.75,179.75,0.5), lat=seq(-89.75,89.75,0.5), probability=NA)
    crds$speciesid<-'aaa'
    
    # Remove the 'cell' column from nr
    nr<-nr[,!(names(nr) %in% c('cell'))]
    
    # Split species IDs into smaller groups for batch processing
    spids<-sort(unique(nr$speciesid))
    spidsl<-split(spids, ceiling(seq_along(spids)/500))
    
    # Function to create rasters at native resolution and interpolate to 1-degree resolution
    f<-function(a){
        d<-subset(nr,speciesid %in% a)
        d<-rbind.fill(crds,data.frame(d))
        d<-spread(d,key=speciesid,value=probability)
        d<-SpatialPixelsDataFrame(points=d[c('lon','lat')],data=d)
        r<-stack(d)
        rm(d)
        r<-dropLayer(r,c(1,2,3))
        crs(r) <- mct
        
        # Resample to 1-degree resolution using nearest neighbor method
        r1<-resample(r, r1deg, method='ngb')
        writeRaster(r1,paste(datadir,'analysis_outputs/one_deg/temp/1deg/chunksnn/',stri_rand_strings(1,length=5),'.grd',sep=''),format='raster')
        
        # Resample to 1-degree resolution using bilinear interpolation
        r2<-resample(r, r1deg, method='bilinear')
        writeRaster(r1,paste(datadir,'analysis_outputs/one_deg/temp/1deg/chunksbl/',stri_rand_strings(1,length=5),'.grd',sep=''),format='raster')
        
        rm(list=c('r','r1','r2'))
        gc()
    }
    
    # Apply function to each species subset
    l<-llply(spidsl,.fun=f,.progress='text')
    
    ## STACK NATIVE 1-DEGREE RASTERS AND SAVE OUTPUTS
    
    # Process nearest neighbor interpolated rasters
    setwd(paste(datadir,'analysis_outputs/one_deg/temp/1deg/chunksnn/',sep=''))
    lf<-list.files()
    lf<-str_subset(lf,'grd')
    system.time(rdata<-stack(lf))
    writeRaster(rdata,filename=paste(datadir,'analysis_outputs/rasters/hcaf_species_native_1deg_rstack_NN.grd',sep=''),format='raster',overwrite=TRUE,progress='text')
    save(rdata,file=paste(datadirout,'hcaf_species_native_1deg_rstack_NN.RData',sep=''))
    print('NN write done')
    
    # Process bilinear interpolated rasters
    setwd(paste(datadir,'analysis_outputs/one_deg/temp/1deg/chunksbl/',sep=''))
    lf<-list.files()
    lf<-str_subset(lf,'grd')
    system.time(rdata<-stack(lf))
    writeRaster(rdata,filename=paste(datadir,'analysis_outputs/one_deg/rasters/hcaf_species_native_1deg_rstack.grd',sep=''),format='raster',overwrite=TRUE,progress='text')
    save(rdata,file=paste(datadirout,'hcaf_species_native_1deg_rstack.RData',sep=''))
    print('BL write done')
    
    ## CONVERT PROBABILITIES TO BINARY VALUES (0=NA, >0=1) FOR RANGE FRAGMENTATION ANALYSIS
    replace0 <- function(x){ return(ifelse(x==0,NA,1)) }
    system.time(rdatab<-calc(rdata,fun=replace0))
    names(rdatab)<-names(rdata)
    writeRaster(rdatab,paste(datadir,'analysis_outputs/one_deg/rasters/hcaf_species_native_1deg_rstack_BINARY.grd',sep=''),format='raster',overwrite=TRUE,progress='text')
    save(rdatab,file=paste(datadirout,'hcaf_species_native_1deg_rstack_BINARY.RData',sep=''))
    print('binary done')
    
    ## CREATE LONG-FORM DATABASE
    rdata<-stack(paste(datadir,'analysis_outputs/one_deg/rasters/hcaf_species_native_1deg_rstack.grd',sep=''))
    nr<-list()
    for(i in seq(1,dim(rdata)[3],1)){
        d<-subset(rdata,i)
        nm<-names(d)
        xy<-data.frame(coordinates(d))
        d<-data.frame(speciesid=nm, lon=xy$x, lat=xy$y, probability=getValues(d))
        nr[[i]]<-subset(d,is.na(probability)==FALSE & probability>0)
    }
    nr<-data.frame(rbindlist(nr))
    nr$cell<-gsub(' ','',paste(nr$lon,'_',nr$lat))
    nr$speciesid<-gsub('\.','\-',nr$speciesid)
    save(nr,file=paste(datadirout,'hcaf_species_native_1deg.RData',sep=''))
    
    rm(list=ls(all=TRUE))
    gc()
    print('Aquamaps_format_done')
}


