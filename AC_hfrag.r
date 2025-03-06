################################################################################################
## R CODE FOR: BOYCE, D. ET AL. 2021. A CLIMATE RISK INDEX FOR MARINE LIFE
## CODE AND ANALYSES DEVELOPED BY DANIEL BOYCE (DBOYCE@DAL.CA); PLEASE REFERENCE USE ACCORDINGLY
## IMPLEMENTED IN R VERSION 4.04 AND RUN ON COMPUTE CANADA COMPUTING CLUSTER
## LAST MODIFIED: AUGUST 2021
################################################################################################

## LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')  # Load parameter settings and dependencies

## FUNCTION TO SET CODE DIRECTORY PATH BASED ON WORKING DIRECTORY
cdir<-function(){
    codedir<-ifelse(grepl('sailfish', getwd(), fixed = TRUE),  # Check if 'sailfish' is in the working directory path
                    'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',  # Set path for Sailfish user
                    'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')  # Set path for Daniel's user
    source(paste(codedir,'GPARAMs_CC_Vuln.r',sep=''))  # Source the parameter file
}

##############################################################
## CHANGE IN HABITAT FRAGMENTATION AND GEOGRAPHIC RANGE
##############################################################

## LOAD RASTER DATA CONTAINING ALL NATIVE RANGES
rdatab <- stack(paste(rasterdir,'hcaf_species_native_1deg_rstack_BINARY.grd',sep=''))  # Load multi-layer raster data (species native ranges)

## SET PROJECTION TO WGS84 (MERCATOR PROJECTION)
crs(rdatab) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

############################################
## CALCULATE HABITAT FRAGMENTATION FOR EACH LAYER
############################################

## CONFIGURE RASTER PROCESSING OPTIONS
## rasterOptions(tmpdir="C:/Users/sailfish/Downloads", tmptime = 24, progress="text", timer=TRUE, overwrite = T, chunksize=2e8, maxmemory=1e8)  # Alternative settings
rasterOptions(progress="text", timer=TRUE)  # Enable progress updates and timing for raster operations

## COMPUTE HABITAT FRAGMENTATION
system.time(hfrag <- lsm_l_np(rdatab, directions=8))  # Compute number of patches (landscape metric), takes ~17230 sec
hfrag$speciesid <- names(rdatab)  # Assign species names to the fragmentation data
hfrag$speciesid <- gsub('\\.', '\\-', hfrag$speciesid)  # Replace dots with dashes in species IDs
hfrag <- data.frame(hfrag)  # Convert to dataframe

## CALCULATE TOTAL CELLS COMPRISING THE RANGE OF EACH SPECIES
dcell <- cellStats(rdatab, 'sum')  # Compute the total number of occupied cells per species (~819 sec runtime)
dcell <- data.frame(speciesid = names(dcell), ncells = unlist(dcell))  # Convert results to dataframe
dcell$speciesid <- gsub('\\.', '\\-', dcell$speciesid)  # Ensure species ID format is consistent

## MERGE FRAGMENTATION DATA WITH CELL COUNT DATA
hfrag <- left_join(hfrag, dcell, by=c('speciesid'))  # Merge fragmentation data with cell count data
hfrag$hfrag <- (hfrag$value / hfrag$ncells)  # Compute habitat fragmentation index

## SELECT AND RENAME COLUMNS FOR FINAL OUTPUT
hfrag <- subset(hfrag, select = c('speciesid', 'value', 'ncells', 'hfrag'))  # Retain relevant columns
names(hfrag) <- c('speciesid', 'nfrags', 'ncells', 'hfrag')  # Rename columns for clarity

## SAVE OUTPUT AS RDATA FILE
save(hfrag, file = paste(datadirout, 'AC_hfrag_1deg.RData', sep=''))  # Save processed data

## CLEAN UP MEMORY AND PRINT COMPLETION MESSAGE
rm(list=ls(all=TRUE))  # Clear all objects from the environment
gc()  # Perform garbage collection to free memory
print('AC_hfrag_done')  # Print message to indicate successful completion




