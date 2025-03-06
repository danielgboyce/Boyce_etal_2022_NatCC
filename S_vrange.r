################################################################################################
## R CODE FOR: BOYCE, D. ET AL. 2021. A CLIMATE RISK INDEX FOR MARINE LIFE
## CODE AND ANALYSES DEVELOPED BY DANIEL BOYCE (DBOYCE@DAL.CA); PLEASE REFERENCE USE ACCORDINGLY
## IMPLEMENTED IN R VERSION 4.04 AND RUN ON COMPUTE CANADA COMPUTING CLUSTER
## LAST MODIFIED: AUGUST 2021
################################################################################################

# LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

# Function to set directory path dynamically based on system being used
cdir <- function(){
  codedir <- ifelse(grepl('sailfish', getwd(), fixed = TRUE),
                    'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                    'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
  source(paste(codedir, 'GPARAMs_CC_Vuln.r', sep = ''))
}
## cdir()  # Uncomment if needed to execute

print('S_vrange_start')  # Indicate the start of vertical range processing

#############################################################
## PROCESS BATHYMETRY DATA FROM GEBCO (GLOBAL BATHYMETRY)
#############################################################

# Load bathymetry data
bath <- read.table(paste(datadir, 'bathymetry/crds4km.bathy.txt', sep = ''), sep = '\t')
names(bath) <- c('lon', 'lat', 'x', 'bathy')
bath <- subset(bath, select = c('lon', 'lat', 'bathy'))
bath <- data.table(bath)

# Round lat/lon values to nearest half-degree grid
bath <- bath[, lat := ceiling(lat) - 0.5]
bath <- bath[, lon := ceiling(lon) - 0.5]

# Compute minimum bathymetry value for each grid cell
system.time(dbath <- bath[, .SD[, min(bathy, na.rm = TRUE)], by = .(lon, lat)])
gc()  # Run garbage collection to free memory

# Rename column and adjust bathymetry values (convert to positive depth values)
names(dbath) <- ifelse(names(dbath) %in% c('V1'), 'mbathy', names(dbath))
dbath$mbathy <- ifelse(dbath$mbathy > 0, 0, dbath$mbathy)
dbath$mbathy <- abs(dbath$mbathy)

#############################################################
## LOAD SPECIES ENVELOPE INFORMATION
#############################################################

# Load species environmental envelopes from AquaMaps
hspen <- fread(paste(datadir, 'Aquamaps/2020/raw_data/hspen.csv', sep = ''), header = TRUE)
names(hspen) <- tolower(names(hspen))  # Standardize column names

# Compute latitudinal range of species
hspen$latrange <- abs(hspen$nmostlat - hspen$smostlat)

# Assign vertical habitat zones based on depth range
hspen$vertzone <- ifelse(hspen$depthmax <= 200, 'Epi', NA)
hspen$vertzone <- ifelse(hspen$depthmax > 200 & hspen$depthmin <= 1000, 'Meso', hspen$vertzone)
hspen$vertzone <- ifelse(hspen$depthmax > 1000, 'Bath', hspen$vertzone)

# Compute vertical range for each species
hspen$vertrange <- hspen$depthmax - hspen$depthmin

# Keep relevant columns only
hspen <- subset(hspen, select = c('speciesid', 'depthmin', 'depthmax', 'vertrange', 'latrange'))

#############################################################
## MERGE SPECIES DATA WITH BATHYMETRY INFORMATION
#############################################################

# Load species range data
load(paste(datadirout, 'hcaf_species_native_1deg.RData', sep = ''))
nr$speciesid <- gsub('\.', '\-', nr$speciesid)  # Standardize species ID format
nr <- subset(nr, speciesid %in% spinclude$speciesid)  # Filter to included species

# Merge species range data with bathymetry and species envelope data
nr <- left_join(nr, dbath, by = c('lon', 'lat'))  # Merge with bathymetry
nr <- left_join(nr, hspen, by = c('speciesid'))  # Merge with species depth range data

# Compute realized vertical range based on maximum possible depth in each grid cell
nr$depthmaxr <- ifelse(nr$depthmax >= nr$mbathy, nr$mbathy, nr$depthmax)
nr$vertranger <- nr$depthmaxr - nr$depthmin  # Calculate vertical range
nr$vertranger <- ifelse(nr$vertranger < 1, 1, nr$vertranger)  # Ensure minimum range of 1

# Select relevant columns and rename for output
vdata <- subset(nr, select = c('speciesid', 'lon', 'lat', 'depthmaxr', 'vertranger', 'latrange'))
names(vdata)[4:5] <- c('depthmax', 'vertrange')

# Save processed data to output file
save(vdata, file = paste(datadirout, 'S_vrange_1deg.RData', sep = ''))

# Clean up environment and free memory
rm(list = ls(all = TRUE))
gc()

print('S_vrange_done')  # Indicate completion of vertical range processing
