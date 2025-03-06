################################################################################################
# R CODE FOR: BOYCE, D. ET AL. 2021. A CLIMATE RISK INDEX FOR MARINE LIFE
# CODE AND ANALYSES DEVELOPED BY DANIEL BOYCE (DBOYCE@DAL.CA); PLEASE REFERENCE USE ACCORDINGLY
# IMPLEMENTED IN R VERSION 4.04 AND RUN ON COMPUTE CANADA COMPUTING CLUSTER
# LAST MODIFIED: AUGUST 2021
################################################################################################

# LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')
rasterOptions(progress="text", timer=TRUE)

# Function to set the working directory based on system name
cdir <- function() {
    codedir <- ifelse(grepl('sailfish', getwd(), fixed = TRUE),
                      'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                      'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
    source(paste(codedir, 'GPARAMs_CC_Vuln.r', sep = ''))
}
# cdir()  # Uncomment to use

#############################################################
# THERMAL SAFETY MARGINS FOR INDIVIDUAL SPECIES (SPECIES-SPATIAL)
#############################################################

# IMPORTS DAILY MAX SEA SURFACE TEMPERATURE (SST) DATA FROM PATHFINDER AND REYNOLDS
# GETS MAX VALUE PER 1-DEGREE CELL AND MERGES DATA

# Load Pathfinder SST data
(load(paste(datadir, 'analysis_outputs_local/one_deg/rfiles/SST/Pathfinder_SST_daily_1deg_MaxSST_2010_2020.RData', sep = '')))
names(pdatap) <- c('lon', 'lat', 'sstp', 'cell')
print('Pathfinder data loaded')

# Load Reynolds SST data
(load(paste(datadir, 'analysis_outputs_local/one_deg/rfiles/SST/Reynolds_OISST_daily_1deg_MaxSST_2010_2020.RData', sep = '')))
names(pdatar) <- c('lon', 'lat', 'sstr', 'cell')
print('Reynolds data loaded')

# Merge both SST datasets
pdata <- merge(pdatap, pdatar, by = c('lon', 'lat'), all = TRUE)

##############################################################
# THERMAL PREFERENCES FOR ALL SPECIES - FROM AQUAMAPS ENVELOPES
##############################################################

# Load environmental envelope data for species temperature preferences
hspen <- fread(paste(datadir, 'Aquamaps/2020/raw_data/hspen.csv', sep = ''), header = TRUE)
print('Environmental data loaded')
names(hspen) <- tolower(names(hspen))

# Filter species to include only those of interest
hspen <- subset(hspen, speciesid %in% spinclude$speciesid)

#############################################################
# ESTIMATE EFFECT OF WARMING BASED ON TEMPERATURE PROBABILITIES
#############################################################

# Extract relevant temperature preference values
z <- subset(hspen, select = c('speciesid', 'tempmin', 'tempprefmin', 'tempprefmax', 'tempmax'))
tdat <- gather(z, key = sstc, y, -speciesid)

# Assign probability values based on temperature range
tdat$prob <- ifelse(tdat$sstc == 'tempmin', 0, 1)
tdat$prob <- ifelse(tdat$sstc == 'tempprefmin', 0.1, tdat$prob)
tdat$prob <- ifelse(tdat$sstc == 'tempprefmax', 0.9, tdat$prob)

# Function to estimate temperature optima
f <- function(d) {
    mod <- glm(prob ~ y, family = quasibinomial('logit'), data = d)
    pdat <- data.frame(y = seq(min(d$y) - 5, max(d$y) + 5, 0.01))
    pdat$p <- predict(mod, newdata = pdat, type = 'response')
    d2 <- subset(pdat, p >= 0.4 & p <= 0.6)
    dout <- data.frame(speciesid = unique(d$speciesid),
                       tempoptmin = min(d2$y),
                       tempoptmax = max(d2$y),
                       tempopt = median(d2$y))
    dout$tempoptmin <- ifelse(dout$tempoptmin %in% c(Inf, -Inf), NA, dout$tempoptmin)
    dout$tempoptmax <- ifelse(dout$tempoptmax %in% c(Inf, -Inf), NA, dout$tempoptmax)
    return(dout)
}

# Apply function across species dataset
l <- dlply(tdat, .(speciesid), .fun = f, .progress = 'text')
tdata <- data.frame(rbindlist(l))

# Merge temperature preference data with species IDs
a <- subset(hspen, select = c('speciesid', 'tempmin', 'tempprefmin', 'tempprefmax', 'tempmax'))
tdata <- merge(tdata, a, by = c('speciesid'), all.x = TRUE, all.y = FALSE)
tdata$speciesid <- as.character(tdata$speciesid)

# Save processed temperature preference data
save(tdata, file = paste(datadirout, 'thermal_prefs_byspecies.RData', sep = ''))
print('Thermal preferences written')

########################################################
# ALL SPECIES RANGES AT 1 DEGREE
########################################################

# Load species distribution data at 1-degree resolution
(load(paste(datadir, 'Aquamaps/2020/processed_data/hcaf_species_native_1deg.RData', sep = '')))
print('Species ranges loaded')

# Filter dataset to include only selected species
nr <- subset(nr, speciesid %in% spinclude$speciesid)

#######################################################
# MERGE SPECIES RANGES WITH TEMPERATURE TOLERANCES
#######################################################

# First join species ranges with estimated thermal preferences
smdata <- left_join(nr, tdata, by = c('speciesid'))
print("First join completed")

# Merge with daily maximum SST data
smdata <- left_join(smdata, pdata, by = c('lon', 'lat'))
print("Second join completed")

# Compute Thermal Safety Margins (TSM) - positive values indicate species are at risk
smdata$TSMp <- smdata$tempmax - smdata$sstp
smdata$TSMp <- ifelse(smdata$TSMp < 0, 0, smdata$TSMp)
smdata$TSMr <- smdata$tempmax - smdata$sstr
smdata$TSMr <- ifelse(smdata$TSMr < 0, 0, smdata$TSMr)

# Compute Warming Effect Flags (WEFF) indicating warming impact
smdata$WEFFr <- ifelse(smdata$sstr < smdata$tempoptmin, 1, 0)
smdata$WEFFr <- ifelse(smdata$sstr > smdata$tempoptmax, -1, smdata$WEFFr)
smdata$WEFFp <- ifelse(smdata$sstp < smdata$tempoptmin, 1, 0)
smdata$WEFFp <- ifelse(smdata$sstp > smdata$tempoptmax, -1, smdata$WEFFp)

# Save final dataset
save(smdata, file = paste(datadirout, 'S_TSM_1deg.RData', sep = ''))
print("Thermal safety margins computed and saved")

