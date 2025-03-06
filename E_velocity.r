################################################################################################
## R CODE FOR: BOYCE, D. ET AL. 2021. A CLIMATE RISK INDEX FOR MARINE LIFE
## CODE AND ANALYSES DEVELOPED BY DANIEL BOYCE (DBOYCE@DAL.CA); PLEASE REFERENCE USE ACCORDINGLY
## IMPLEMENTED IN R VERSION 4.04 AND RUN ON COMPUTE CANADA COMPUTING CLUSTER
## LAST MODIFIED: AUGUST 2021
################################################################################################

# Load required parameters and configurations
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

# Set the directory path dynamically based on the working directory
codedir <- ifelse(grepl('sailfish', getwd(), fixed = TRUE),
                  'C:/Users/sailfish/Documents/aalldocuments/literature/research/complete/2022_CC_risk/code/',
                  'C:/Users/danie/Documents/aalldocuments/literature/research/complete/2022_CC_risk/code/')
source(paste(codedir, 'GPARAMs_CC_Vuln.r', sep = ''))

print('E velocity start')

#######################################################
# Process monthly SST projections and calculate temperature change and velocity for SSP 585 scenario
#######################################################

# Define data directory
pn <- paste(datadir, 'henson_TOE/Raw_CMIP_data/2021/SSP585/matlab/', sep = '')

# List files in the directory, excluding certain models
fls <- list.files(pn)
fls <- fls[!(fls %in% c('sstCESMW85.mat', 'sstCESMW26.mat', "sstIPSL85.mat", "sstCanESM85.mat", "sstIPSL26.mat", "sstCanESM26.mat"))]
print(fls)

# Initialize list to store data
l <- list()

# Loop through each file and process SST data
for (i in 1:length(fls)) {
    pathname <- paste(pn, fls[i], sep = '')
    data <- readMat(pathname)  # Read .mat file
    names(data) <- 'x'
    data <- data$x
    print(fls[i])
    print(dim(data))
    
    # Convert data to raster format
    bdata <- brick(data, xmn = 0, xmx = 360, ymn = -90, ymx = 90, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    print(bdata)
    
    # Calculate yearly averages from monthly data
    r <- sumSeries(bdata, p = "2015-01/2100-12", yr0 = "2015-01-01", l = nlayers(bdata), fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years")
    
    # Compute temporal trend of SST
    vt <- tempTrend(r, th = 10)
    
    # Compute spatial gradient (degrees per km and angle)
    vg <- spatGrad(r, th = 0.0001, projected = FALSE)
    
    # Calculate climate velocity and its magnitude (km/year) and angle (degrees)
    gv <- gVoCC(vt, vg)
    
    # Convert raster data to data frames
    dfvel <- data.frame(rasterToPoints(gv))
    df <- data.frame(rasterToPoints(vt))
    df$vel <- dfvel$voccMag
    df$ang <- dfvel$voccAng
    df$y <- df$y * -1  # Invert latitude
    df$model <- gsub('.mat', '', paste(fls[i]))
    l[[i]] <- df
}

# Combine data from all models
dout <- data.frame(rbindlist(l))

# Compute multi-model average
mmfun <- function(d) {
    l <- data.frame(sstb = mean(d$slpTrends, na.rm = TRUE),
                    sstbse = sd(d$slpTrends, na.rm = TRUE),
                    sstvel = mean(d$vel, na.rm = TRUE),
                    sstvelse = sd(d$vel, na.rm = TRUE),
                    sstangd = mean(circular(d$ang, units = 'degrees', zero = 0, rotation = 'clock')),
                    nmods = length(unique(d$mod)))
    l$ang <- abs(abs(l$sstangd) - 360)
    return(l)
}

dvel <- ddply(dout, .(x, y), .fun = mmfun, .progress = 'text')
names(dvel)[1:2] <- c('lon', 'lat')
dvel <- subset(dvel, select = c('lon', 'lat', 'sstvel', 'sstvelse'))
dvel$lon <- ifelse(dvel$lon > 180, -360 + dvel$lon, dvel$lon)  # Convert longitudes

# Save processed data
save(dvel, file = paste(datadirout, 'E_SSTvel_1deg_85.RData', sep = ''))
print(summary(dvel))
rm(dvel)
gc()

#############################################################
# Repeat the process for SSP 126 scenario
#############################################################

pn <- paste(datadir, 'henson_TOE/Raw_CMIP_data/2021/SSP126/matlab/', sep = '')
fls <- list.files(pn)
fls <- fls[!(fls %in% c('sstCESMW85.mat', 'sstCESMW26.mat', "sstIPSL85.mat", "sstCanESM85.mat", "sstIPSL26.mat", "sstCanESM26.mat"))]
print(fls)

l <- list()

for (i in 1:length(fls)) {
    pathname <- paste(pn, fls[i], sep = '')
    data <- readMat(pathname)
    names(data) <- 'x'
    data <- data$x
    print(fls[i])
    print(dim(data))
    
    bdata <- brick(data, xmn = 0, xmx = 360, ymn = -90, ymx = 90, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    print(bdata)
    
    r <- sumSeries(bdata, p = "2015-01/2100-12", yr0 = "2015-01-01", l = nlayers(bdata), fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years")
    
    vt <- tempTrend(r, th = 10)
    vg <- spatGrad(r, th = 0.0001, projected = FALSE)
    gv <- gVoCC(vt, vg)
    
    dfvel <- data.frame(rasterToPoints(gv))
    df <- data.frame(rasterToPoints(vt))
    df$vel <- dfvel$voccMag
    df$ang <- dfvel$voccAng
    df$y <- df$y * -1
    df$model <- gsub('.mat', '', paste(fls[i]))
    l[[i]] <- df
}

dout <- data.frame(rbindlist(l))

save(dvel, file = paste(datadirout, 'E_SSTvel_1deg_26.RData', sep = ''))
print(summary(dvel))
rm(list = ls(all = TRUE))
gc()
print('E_velocity_done')