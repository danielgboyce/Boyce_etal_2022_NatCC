## LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
# Load external R script containing parameter settings
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

# Function to set code directory dynamically based on the working directory
cdir <- function() {
  codedir <- ifelse(grepl('sailfish', getwd(), fixed = TRUE),
                    'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                    'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
  source(paste(codedir, 'GPARAMs_CC_Vuln.r', sep = ''))
}
## cdir()  # Uncomment to execute the function if needed


## FUNCTION TO CALCULATE TIME OF EMERGENCE (ToE) FOR SPECIES TEMPERATURE THRESHOLDS
# This function determines the first year when a species experiences sustained exposure
# to temperatures exceeding its maximum tolerance for at least 'runLength' years
findToE <- function(focmaxT, futTs, yearsFut, runLength = 5) {
  # Default ToE value (if no emergence occurs by the end of the projection period)
  ToEFirst <- max(yearsFut) + 1
  
  # Identify years where future temperatures exceed the species' maximum tolerance
  emYears <- which(futTs > focmaxT)
  
  if (length(emYears) > 0) {
    x <- rep(0, length(yearsFut))  # Initialize binary vector (1 = emergent, 0 = non-emergent)
    x[emYears] <- 1  # Mark years of emergence
    
    # Identify runs of consecutive '1's that meet the threshold for emergence
    diffs <- x[-1L] != x[-length(x)]
    idx <- c(which(diffs), length(x))
    runs <- diff(c(0, idx))
    x2 <- x[-which(c(NA, diff(x)) == 0)]
    yearsRunStart <- yearsFut[-which(c(NA, diff(x)) == 0)]
    runsTab <- data.frame(year = yearsRunStart, val = x2, run = runs)
    runsTab <- runsTab[runsTab$val == 1, ]  # Keep only runs of '1's
    runStart <- which(runsTab$run >= runLength)
    
    if (length(runStart) > 0) {
      ToEFirst <- min(runsTab$year[runStart])  # Earliest emergence year meeting criteria
    }
  }
  return(ToEFirst)
}


## LOAD SPECIES TEMPERATURE TOLERANCE DATA
# Read in species-specific temperature tolerance data
hspen <- fread(paste(datadir, 'Aquamaps/2020/raw_data/hspen.csv', sep = ''), header = TRUE)
names(hspen) <- tolower(names(hspen))
hspen <- subset(hspen, speciesid %in% spinclude$speciesid)
hspen <- subset(hspen, select = c('speciesid', 'tempmax'))


## LOAD SPECIES RANGE DATA
# Load habitat suitability data and merge with species temperature tolerance data
load(paste(datadirout, 'hcaf_species_native_1deg.RData', sep = ''))
nr <- subset(nr, speciesid %in% spinclude$speciesid)
nr <- left_join(nr, hspen, by = c('speciesid'))

# Define the global domain of interest
# Extract unique cells (grid locations) where species are present
domain <- unique(subset(nr, select = c('cell', 'lon', 'lat')))

# Subset data to retain only relevant columns
nr <- subset(nr, select = c('speciesid', 'cell', 'tempmax'))

# Function to structure range data for each species
f <- function(d) {
  return(list('speciesid' = d$speciesid[1], 'cell' = d$cell, 'tmax' = d$tempmax[1]))
}

# Apply function to all species and store in rangeCells
rangeCells <- dlply(nr, .(speciesid), .fun = f, .progress = 'text')

# Clean up workspace
rm(nr)
gc()


## PARALLEL COMPUTATION SETUP
# Load parallel processing libraries
library(doSNOW)
library(foreach)
library(doParallel)

# Detect and set the number of cores for parallel execution
UseCores <- detectCores()
print(UseCores)
UseCores <- 32  # Set to 32 cores for computation
print(UseCores)

# Create and register parallel computing cluster
cl <- makeCluster(UseCores, type = "SOCK")
registerDoSNOW(cl)

## DEFINE DIRECTORY CONTAINING CLIMATE PROJECTIONS
projdir <- paste(datadir, '/analysis_outputs_local/one_deg/rfiles/cmodels/SSP585/annual_mean/', sep = '')

# List available climate projection files
fls <- paste0(projdir, list.files(projdir))
nms <- list.files(projdir)
setwd(projdir)


## PROCESS EACH CLIMATE MODEL PROJECTION IN A LOOP
l <- list()
for (x in seq(n, length(fls), 1)) {
  print(fls[x])
  
  # Load climate projection data
  load(fls[x])
  db <- subset(db, select = c(paste('X', c(seq(2015, 2100, 1)), sep = '')))
  
  # Function to process each species' temperature exposure across its range
  funspec <- function(y) {
    sstts <- subset(db, rownames(db) %in% y$cell)
    l <- list()
    
    for (i in seq(1:dim(sstts)[1])) {
      z <- sstts[i,]
      l[[i]] <- data.frame(cell = rownames(z),
                            ToE = findToE(focmaxT = y$tmax, futTs = z, yearsFut = 2015:2100, runLength = 1))
    }
    return(list(data.frame(rbindlist(l))))
  }
  
  # Execute processing in parallel
  print(system.time(
    toedat <- llply(rangeCells,
                     .fun = funspec,
                     .parallel = TRUE,
                     .progress = 'text',
                     .paropts = list(.packages = c('data.table'),
                                     .export = c('db', 'rangeCells', 'findToE')))
  ))
  
  # Save results to local storage (optional)
  ## save(toedat, file = paste('C:/Users/sailfish/Downloads/rfiles/', 'Sp_ToE_', nms[x], sep = ''))
  ## save(toedat, file = paste(datadir, 'analysis_outputs/rfiles/projected_ToEs/', 'Sp_ToE_', nms[x], sep = ''))
  
  l[[x]] <- toedat
  rm(toedat)
  gc()
}

# Assign names to list elements
names(l) <- nms

# Save final results
save(l, file = paste(datadirout, 'E_ToEdata_1deg_list_85.RData', sep = ''))

# Stop parallel processing cluster
stopCluster(cl)
gc()



## LISTS TO DATAFRAME
## Converts a list `l` into a dataframe format
## Uses nested functions to achieve this transformation
f <- function(d) {
    f2 <- function(d2) { return(data.frame(d2)) }
    return(ldply(d, .fun = f2, .id = 'speciesid'))
}

## Apply function `f` across list `l`, extracting data with progress display
## Replaces `.rda` extension from `model` column and restructures the data
## Spreads data into a wide format based on the model name
## Converts the time of emergence (ToE) data for further analysis
toedat <- ldply(l, .fun = f, .progress = 'text', .id = 'model')
toedat$model <- gsub('\\.rda', '', toedat$model)
toedat <- spread(toedat, key = model, value = ToE)
gc()  # Garbage collection

## CONVERT TO YEARS FROM PRESENT
## Ensures values are within a valid range (not exceeding 2100 and not before 2020)
for (i in seq(1, length(nms), 1)) {
    toedat[, i + 2] <- ifelse(toedat[, i + 2] >= 2100, 2100, toedat[, i + 2])
    toedat[, i + 2] <- ifelse(toedat[, i + 2] < 2020, 2020, toedat[, i + 2])
    toedat[, i + 2] <- toedat[, i + 2] - 2020  # Convert to years from 2020
}

################################################
## EFFECT OF CLIMATE CHANGE ON FUTURE RANGE OF EACH SPECIES
################################################
## Calculates the percentage of native range lost by 2100 for each model
rcdata0 <- toedat[, -2]  # Remove second column

f <- function(a) {
    f2 <- function(b) {
        b <- b[b < 200]  # Filter values below 200
        em <- b[b <= 79]  # Define range loss threshold
        return(data.frame(nrchng = length(em) / length(b)))  # Proportion lost
    }
    return(adply(a[, -1], 2, .fun = f2, .expand = TRUE, .id = 'model'))
}

## Apply function `f` to calculate percentage loss across species
rcdat <- ddply(rcdata0, .(speciesid), .fun = f, .progress = 'text')
rcdat <- spread(rcdat, key = 'model', value = 'nrchng')

## Compute the average range change across models
rcdat$nrchng <- rowMeans(rcdat[, 2:dim(rcdat)[2]], na.rm = TRUE)
rcdata <- subset(rcdat, select = c('speciesid', 'nrchng'))
save(rcdata, file = paste(datadirout, 'E_ToEdata_NR_1deg_85.RData', sep = ''))
rm(rcdata)
gc()  # Garbage collection

################################################
## EFFECT OF CLIMATE CHANGE ON ECOSYSTEM INDEX
################################################
## Calculates the percentage of all species lost in each cell

toedataeco0 <- toedat

f <- function(a) {
    f2 <- function(d) {
        d <- d[d <= 200]  # Filter values below 200
        dlost <- d[d <= 79]  # Define loss threshold
        return(data.frame(plost = length(dlost) / length(d)))  # Proportion lost
    }
    d2 <- adply(a[, -c(1,2)], 2, .fun = f2, .id = 'model')
    return(spread(d2, key = 'model', value = 'plost'))
}

## Apply function `f` across spatial cells
toedateco <- ddply(toedataeco0, .(cell), .fun = f, .progress = 'text')

## Compute the average ecosystem loss across models
toedateco$plost <- rowMeans(toedateco[, 2:dim(toedateco)[2]], na.rm = TRUE)

toedataeco <- subset(toedateco, select = c('cell', 'plost'))
save(toedataeco, file = paste(datadirout, 'E_ToEdataeco_1deg_85.RData', sep = ''))
rm(toedataeco)
gc()

## Compute an adjusted loss measure (commented out)
## toedataeco$e.plost <- (1 - exp(-toedataeco$plost / 0.25))

################################################
## TIME OF EMERGENCE (ToE) CALCULATION
################################################

toedata <- toedat

## Compute average time of emergence (ToE) across models
toedata$ToE <- rowMeans(toedata[3:length(nms) + 2], na.rm = TRUE)
toedata <- subset(toedata, select = c('speciesid', 'cell', 'ToE'))

## Save processed ToE data
save(toedata, file = paste(datadirout, 'E_ToEdata_1deg_85.RData', sep = ''))
rm(toedat)
gc()

## Final print statement indicating completion
print('E_ToEs_species_done 85')

