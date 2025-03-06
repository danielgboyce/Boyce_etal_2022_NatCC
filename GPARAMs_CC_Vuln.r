##########################################################################################

##R CODE FOR: BOYCE, D. ET AL. 2021. A CLIMATE RISK INDEX FOR MARINE LIFE
##CODE AND ANALYSES DEVELOPED BY DANIEL BOYCE (DBOYCE@DAL.CA); PLEASE REFERENCE USE ACCORDINGLY
##IMPLEMENTED IN R VERSION 4.04 AND RUN ON COMPUTE CANADA COMPUTING CLUSTER
##LAST MODIFIED: AUGUST 2021

##########################################################################################
library(Hmisc)
library(grid)
library(fossil)
library(gstat)
library(sp)
library(maptools)
library(rgeos)
library(mapdata)
library(rgdal)
library(landscapemetrics)
library(viridis)
library(tidyr)
library(gridExtra)
library(plotrix)
library(ggplot2)
library(plyr)
library(data.table)
library(pals)
library(dplyr)
library(maps)
library(VoCC)
library(repmis)
library(circular)
library(abind)
library(DescTools)
library(ncdf4)
library(foreign)
library(stringdist)
library(fuzzyjoin)
library(ggrepel)
library(ggExtra)
library(stringi)
library(stringr)
library(rcompanion)
library(bestNormalize)
library(R.matlab)
library(raster)
library(mice)

##THIS IS FOR 'BASE' ANALYSIS - USING ALL SPECIES PROBABILTIY OF OCCURANCES
datadir<-'/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/data/'
datadirout<-'/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/data/analysis_outputs/one_deg/rfiles/'
rasterdir<-'/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/data/analysis_outputs/one_deg/rasters/'

##CRS
mct<- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
##GLOBAL 1DEG COORDS
crds<-expand.grid(lon=seq(-179.75,179.75,0.5),
                  lat=seq(-89.75,89.75,0.5))


## IDENTIFY SPECIES TO REMOVE
## Load species presence data for upper 100m depth
hspen <- fread(paste(datadir, 'Aquamaps/2020/raw_data/hspen.csv', sep=''), header=TRUE)
names(hspen) <- tolower(names(hspen))  # Convert column names to lowercase for consistency

## Load species occurrence summary data
tdat <- read.csv(paste(datadir, 'Aquamaps/2020/raw_data/speciesoccursum.csv', sep=''), header=TRUE)
names(tdat) <- tolower(names(tdat))  # Convert column names to lowercase

## Remove bird species (Class Aves)
tdat <- subset(tdat, !(class %in% c('Aves')))

## Filter hspen dataset: 
## - Species with minimum depth ≤ 100m and maximum depth ≤ 1000m
## - Species present in tdat dataset
hspen <- subset(hspen, depthmin <= 100 & depthmax <= 1000 & speciesid %in% unique(tdat$speciesid))

## Keep only species ID column
spinclude <- subset(hspen, select=c('speciesid'))

## Remove hspen object to free memory
rm(hspen)

########################################################################

## FUNCTIONS THAT TRANSFORM AND STANDARDIZE RAW INDICES

########################################################################

## HABITAT FRAGMENTATION
tf.AC.hfrag <- function(x) {
  ## Converts habitat fragmentation index using an exponential decay function
  gam <- 8  # Steepness parameter
  xout <- exp(-x * gam)  # Transformation
  print(mean(xout, na.rm=TRUE))  # Print mean of transformed values
  return(xout)
}

## RANGE AREA
tf.AC.rarea <- function(x) {
  ## Converts range area index using an exponential decay function
  gam <- 20  # Steepness parameter
  xout <- 1 - exp(-x * gam)  
  print(mean(xout, na.rm=TRUE))  
  return(xout)
}

## LATITUDE RANGE
tf.AC.lrange <- function(x) {
  ## Normalizes latitude range by dividing by 90 (max latitude)
  xout <- x / 90
  print(mean(xout, na.rm=TRUE))  
  return(xout)
}

## HABITAT RANGE
tf.AC.hrange <- function(x1, x2) {
  ## Computes the average of two habitat range values
  xout <- (x1 + x2) / 2
  print(mean(xout, na.rm=TRUE))  
  return(xout)
}

## BODY LENGTH TRANSFORMATION
tf.AC.lmax <- function(x) {
  lmax <- 3000  ## Max length (approx. 3300 for a blue whale)
  xout <- log10(x + 1) / log10(lmax + 1)  # Normalize on a log scale
  xout <- 1 - xout  # Invert scaling
  print(summary(xout))  
  return(xout)
}

## AC TAUC TRANSFORMATION
tf.AC.tauc <- function(x) {
  gam1 <- 4  # Exponential growth factor
  xout <- exp(x * gam1) / exp(gam1)  # Normalize with exponential function
  print(summary(xout))  
  return(xout)
}

## AC TRNG TRANSFORMATION
tf.AC.trng <- function(x) {
  ## Normalizes temperature range (TRNG) by dividing by 30
  xout <- x / 30
  xout <- ifelse(xout > 1, 1, xout)  # Cap values at 1
  print(summary(xout))  
  return(xout)
}

## TEMPERATURE VARIATION
tf.AC.tvar <- function(x1, x2) {
  ## Computes the average of two temperature variation values
  xout <- (x1 + x2) / 2
  print(summary(xout))  
  return(xout)
}

## HUMAN INFLUENCE INDEX TRANSFORMATION
tf.S.HII <- function(x) {
  hiimax <- 12  ## Max possible value
  gam <- 0.6
  xout <- (1 - exp(-x * gam))  # Apply exponential decay function
  print(summary(xout))  
  return(xout)
}

## VERTICAL RANGE TRANSFORMATION
tf.S.vertrange <- function(x) {
  gam <- 0.007  # Exponential decay factor
  xout <- exp(-x * gam)
  print(summary(xout))  
  return(xout)
}

## MAXIMUM DEPTH TRANSFORMATION
tf.S.depthmax <- function(x) {
  gam <- 0.007
  xout <- exp(-x * gam)
  print(summary(xout))  
  return(xout)
}

## COMBINED VERTICAL INDEX
tf.S.vind <- function(x1, x2) {
  ## Computes the average of vertical range and max depth
  xout <- (x1 + x2) / 2
  print(summary(xout))  
  return(xout)
}

## TEMPERATURE SEASONALITY INDEX
tf.S.TSMr <- function(x) {
  gam <- 0.33
  xout <- exp(-x * gam)
  print(summary(xout))  
  return(xout)
}

## ENVIRONMENTAL VELOCITY TRANSFORMATION
tf.E.vel <- function(x) {
  gam <- 0.02
  x <- ifelse(x < 0, 0, x)  # Ensure non-negative values
  xout <- 1 - exp(-x * gam)  # Apply exponential decay function
  print(summary(xout))  
  return(xout)
}

## HABITAT LOSS TRANSFORMATION
tf.E.plost <- function(x) {
  gam <- 5
  xout <- 1 - exp(-x * gam)
  print(summary(xout))  
  return(xout)
}

## TEMPERATURE OVERLAP INDEX
tf.E.toe <- function(x) {
  gam <- 30
  xout <- exp(-x / gam)
  print(summary(xout))  
  return(xout)
}

## ENVIRONMENTAL CHANGE INDEX
tf.E.nrchng <- function(x) {
  gam <- 4
  xout <- 1 - exp(-x * gam)
  print(summary(xout))  
  return(xout)
}

########################################################

## DISCOUNT RATE FUNCTION

########################################################
discrate <- function(nyrs, nmods, theta, epsilon) {
  ## Computes a discount rate using a combination of years, model count, and parameters
  return(round(((nyrs / (100 * theta)) + ((nmods / (-20 * theta)) + epsilon)), digits=3))
}

print('GPARAMS LOADED')  ## Confirmation message





