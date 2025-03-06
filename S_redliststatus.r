################################################################################################

## R CODE FOR: BOYCE, D. ET AL. 2021. A CLIMATE RISK INDEX FOR MARINE LIFE
## CODE AND ANALYSES DEVELOPED BY DANIEL BOYCE (DBOYCE@DAL.CA); PLEASE REFERENCE USE ACCORDINGLY
## IMPLEMENTED IN R VERSION 4.04 AND RUN ON COMPUTE CANADA COMPUTING CLUSTER
## LAST MODIFIED: AUGUST 2021

################################################################################################

## LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS

# Source an external R script containing global parameters and settings
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

# Define a function to set the code directory dynamically based on the current working directory
cdir <- function(){
    codedir <- ifelse(grepl('sailfish', getwd(), fixed = TRUE),
                      'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                      'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
    
    # Source the parameter file from the determined directory
    source(paste(codedir, 'GPARAMs_CC_Vuln.r', sep=''))
}
## cdir() # This function is currently not executed

#############################################################

## RED LIST OF SPECIES (SPECIES-LEVEL)

#############################################################

# Read in the IUCN Red List species assessment data
rl <- read.csv(paste(datadir, 'redlist/marine_species_data_2022/assessments.csv', sep=''), 
               header = TRUE, encoding = "UTF-8")

# Remove unnecessary columns to streamline the dataset
rl <- rl[, !(names(rl) %in% c('yearPublished', 'assessmentDate', 'criteriaVersion', 'language',
                              'rationale', 'habitat', 'population', 'range', 'yearLastSeen',
                              'possiblyExtinct', 'possiblyExtinctInTheWild', 'threats',
                              'conservationActions', 'redlistCriteria'))]

# Convert column names to lowercase for consistency
names(rl) <- tolower(names(rl))

# Standardize scientific name format (lowercase, remove extra spaces, and trim whitespace)
rl$scientificname <- tolower(as.character(rl$scientificname))
rl$scientificname <- gsub('  ', ' ', rl$scientificname)  # Replace double spaces with single space
rl$scientificname <- trimws(rl$scientificname, which = 'both')  # Trim leading and trailing spaces

# Convert relevant columns to character type
rl$redlistcategory <- as.character(rl$redlistcategory)
rl$populationtrend <- as.character(rl$populationtrend)
rl$systems <- as.character(rl$systems)
rl$realm <- as.character(rl$realm)
rl$scopes <- as.character(rl$scopes)

## FORMAT STATUS

# Remove additional unnecessary columns
rl <- rl[, !(names(rl) %in% c('usetrade', 'realm', 'scopes', 'systems'))]

# Exclude species classified as Extinct or Extinct in the Wild
rl <- subset(rl, !(redlistcategory %in% c('Extinct', 'Extinct in the Wild')))

# Assign numerical risk scores to Red List categories
rl$s.rlstatus <- ifelse(rl$redlistcategory %in% c('Critically Endangered'), 1, NA)
rl$s.rlstatus <- ifelse(rl$redlistcategory %in% c('Endangered'), 0.75, rl$s.rlstatus)
rl$s.rlstatus <- ifelse(rl$redlistcategory %in% c('Vulnerable'), 0.5, rl$s.rlstatus)
rl$s.rlstatus <- ifelse(rl$redlistcategory %in% c('Near Threatened', 'Lower Risk/near threatened'), 0.25, rl$s.rlstatus)
rl$s.rlstatus <- ifelse(rl$redlistcategory %in% c('Least Concern', 'Lower Risk/least concern', 'Lower Risk/conservation dependent'), 0, rl$s.rlstatus)

## FORMAT POPULATION TREND

# Assign numerical trend scores to population trends
rl$s.rltrend <- ifelse(rl$populationtrend %in% c('Decreasing'), 1, NA)
rl$s.rltrend <- ifelse(rl$populationtrend %in% c('Stable'), 0.5, rl$s.rltrend)
rl$s.rltrend <- ifelse(rl$populationtrend %in% c('Increasing'), 0, rl$s.rltrend)

# Save the processed dataset as an RData file for future use
save(rl, file = paste(datadirout, 'S_redlist_1deg_unimputed.RData', sep=''))





##########################################################
## IMPUTE MISSING RL STATUS
##########################################################

######################################################
## LOAD ENVIRONMENTAL ENVELOPE DATA
hspen <- fread(paste(datadir, 'Aquamaps/2020/raw_data/hspen.csv', sep=''), header=TRUE)
hspen <- data.frame(hspen)
names(hspen) <- tolower(names(hspen))  # Convert column names to lowercase

# Calculate environmental ranges
hspen$latrange <- abs(hspen$nmostlat - hspen$smostlat)
hspen$depthrange <- hspen$depthmax - hspen$depthmin

# Subset relevant columns
hspen <- subset(hspen, select=c('speciesid', 'depthmax', 'depthmin', 'tempmin', 'tempmax', 
                                'salinitymin', 'salinitymax', 'primprodmin', 'primprodmax', 
                                'iceconmax', 'oxymin', 'oxymax', 'landdistmax', 'landdistmin'))
print(length(unique(hspen$speciesid)))  # Print number of unique species

## NORMALIZE TO IMPROVE IMPUTATION
# Apply log transformation to selected variables to normalize distributions
hspen$depthmax <- log10(hspen$depthmax + 1)
hspen$primprodmin <- log10(hspen$primprodmin + 1)
hspen$primprodmax <- log10(hspen$primprodmax + 1)
hspen$iceconmax <- log10(hspen$iceconmax + 1)
hspen$landdistmax <- log10(hspen$landdistmax + 1)
hspen$landdistmin <- log10(hspen$landdistmin + 1)

########################################################
## LOAD GLOBAL AREA DATA FOR EACH SPECIES
load(paste(datadirout, 'AC_hrange_1deg.RData', sep=''))
hrange <- subset(hrange, select=c('speciesid', 'areakm2'))
hrange$areakm2 <- log10(hrange$areakm2 + 1)  # Normalize area size
print(summary(hrange))
print(length(unique(hrange$speciesid)))

########################################################
## LOAD GLOBAL HABITAT FRAGMENTATION DATA
load(paste(datadirout, 'AC_hfrag_1deg.RData', sep=''))
hfrag <- subset(hfrag, select=c('speciesid', 'hfrag'))
hfrag$hfrag <- log10(hfrag$hfrag + 1)  # Normalize habitat fragmentation measure
print(summary(hfrag))
print(length(unique(hfrag$speciesid)))

########################################################
## LOAD TAXONOMIC INFORMATION
tdat <- read.csv(paste(datadir, 'Aquamaps/2020/raw_data/speciesoccursum.csv', sep=''), header=TRUE)
names(tdat) <- tolower(names(tdat))  # Convert column names to lowercase
tdat <- subset(tdat, select=c('speciesid', 'kingdom', 'phylum', 'order', 'family', 'genus'))
print(head(tdat))

#######################################################
## LOAD WORMS TAXONOMIC BACKBONE
load(paste(datadirout, 'WoRMS_mlist.RData', sep=''))
mlist$cname <- ifelse(is.na(mlist$cname.am), mlist$cname.rl, mlist$cname.am)

# Remove unnecessary taxonomic columns
mlist <- mlist[, !(names(mlist) %in% c('taxonrank', 'internaltaxonid.rl', 'spec.rl', 'speccode.am', 
                                       'spec.am', 'sdist.am', 'sdist.rl', 'cname.rl', 'cname.am', 
                                       'genus', 'acceptedname', 'kingdom', 'phylum', 'order', 'family', 'cname'))]

# Remove instances where there is no RL or AM species data
mlist <- unique(subset(mlist, !is.na(speciesid.am) & !is.na(assessmentid.rl)))

# Identify cases where multiple Red List (RL) species correspond to a single Aquamaps (AM) species
dt <- data.frame(speciesid.am=sort(unique(mlist$speciesid.am)),
                 n=tapply(mlist$assessmentid.rl, mlist$speciesid.am, function(x) length(unique(x))))
dt <- dt[order(dt$n, decreasing=TRUE),]
a <- subset(dt, n > 1)  # Find species with multiple RL assessments

# For species with multiple RL assessments, set RL assessment ID to NA
mlist$assessmentid.rl <- ifelse(mlist$speciesid.am %in% a$speciesid.am, NA, mlist$assessmentid.rl)
mlist <- subset(mlist, !is.na(assessmentid.rl))
names(mlist) <- c('assessmentid', 'speciesid')

#######################################################
## COMBINE ALL DATA LAYERS
dat <- left_join(hspen, hrange, by='speciesid')
dat <- left_join(dat, hfrag, by='speciesid')
dat <- left_join(dat, tdat, by='speciesid')
dat <- left_join(dat, mlist, by='speciesid')
dat <- left_join(dat, subset(rl, select=c('assessmentid', 's.rlstatus')), by='assessmentid')
print(length(unique(dat$speciesid)))

## ENSURE RL STATUS IS CATEGORICAL
dat$s.rlstatus <- factor(dat$s.rlstatus, levels=as.character(sort(unique(dat$s.rlstatus))))
dat <- dat[, !(names(dat) %in% c('assessmentid'))]  # Remove assessment ID

# Convert taxonomic info into numerical factors
dat$kingdom <- as.numeric(as.factor(dat$kingdom))
dat$phylum <- as.numeric(as.factor(dat$phylum))
dat$order <- as.numeric(as.factor(dat$order))
dat$family <- as.numeric(as.factor(dat$family))
dat$genus <- as.numeric(as.factor(dat$genus))
print(dim(dat))
print(summary(dat))

#######################################################
## PERFORM MISSING DATA IMPUTATION USING MICE PACKAGE
system.time(tempdata <- mice(dat[, -1], m=5, maxit=50, meth='sample', seed=500, print=FALSE))

## GET COMPLETE IMPUTED DATASET
impdat <- complete(tempdata)

## CALCULATE AVERAGE RL STATUS ACROSS 5 IMPUTATIONS
dt <- data.frame(id1=complete(tempdata, 1)$s.rlstatus,
                 id2=complete(tempdata, 2)$s.rlstatus,
                 id3=complete(tempdata, 3)$s.rlstatus,
                 id4=complete(tempdata, 4)$s.rlstatus,
                 id5=complete(tempdata, 5)$s.rlstatus)

f <- function(d) {
    d2 <- data.frame(s.rlstatus=sort(d[1,])[3])  # Select median RL status
    names(d2) <- c('s.rlstatus')
    return(d2)
}

impdatav <- adply(dt, 1, .fun=f)
impdatav$speciesid <- dat$speciesid
rl <- subset(impdatav, select=c('speciesid', 's.rlstatus'))
rl$s.rlstatus <- as.numeric(as.character(rl$s.rlstatus))
print(length(unique(rl$speciesid)))

# SAVE FINAL RED LIST DATA
save(rl, file=paste(datadirout, 'S_redlist_1deg.RData', sep=''))
print('Redlist imputation done')
print(summary(rl))

# CLEAN UP WORKSPACE
rm(list=ls(all=TRUE))
gc()  # Run garbage collection


