## THIS CODE WILL: CREATE A SPECIES 'MASTER LIST' CONTAINING ALL WORMS NAMES AND THOSE THAT CONTAIN FUZZY MATCHES 
## IN THE IUCN RED LIST, AND AQUAMAPS. THIS MASTERLIST CAN THEN BE USED TO COMBINE THE DIFFERENT DATA SOURCES (REDLIST, AQUAMAPS, ETC.) BY SPECIES

## LOAD REQUIRED PACKAGES, PATHNAMES, PARAMETERS
source('/home/dgboyce/projects/def-bworm/dgboyce/CC_vulnerability/code/GPARAMs_CC_Vuln.r')

## FUNCTION TO SET DIRECTORY PATH BASED ON THE WORKING DIRECTORY
cdir<-function(){
    codedir<-ifelse(grepl('sailfish', getwd(), fixed = TRUE),
                    'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/',
                    'C:/Users/danie/Documents/aalldocuments/literature/research/active/CC_vulnerability/code/')
    source(paste(codedir, 'GPARAMs_CC_Vuln.r', sep=''))
}
## cdir() ## (Commented out: manually sourcing parameter file instead)

print('WoRMS_species_masterlist_start')

## LOAD PREVIOUSLY FORMATTED WORMS TAXONOMY DATA
load(paste(datadir, 'worms_taxonomy/WoRMS_formatted.RData', sep=''))

## SELECT RELEVANT COLUMNS FROM WORMS DATAFRAME
worms <- subset(worms, select = c('acceptedname', 'synonym', 'kingdom', 'phylum', 'order', 'family', 'genus', 'taxonrank', 'taxonomicstatus'))

## PRINT CONFIRMATION MESSAGE
print('worms read done')

#####################################################
# PROCESS AQUAMAPS DATA
#####################################################

## LOAD AQUAMAPS SPECIES OCCURRENCE DATA
am <- read.csv(paste(datadir, 'Aquamaps/2020/raw_data/speciesoccursum.csv', sep=''), header=TRUE)
names(am) <- tolower(names(am))

## EXTRACT PHYLA FROM AQUAMAPS DATA
aphylum <- subset(am, select=c('phylum'))
aphylum$phylum <- as.character(aphylum$phylum)

## LOAD ENVIRONMENTAL ENVELOPE DATA FOR AQUAMAPS
hspen <- fread(paste(datadir, 'Aquamaps/2020/raw_data/hspen.csv', sep=''), header=TRUE)
names(hspen) <- tolower(names(hspen))

## FILTER AQUAMAPS DATA TO SPECIES FOUND IN THE UPPER 100M DEPTH
am <- subset(am, speciesid %in% spinclude$speciesid)
gc()

## STANDARDIZE SPECIES NAMES IN AQUAMAPS
am$spec.am <- tolower(paste(am$genus, am$species))
am$spec.am <- trimws(am$spec.am, which='both')

## SELECT RELEVANT COLUMNS FROM AQUAMAPS DATASET
am <- subset(am, select = c('speciesid', 'speccode', 'fbname', 'spec.am'))
names(am) <- c('speciesid.am', 'speccode.am', 'cname.am', 'spec.am')
am$speciesid.am <- as.character(am$speciesid.am)
am$cname.am <- as.character(am$cname.am)
am$spec.am <- as.character(am$spec.am)

## RETAIN ONLY WORMS ENTRIES THAT ARE ALSO PRESENT IN AQUAMAPS (FOR VULNERABILITY ANALYSIS)
## ALSO INCLUDE SPECIES WITH UNIDENTIFIED PHYLA ('')
worms <- subset(worms, phylum %in% c('', aphylum$phylum))

########################################################
# PROCESS IUCN RED LIST DATA
########################################################

## LOAD IUCN RED LIST SPECIES ASSESSMENTS
rl <- read.csv(paste(datadir, 'redlist/marine_species_data_2022/assessments.csv', sep=''), header=TRUE, encoding="UTF-8")

## REMOVE UNNECESSARY COLUMNS FROM RED LIST DATA
rl <- rl[, !(names(rl) %in% c('yearPublished', 'assessmentDate', 'criteriaVersion', 'language', 'rationale', 'habitat',
                              'population', 'range', 'yearLastSeen', 'possiblyExtinct', 'possiblyExtinctInTheWild',
                              'threats', 'conservationActions', 'redlistCriteria', 'redlistCategory', 'populationTrend',
                              'useTrade', 'systems', 'realm', 'scopes'))]

## FORMAT SPECIES NAMES IN RED LIST DATA
names(rl) <- tolower(names(rl))
rl$scientificname <- tolower(as.character(rl$scientificname))
rl$scientificname <- gsub('  ', ' ', rl$scientificname)
rl$scientificname <- trimws(rl$scientificname, which='both')

#########################################################
# ADD COMMON NAMES TO RED LIST DATA
#########################################################

## LOAD RED LIST COMMON NAMES DATASET
cn <- read.csv(paste(datadir, 'redlist/marine_species_data/common_names.csv', sep=''), header=TRUE, encoding="UTF-8")
names(cn) <- tolower(names(cn))

## FILTER TO INCLUDE ONLY ENGLISH COMMON NAMES MARKED AS MAIN
cn <- subset(cn, language == 'English' & main == 'true')
cn$name <- as.character(cn$name)

## FUNCTION TO ADD UNIQUE IDENTIFIERS TO COMMON NAMES DATAFRAME
f <- function(d){
    d$id <- seq(1, dim(d)[1], 1)
    return(d)
}

## APPLY FUNCTION TO GROUP DATA BY INTERNAL TAXON ID
l <- dlply(cn, .(internaltaxonid), .fun = f, .progress = 'text')
cn2 <- rbindlist(l)

## RETAIN ONLY THE FIRST ENTRY FOR EACH SPECIES
cn <- subset(cn2, id == 1)
cn <- subset(cn, select = c('internaltaxonid', 'name'))

## MERGE COMMON NAMES WITH RED LIST DATA
rl <- merge(rl, cn, by = c('internaltaxonid'), all.x = TRUE, all.y = FALSE)
names(rl) <- c('internaltaxonid.rl', 'assessmentid.rl', 'spec.rl', 'cname.rl')

## CLEAN UP UNUSED OBJECTS TO FREE MEMORY
rm(list = c('hspen', 'cn', 'cn2'))
gc()




######################################################
# JOIN WORMS SYNONYMS TO REDLIST USING FUZZY MATCHING
# FUNCTION PERFORMS THE FUZZY MATCHING
######################################################

# Define a function to perform fuzzy matching between WoRMS synonyms and Red List species
mfun <- function(d) {
    rl$sdist.rl <- stringdist(d$synonym, rl$spec.rl, method = 'osa')  # Compute string distance
    return(cbind(d, subset(rl, sdist.rl == min(rl$sdist.rl, na.rm = TRUE))[1,]))  # Join closest match
}

# Load parallel computing libraries
library(doSNOW)
library(foreach)
library(doParallel)

# Detect and use all available cores for parallel processing
UseCores <- detectCores()
print(UseCores)
cl <- makeCluster(UseCores, type = "SOCK")  # Create cluster
registerDoSNOW(cl)  # Register cluster
## stopCluster(cl)  # Uncomment to stop cluster manually if needed

# Perform fuzzy matching in parallel and measure execution time
system.time(
    mlist <- adply(
        worms, 1, .fun = mfun, .progress = 'text', .parallel = TRUE,
        .paropts = list(.packages = c('stringdist'), .export = c('worms', 'rl'))
    )
)

# SET INSTANCES WHERE STRING DISTANCE >1 TO NA TO FILTER OUT BAD MATCHES
mlist$internaltaxonid.rl <- ifelse(mlist$sdist.rl > 1, NA, mlist$internaltaxonid.rl)
mlist$assessmentid.rl <- ifelse(mlist$sdist.rl > 1, NA, mlist$assessmentid.rl)
mlist$spec.rl <- ifelse(mlist$sdist.rl > 1, NA, mlist$spec.rl)
mlist$cname.rl <- ifelse(mlist$sdist.rl > 1, NA, mlist$cname.rl)

######################################################
# JOIN WORMS SYNONYMS TO AQUAMAPS USING FUZZY MATCHING
######################################################

# Define function for fuzzy matching with AquaMaps species
mfun <- function(d) {
    am$sdist.am <- stringdist(d$synonym, am$spec.am, method = 'osa')  # Compute string distance
    return(cbind(d, subset(am, sdist.am == min(am$sdist.am, na.rm = TRUE), 
                select = c('speciesid.am', 'speccode.am', 'cname.am', 'spec.am', 'sdist.am'))[1,]))  # Join closest match
}

# Perform fuzzy matching with AquaMaps in parallel and measure execution time
system.time(
    mlist2 <- adply(
        mlist, 1, .fun = mfun, .progress = 'text', .parallel = TRUE,
        .paropts = list(.packages = c('stringdist'), .export = c('mlist', 'am'))
    )
)

# SET ALL INSTANCES WITH POOR FUZZY MATCHES (DISTANCE >1) TO NA
mlist2$speciesid.am <- ifelse(mlist2$sdist.am > 1, NA, mlist2$speciesid.am)
mlist2$speccode.am <- ifelse(mlist2$sdist.am > 1, NA, mlist2$speccode.am)
mlist2$spec.am <- ifelse(mlist2$sdist.am > 1, NA, mlist2$spec.am)
mlist2$cname.am <- ifelse(mlist2$sdist.am > 1, NA, mlist2$cname.am)

########################################################
# FILTER MATCHES FOR FINAL SPECIES MASTER LIST
#########################################################

# Subset relevant columns for final analysis
mlist <- subset(mlist2, select = c(
    'acceptedname', 'kingdom', 'phylum', 'order', 'family', 'genus', 'taxonrank',
    'assessmentid.rl', 'internaltaxonid.rl', 'spec.rl', 'cname.rl', 'sdist.rl',
    'speciesid.am', 'speccode.am', 'cname.am', 'spec.am', 'sdist.am'
))

# Function to select the best matching records per species
f <- function(d) {
    d2 <- subset(d, select = c('kingdom', "phylum", "order", "family", "genus", "taxonrank"))[1,]
    drl <- subset(d, sdist.rl == min(d$sdist.rl, na.rm = TRUE), select = c('sdist.rl', 'assessmentid.rl', 'internaltaxonid.rl', 'spec.rl', 'cname.rl'))
    dam <- subset(d, sdist.am == min(d$sdist.am, na.rm = TRUE), select = c('sdist.am', 'speciesid.am', 'speccode.am', 'cname.am', 'spec.am'))
    
    # Assign the best Red List match
    d2$assessmentid.rl <- unique(na.omit(drl$assessmentid.rl))[1]
    d2$internaltaxonid.rl <- unique(na.omit(drl$internaltaxonid.rl))[1]
    d2$spec.rl <- unique(na.omit(drl$spec.rl))[1]
    d2$cname.rl <- unique(na.omit(drl$cname.rl))[1]
    
    # Assign the best AquaMaps match
    d2$speciesid.am <- unique(na.omit(dam$speciesid.am))[1]
    d2$speccode.am <- unique(na.omit(dam$speccode.am))[1]
    d2$cname.am <- unique(na.omit(dam$cname.am))[1]
    d2$spec.am <- unique(na.omit(dam$spec.am))[1]
    
    # Store minimum distances
    d2$sdist.am <- min(dam$sdist.am, na.rm = TRUE)
    d2$sdist.rl <- min(drl$sdist.rl, na.rm = TRUE)
    return(d2)
}

# Remove blank accepted names
mlist <- subset(mlist, !(acceptedname %in% c('')))

# Aggregate final matched data by accepted species names
system.time(
    mlist <- ddply(
        mlist, .(acceptedname), .fun = f, .progress = 'text', .parallel = TRUE,
        .paropts = list(.export = c('mlist'))
    )
)

# Save the final master list to an RData file
save(mlist, file = paste(datadirout, 'WoRMS_mlist.RData', sep = ''))

# Clean up environment and free memory
rm(list = ls(all = TRUE))
gc()

# Print completion message
print('WoRMS_species_masterlist_done')
