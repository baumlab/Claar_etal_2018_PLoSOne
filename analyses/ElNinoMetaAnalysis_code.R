# Global patterns and impacts of El Niño events on coral reefs: a meta-analysis

# Authors: Danielle C. Claar $^1$, Lisa Szostek$^1$, Jamie M. McDevitt-Irwin $^1$, Julian J.Schanze $^2$, Julia K. Baum $^1$
# Institute: $^1$ Department of Biology, University of Victoria, PO BOX 1700 Station CSC, Victoria, British Columbia, V8W 2Y2, Canada  
# $^2$ Earth and Space Research, 2101 Fourth Ave., Suite 1310, Seattle WA 98121-2350 USA
# Corresponding Author: Danielle C. Claar, Tel: (208) 250-0161, Email: dclaar@uvic.ca



## Install and source necessary packages
rm(list = ls())

## Load and attach necessary add-on packages
library(metafor)
library(reshape)
library(glmulti)

## You need to set your working directory to be within the GitHub repo "ElNino_MetaAnalysis" to proceed with the following filepaths.

################################################
# Import full data set and configure factors
ElNino_FULL <- read.csv("data/ElNino_FULL_dhw.csv") # Import full data set
ElNino_afteronly <- read.csv("data/ElNino_afteronly_dhw.csv") # Import "after only" bleaching data set. These are data sets which have no before-El Niño/La Niña metric of bleaching, but provide all necessary information to evaluate coral bleaching either during or after the El Niño/La Niña event.

ElNino_FULL$ACC<-as.factor(ElNino_FULL$ACC) # Set ACC as a factor, to account for differences among studies
ElNino_afteronly$ACC<-as.factor(ElNino_afteronly$ACC) # Set ACC as a factor, to account for differences among studies

ElNino_FULL$El.Nino.simple<-as.factor(ElNino_FULL$El.Nino.simple) # Set El Niño year as a factor
## This is saved in the csv as the second year in the El Niño, e.g. for 97/98 it is 1998. (This is because the csv kept messing up the formatting of any other name/date structure)
ElNino_afteronly$El.Nino.simple<-as.factor(ElNino_afteronly$El.Nino.simple) # Set El Niño year as a factor

names(ElNino_FULL)[names(ElNino_FULL)=="V8"] <- "DHWnow" # Current DHW calculated at the queried time
names(ElNino_FULL)[names(ElNino_FULL)=="V9"] <- "MaxDHW_beforeafter" # Maximum DHW during the current year and the past 2 years. This value may be in the future (refer to TimeLag_beforeafter value, if negative, this is in the future)
names(ElNino_FULL)[names(ElNino_FULL)=="V10"] <- "TimeLag_beforeafter" # Time (in days) between queried time and maximum DHW, can be positive (happened before the queried date) or negative (happened after the queried date)
names(ElNino_FULL)[names(ElNino_FULL)=="V11"] <- "MaxDHW" # Same as MaxDHW_beforeafter, except it only considers values that occured before the queried date (i.e. shows only what the corals at that location have experienced so far)
names(ElNino_FULL)[names(ElNino_FULL)=="V12"] <- "TimeLag" # Same as TimeLag_beforeafter, and paired with MaxDHW. Shows only the time since the maxDHW the corals have experienced at that location so far
names(ElNino_FULL)[names(ElNino_FULL)=="V13"] <- "SSTnow" # Current temperature at the queried time and location
names(ElNino_FULL)[names(ElNino_FULL)=="V14"] <- "SSTmean" # Mean temperature between 1982 and 2014 (AVHRR Reynolds OI SST)
names(ElNino_FULL)[names(ElNino_FULL)=="V15"] <- "SSTvar" # Temperature variance between 1982 and 2014
names(ElNino_FULL)[names(ElNino_FULL)=="V16"] <- "SSTstd" # Temperature standard deviation between 1982 and 2014
names(ElNino_FULL)[names(ElNino_FULL)=="V17"] <-  "MMM" # Maximum monthly mean for the queried location
names(ElNino_FULL)[names(ElNino_FULL)=="V18"] <- "MMMind" # Month in which the MMM occurs, ranges from January=1 to December=12.

names(ElNino_afteronly)[names(ElNino_afteronly)=="V8"] <- "DHWnow" # Current DHW calculated at the queried time
names(ElNino_afteronly)[names(ElNino_afteronly)=="V9"] <- "MaxDHW_beforeafter" # Maximum DHW during the current year and the past 2 years. This value may be in the future (refer to TimeLag_beforeafter value, if negative, this is in the future)
names(ElNino_afteronly)[names(ElNino_afteronly)=="V10"] <- "TimeLag_beforeafter" # Time (in days) between queried time and maximum DHW, can be positive (happened before the queried date) or negative (happened after the queried date)
names(ElNino_afteronly)[names(ElNino_afteronly)=="V11"] <- "MaxDHW" # Same as MaxDHW_beforeafter, except it only considers values that occured before the queried date (i.e. shows only what the corals at that location have experienced so far)
names(ElNino_afteronly)[names(ElNino_afteronly)=="V12"] <- "TimeLag" # Same as TimeLag_beforeafter, and paired with MaxDHW. Shows only the time since the maxDHW the corals have experienced at that location so far
names(ElNino_afteronly)[names(ElNino_afteronly)=="V13"] <- "SSTnow" # Current temperature at the queried time and location
names(ElNino_afteronly)[names(ElNino_afteronly)=="V14"] <- "SSTmean" # Mean temperature between 1982 and 2014 (AVHRR Reynolds OI SST)
names(ElNino_afteronly)[names(ElNino_afteronly)=="V15"] <- "SSTvar" # Temperature variance between 1982 and 2014
names(ElNino_afteronly)[names(ElNino_afteronly)=="V16"] <- "SSTstd" # Temperature standard deviation between 1982 and 2014
names(ElNino_afteronly)[names(ElNino_afteronly)=="V17"] <-  "MMM" # Maximum monthly mean for the queried location
names(ElNino_afteronly)[names(ElNino_afteronly)=="V18"] <- "MMMind" # Month in which the MMM occurs, ranges from January=1 to December=12.

################################################

# For manuscript section: El Niño/La Niña Effects on Coral Bleaching and Coral Cover
# Calculate absolute and proportional changes in coral cover
change <- cbind(ElNino_FULL$ACC, ElNino_FULL$Mean_Cover_Before, ElNino_FULL$Mean_Cover_Resp, ElNino_FULL$Mean_Cover_Before - ElNino_FULL$Mean_Cover_Resp) # Create a new dataframe to calculate absolute change in cover
change <- change[which(!is.na(change[,2])), ] # Remove NAs in the cover column, specifically catches and removes any bleaching data
abs_max_change <- max(change[,4]) # Calculate the maximum absolute loss of coral cover (e.g if coral cover drops from 50% to 10%, this is 40%)
p_change <- change[,4]/change[,2]*100 # Calculate the proportional loss of coral cover (e.g. if coral cover drops from 50% to 10%, this is 80%)
change <- cbind(change, p_change) # add proportional change to dataframe
p_max_change <- max(change[,5]) # Calculate the maximum proportional loss of coral cover

max_bleach <- max(ElNino_afteronly$Mean_Bleaching_Resp) # Calculate the maximum coral bleaching

################################################

# Add small number (0.1) to standard deviations that == 0
# Because some papers reported mean/error on graphs, but the error was so small it was not visible. Note: this does not add a small number to studies that did not report error (these studies were not included)
ElNino_FULL$SD_Cover_Before[ElNino_FULL$SD_Cover_Before == "0"] <- "0.1" # If SD was 0, change it to 0.1
ElNino_FULL$SD_Cover_Before <- as.numeric(ElNino_FULL$SD_Cover_Before) # Make sure it's numeric
# str(ElNino_FULL$SD_Cover_Before) #Check structure to make sure it's right
ElNino_FULL$SD_Bleaching_Before[ElNino_FULL$SD_Bleaching_Before == "0"] <- "0.1" # If SD was 0, change it to 0.1
ElNino_FULL$SD_Bleaching_Before <- as.numeric(ElNino_FULL$SD_Bleaching_Before) # Make sure it's numeric
# str(ElNino_FULL$SD_Bleaching_Before) # Check strucutre to make sure it's right

ElNino_FULL$SD_Cover_Resp[ElNino_FULL$SD_Cover_Resp == "0"] <- "0.1" # If SD was 0, change it to 0.1
ElNino_FULL$SD_Cover_Resp <- as.numeric(ElNino_FULL$SD_Cover_Resp) # Make sure it's numeric
# str(ElNino_FULL$SD_Cover_Resp) # Check structure to make sure it's right
ElNino_FULL$SD_Bleaching_Resp[ElNino_FULL$SD_Bleaching_Resp == "0"] <- "0.1" # If SD was 0, change it to 0.1
ElNino_FULL$SD_Bleaching_Resp <- as.numeric(ElNino_FULL$SD_Bleaching_Resp) # Make sure it's numeric
# str(ElNino_FULL$SD_Bleaching_Resp) # Check structure to make sure it's right

ElNino_afteronly$SD_Bleaching_Resp[ElNino_afteronly$SD_Bleaching_Resp == "0"] <- "0.1" # If SD was 0, change it to 0.1
ElNino_afteronly$SD_Bleaching_Resp <- as.numeric(ElNino_afteronly$SD_Bleaching_Resp) # Make sure it's numeric

################################################

# Subset data set into bleaching/cover and check that studies without standard deviation (SD) are removed
cover<-ElNino_FULL[which(ElNino_FULL$Parameter=="Cover"), ] # Subset full data set to separate coral cover
cover<-cover[which(!is.na(cover$SD_Cover_Before)), ] # Remove any coral cover datapoints which do not have SD in "Before"
cover<-cover[which(!is.na(cover$SD_Cover_Resp)), ] # Remove any coral cover datapoints which do not have SD in "Response"

bleaching<-ElNino_FULL[which(ElNino_FULL$Parameter=="Bleaching"), ] # Subset full data set to separate coral bleaching
bleaching<-bleaching[which(!is.na(bleaching$SD_Bleaching_Resp)), ] # Remove any coral bleaching data points which do not have SD in "Response"

################################################

# Subset/Combine all After-only bleaching data, simulate "before" mean & error
b <- ElNino_afteronly[which(!is.na(ElNino_afteronly$SD_Bleaching_Resp)), ] # Create temp variables to merge bleaching afteronly and only the after values of whole bleaching data set
temp.b <- bleaching # Create temp variables to merge bleaching afteronly and only the after values of whole bleaching data set

rbind.match.columns <- function(input1, input2) { # Thanks Amy Whitehead for sharing your code online! https://amywhiteheadresearch.wordpress.com/2013/05/13/combining-dataframes-when-the-columns-dont-match/
  n.input1 <- ncol(input1)
  n.input2 <- ncol(input2)
  
  if (n.input2 < n.input1) {
    TF.names <- which(names(input2) %in% names(input1))
    column.names <- names(input2[, TF.names])
  } else {
    TF.names <- which(names(input1) %in% names(input2))
    column.names <- names(input1[, TF.names])
  }
  
  return(rbind(input1[, column.names], input2[, column.names]))
}

ElNino_afteronly<-rbind.match.columns(b, temp.b) # Use function to merge b and temp.b
temp.df <- data.frame(Mean_Bleaching_Before_sim = NA, SD_Bleaching_Before_sim = NA, Mean_Bleaching_Before_sim.plus = NA, SD_Bleaching_Before_sim.plus = NA, Mean_Bleaching_Before_sim2 = NA, SD_Bleaching_Before_sim2 = NA, Mean_Bleaching_Before_sim.plus2 = NA, SD_Bleaching_Before_sim.plus2 = NA) # Initialize a dataframe with columns for simulated before data
ElNino_afteronly<-cbind(ElNino_afteronly,temp.df) # Bind together the temp.data.frame and the current afteronly dataframe
ElNino_afteronly$Mean_Bleaching_Before_sim[is.na(ElNino_afteronly$Mean_Bleaching_Before_sim)] <-"5" # Fill all of Mean_Bleaching_Before_sim with 5. This represents a 5% cover baseline for bleaching 
ElNino_afteronly$Mean_Bleaching_Before_sim <- as.numeric(ElNino_afteronly$Mean_Bleaching_Before_sim) # Make sure it's numeric
ElNino_afteronly$SD_Bleaching_Before_sim[is.na(ElNino_afteronly$SD_Bleaching_Before_sim)] <-"15" # Fill all of SD_Bleaching_Before_sim with 15. This represents a 15% SD from the mean (5%)
ElNino_afteronly$SD_Bleaching_Before_sim <- as.numeric(ElNino_afteronly$SD_Bleaching_Before_sim) # Make sure it's numeric
ElNino_afteronly$Mean_Bleaching_Before_sim.plus <- ElNino_afteronly$Mean_Bleaching_Before # Create column to input mean sim data into NAs
ElNino_afteronly$SD_Bleaching_Before_sim.plus <- ElNino_afteronly$SD_Bleaching_Before # Create column to input SD sim data into NAs
ElNino_afteronly$Mean_Bleaching_Before_sim.plus[is.na(ElNino_afteronly$Mean_Bleaching_Before_sim.plus)] <-"5" # Insert sim data to NAs
ElNino_afteronly$Mean_Bleaching_Before_sim.plus <- as.numeric(ElNino_afteronly$Mean_Bleaching_Before_sim.plus) # Make sure it's numeric
ElNino_afteronly$SD_Bleaching_Before_sim.plus[is.na(ElNino_afteronly$SD_Bleaching_Before_sim.plus)] <-"15" # Insert sim data to NAs
ElNino_afteronly$SD_Bleaching_Before_sim.plus <- as.numeric(ElNino_afteronly$SD_Bleaching_Before_sim.plus) # Make sure it's numeric

################################################

### Calculate the number of papers and the number of datapoints in meta-analysis ###

ElNino_Excel_forcalc <- read.csv("data/ElNino_Data_forcalc.csv", header=TRUE, sep=",") # Read in csv summary file

# Calculate and summarize total number of papers evaluated
NTotal_incl <- length(unique(ElNino_FULL$ACC)) # Total number of papers included

StudiesBleachWBefore <- length(unique(bleaching$ACC)) # Count how many bleaching studies there are with before El Niño/La Niña data
StudiesBleachAfterOnly <- length(unique(ElNino_afteronly$ACC)) # Count how many bleaching stuies there are with no before data

## Calculate and summarize number of bleaching/cover papers evaluated
bleaching_incl <- length(unique(bleaching$ACC)) # Create table to calc number of bleaching papers included
cover_incl <- length(unique(cover$ACC)) # Create table to calc number of coral cover papers included

NCover <- nrow(cover) # Count how many rows are remaining

NBlAfterOnly <- nrow(ElNino_afteronly) # Count the number of after-only coral bleaching data points, includes data points that were included in "bleaching" (i.e. that had before data points)

### Calculate number of data points extracted with GraphClick software
graphclickY1<-cover[which(cover$Graphclick=="Y"), ] # Calc # of data extracted with GraphClick
graphclickY2<-ElNino_afteronly[which(ElNino_afteronly$Graphclick=="Y"), ] # Calc # of data extracted with GraphClick
graphclickY <- nrow(graphclickY1) + nrow(graphclickY2)
graphclickN1<-cover[which(cover$Graphclick=="N"), ] # Calc # of data extracted from text/tables
graphclickN2<-ElNino_afteronly[which(ElNino_afteronly$Graphclick=="N"), ] # Calc # of data extracted from text/tables
graphclickN <- nrow(graphclickN1) + nrow(graphclickN2)

graphclickY<-nrow(ElNino_FULL[which(ElNino_FULL$Graphclick=="Y"), ]) # Calc # of data extracted with GraphClick
graphclickN<-nrow(ElNino_FULL[which(ElNino_FULL$Graphclick=="N"), ]) # Calc # of data extracted from text/tables

################################################

### Calculate number of papers excluded in each step for flowchart
ElNino_Excel_NoDup <- ElNino_Excel_forcalc[which(ElNino_Excel_forcalc$Duplicate.Of=="N"), ] # Subset papers, excluding duplicates
NTotal_NoDup <- nrow(ElNino_Excel_NoDup) # Number of papers

ElNino_Excel_CanAccess <- ElNino_Excel_NoDup[which(ElNino_Excel_NoDup$Cant.Access=="N"), ] # Subset papers, excluding papers which are inaccessible
NTotal_CanAccess <- nrow(ElNino_Excel_CanAccess) # Number of papers

ElNino_Excel_NoRevMA <- ElNino_Excel_CanAccess[which(ElNino_Excel_CanAccess$Review.MetaAnalysis=="N"), ] #Subset papers, excluding reviews and meta-analyses
NTotal_NoRevMA <- nrow(ElNino_Excel_NoRevMA) # Number of papers

ElNino_Excel_PrLit <- ElNino_Excel_NoRevMA[which(ElNino_Excel_NoRevMA$NotPrimaryLit=="N"), ] #Subset papers, excluding papers that are not primary lit or ICRS Proceedings
NTotal_PrLit <- nrow(ElNino_Excel_PrLit) # Number of papers

ElNino_Excel_NoPseudo <- ElNino_Excel_PrLit[which(ElNino_Excel_PrLit$Data.elsewhere=="N"), ] # Subset papers, excluding papers that have the same data as other papers
NTotal_NoPseudo <- nrow(ElNino_Excel_NoPseudo) # Number of papers

ElNino_Excel_PossRel <- ElNino_Excel_NoPseudo[which(ElNino_Excel_NoPseudo$Possibly.Relevant=="Y"), ] # Subset papers, excluding irrelevant papers; i.e. not about El NiÃ±o/La NiÃ±a and coral cover/bleaching responses
NTotal_PossRel <- nrow(ElNino_Excel_PossRel) # Number of papers

ElNino_Excel_ElNino <- ElNino_Excel_PossRel[which(ElNino_Excel_PossRel$Not.El.Nino!="Y"), ] # Subset papers, excluding those that are not during an El NiÃ±o event
NTotal_ElNino <- nrow(ElNino_Excel_ElNino) # Number of papers

ElNino_Excel_Cover <- ElNino_Excel_ElNino[which(ElNino_Excel_ElNino$Cover=="Y"), ] # Subset papers, cover
NTotal_Cover <- nrow(ElNino_Excel_Cover) # Number of papers

ElNino_Excel_Bleaching <- ElNino_Excel_ElNino[which(ElNino_Excel_ElNino$Bleaching=="Y"), ] # Subset papers, bleaching
NTotal_Bleaching <- nrow(ElNino_Excel_Bleaching) # Number of papers

ElNino_Excel_Cover_quant <- ElNino_Excel_Cover[which(ElNino_Excel_Cover$Qualitative!="Y"), ] # Subset papers, remove qualitative studies, keep only quantitative studies
NTotal_Cover_quant <- nrow(ElNino_Excel_Cover_quant) # Number of papers

ElNino_Excel_Bleaching_quant <- ElNino_Excel_Bleaching[which(ElNino_Excel_Bleaching$Qualitative!="Y"), ] # Subset papers, remove qualitative studies, keep only quantitative studies
NTotal_Bleaching_quant <- nrow(ElNino_Excel_Bleaching_quant) # Number of papers

ElNino_Excel_CovWBef <- ElNino_Excel_Cover_quant[which(ElNino_Excel_Cover_quant$Cover.no.before=="N"), ] # Subset papers, remove coral cover studies that do not provide before-impact data
NTotal_CovWBef <- nrow(ElNino_Excel_CovWBef)  # Number of papers

ElNino_Excel_HasNandSD_cov <- ElNino_Excel_CovWBef[which(ElNino_Excel_CovWBef$NoNorSD=="N"), ] # Subset papers, excluding papers that do not have sample size (N) and/or error (SD or SE)
NTotal_HasNandSD_cov <- nrow(ElNino_Excel_HasNandSD_cov) # Number of papers
MissingNorSD.cov <- ElNino_Excel_CovWBef[which(ElNino_Excel_CovWBef$NoNorSD=="Y"), ]
NTotal_MissingNorSD.cov <- nrow(MissingNorSD.cov) # Number of papers

ElNino_Excel_BlWBef <- ElNino_Excel_Bleaching_quant[which(ElNino_Excel_Bleaching_quant$Bleaching.Wbefore=="Y"), ] # Subset papers, studies that have before-El NiÃ±o/La NiÃ±a quantification of bleaching prevalence
NTotal_BlWBef <- nrow(ElNino_Excel_BlWBef)  # Number of papers

ElNino_Excel_HasNandSD_BlWBef <- ElNino_Excel_BlWBef[which(ElNino_Excel_BlWBef$NoNorSD=="N"), ] # Subset papers, excluding papers that do not have sample size (N) and/or error (SD or SE)
NTotal_HasNandSD_BlWBef <- nrow(ElNino_Excel_HasNandSD_BlWBef) # Number of papers
MissingNorSD.bl.wb <- ElNino_Excel_BlWBef[which(ElNino_Excel_BlWBef$NoNorSD=="Y"), ]
NTotal_MissingNorSD.bl.wb <- nrow(MissingNorSD.bl.wb)

ElNino_Excel_BlAfterOnly <- ElNino_Excel_Bleaching_quant[which(ElNino_Excel_Bleaching_quant$Bleaching.Wbefore=="N"), ] # Subset papers, remove qualitative studies, keep only quantitative studies
NTotal_BlAfterOnly <- nrow(ElNino_Excel_BlAfterOnly)  # Number of papers

bl.BeforeandAfter <- NTotal_BlWBef+NTotal_BlAfterOnly

ElNino_Excel_HasNandSD_BlAfterOnly <- ElNino_Excel_BlAfterOnly[which(ElNino_Excel_BlAfterOnly$NoNorSD=="N"), ] # Subset papers, excluding papers that do not have sample size (N) and/or error (SD or SE)
NTotal_HasNandSD_BlAfterOnly <- nrow(ElNino_Excel_HasNandSD_BlAfterOnly) # Number of papers
MissingNorSD.bl.ao <- ElNino_Excel_BlAfterOnly[which(ElNino_Excel_BlAfterOnly$NoNorSD=="Y"), ]
NTotal_MissingNorSD.bl.ao <- nrow(MissingNorSD.bl.ao)

ElNino_Excel_cov_pcovonly <- ElNino_Excel_HasNandSD_cov[which(ElNino_Excel_HasNandSD_cov$Units.are.Density!="Y"), ] # Subset papers, excluding studies which measure effects in density (rather than percent cover or similar metric)
NTotal_ElNino_Excel_cov_pcovonly <- nrow(ElNino_Excel_cov_pcovonly)

ElNino_Excel_cov_LATE <- ElNino_Excel_cov_pcovonly[which(ElNino_Excel_cov_pcovonly$Recovery..2..years.!="Y"), ] # Subset papers, excluding papers quantifying changes >2 years after the El NiÃ±o/La NiÃ±a event.
NTotal_ElNino_Excel_cov_LATE <- nrow(ElNino_Excel_cov_LATE)

ElNino_Excel_BlWBef_pblonly <- ElNino_Excel_HasNandSD_BlWBef[which(ElNino_Excel_HasNandSD_BlWBef$Units.are.Density!="Y"), ] # Subset papers, excluding studies which measure effects in density (rather than percent cover or similar metric)
NTotal_ElNino_Excel_BlWBef_pblonly <- nrow(ElNino_Excel_BlWBef_pblonly)

ElNino_Excel_BlWBef_LATE <- ElNino_Excel_BlWBef_pblonly[which(ElNino_Excel_BlWBef_pblonly$Recovery..2..years.!="Y"), ] # Subset papers, excluding papers quantifying changes >2 years after the El NiÃ±o/La NiÃ±a event.
NTotal_ElNino_Excel_BlWBef_LATE <- nrow(ElNino_Excel_BlWBef_LATE)

ElNino_Excel_BlAfterOnly_pblonly <- ElNino_Excel_HasNandSD_BlAfterOnly[which(ElNino_Excel_HasNandSD_BlAfterOnly$Units.are.Density!="Y"), ] # Subset papers, excluding studies which measure effects in density (rather than percent cover or similar metric)
NTotal_ElNino_Excel_BlAfterOnly_pblonly <- nrow(ElNino_Excel_BlAfterOnly_pblonly)

ElNino_Excel_BlAfterOnly_LATE <- ElNino_Excel_BlAfterOnly_pblonly[which(ElNino_Excel_BlAfterOnly_pblonly$Recovery..2..years.!="Y"), ] # Subset papers, excluding papers quantifying changes >2 years after the El NiÃ±o/La NiÃ±a event.
NTotal_ElNino_Excel_BlAfterOnly_LATE <- nrow(ElNino_Excel_BlAfterOnly_LATE)

################################################

## Calculate the effect size (Hedges D) for each parameter
#Cover
cover.hedges <- escalc(measure="SMD",m2i=Mean_Cover_Before, m1i=Mean_Cover_Resp, sd2i=SD_Cover_Before, sd1i=SD_Cover_Resp, n2i=N_Before, n1i=N_During.After, data=cover, append=TRUE, vtype="UB")
cover.hedges <- cover.hedges[!is.na(cover.hedges$yi), ] # Remove rows with nas in yi because they will be removed by the model anyway
cover.hedges <- cover.hedges[!is.na(cover.hedges$MaxDHW), ] # Remove rows with nas in MaxDHW because cannot be used in the full model
cover.hedges <- cover.hedges[!is.na(cover.hedges$TimeLag), ] # Remove rows with nas in TimeLag because cannot be used in the full model

# Subset datapoints collected <1 year after maximum heat stress - for sensitivity testing
cover.hedges.1yr <- cover.hedges[which(cover.hedges$TimeLag_beforeafter<=365), ]

#Bleaching
#Bleaching After-only, before is all simulated data
bleaching.ao.sim.hedges <- escalc(measure="SMD",m2i=Mean_Bleaching_Before_sim, m1i=Mean_Bleaching_Resp, sd2i=SD_Bleaching_Before_sim, sd1i=SD_Bleaching_Resp, n2i=N_During.After, n1i=N_During.After, data=ElNino_afteronly, append=TRUE, vtype="UB")
bleaching.ao.sim.hedges <-bleaching.ao.sim.hedges[!is.na(bleaching.ao.sim.hedges$yi), ] # Remove rows with nas in yi because they will be removed by the model anyway

#Bleaching After-only, before includes any true "before" data, with any NAs filled with simulated data
bleaching.ao.sim.plus.hedges <- escalc(measure="SMD",m2i=Mean_Bleaching_Before_sim.plus, m1i=Mean_Bleaching_Resp, sd2i=SD_Bleaching_Before_sim.plus, sd1i=SD_Bleaching_Resp, n2i=N_During.After, n1i=N_During.After, data=ElNino_afteronly, append=TRUE, vtype="UB")
bleaching.ao.sim.plus.hedges <-bleaching.ao.sim.plus.hedges[!is.na(bleaching.ao.sim.plus.hedges$yi), ] # Remove rows with nas in yi because they will be removed by the model anyway

################################################

## Calculate Fail-Safe N (FSN) for each parameter to test for publication bias
#Cover
fsn.cover <- fsn(yi, vi, data=cover.hedges, type="Rosenberg") # Calc FSN for overall coral cover
fsn.cover.N <- as.numeric(fsn.cover[2]) # Extract fail-safe number for text

#Bleaching
fsn.bleaching.sim <- fsn(yi, vi, data=bleaching.ao.sim.hedges, type="Rosenberg") # Calc FSN for all sim bleaching data
fsn.bleaching.sim.N <- as.numeric(fsn.bleaching.sim[2]) # Extract fail-safe number for text

fsn.bleaching.sim.plus <- fsn(yi, vi, data=bleaching.ao.sim.plus.hedges, type="Rosenberg") # Calc FSN for partially simulated bleaching data
fsn.bleaching.sim.plus.N <- as.numeric(fsn.bleaching.sim.plus[2]) # Extract fail-safe number for text

################################################

# %-Cover
# Create column called id to include in the random factor (creates unique id)
cover.hedges$id <- 1:nrow(cover.hedges)
cover.hedges.1yr$id <- 1:nrow(cover.hedges.1yr)
# BACKWARDS STEPWISE ANOVAS: Start with the full model, including MaxDHW (maximum degree heat week), TimeLag (time lag between maximum DHW and sampling timepoint), SSTmean (long-term temperature mean) and SSTvar (long-term temperature variance), and all interactions between these moderators. Note that no categorical moderators (e.g. El.Nino.simple, Ocean.Region.simple, etc) are included in the full model because we do not have enough data in order to properly include these moderators.
cover.rma001 <- rma.mv(yi~MaxDHW*TimeLag+MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
cover.rma002 <- rma.mv(yi~MaxDHW*TimeLag+MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean+TimeLag*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
anova(cover.rma001,cover.rma002) # Compare these two models. Better
cover.rma003 <- rma.mv(yi~MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean+TimeLag*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
anova(cover.rma002,cover.rma003) # Compare these two models. Better
cover.rma004 <- rma.mv(yi~MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
anova(cover.rma003,cover.rma004) # Compare these two models. Better
cover.rma005 <- rma.mv(yi~MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
anova(cover.rma004,cover.rma005) # Compare these two models. Better
cover.rma006 <- rma.mv(yi~MaxDHW*SSTmean+MaxDHW*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
anova(cover.rma005,cover.rma006) # Compare these two models. Better
cover.rma007 <- rma.mv(yi~MaxDHW+SSTmean+MaxDHW*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
anova(cover.rma006,cover.rma007) # Compare these two models. WORSE
cover.rma008 <- rma.mv(yi~MaxDHW*SSTmean+SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
anova(cover.rma006,cover.rma008) # Compare these two models. Better 
cover.rma009 <- rma.mv(yi~MaxDHW*SSTmean, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
anova(cover.rma008,cover.rma009) # Compare these two models. Better - BEST MODEL
cover.rma010 <- rma.mv(yi~MaxDHW+SSTmean, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
anova(cover.rma009,cover.rma010) # Compare these two models. Worse 
cover.rma011 <- rma.mv(yi~SSTmean, vi, random=~1|ACC/id, method="ML", data=cover.hedges)
anova(cover.rma010,cover.rma011) # Compare these two models. Worse 

# First year of coral cover data only
cover.1yr.rma001 <- rma.mv(yi~MaxDHW*TimeLag+MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
cover.1yr.rma002 <- rma.mv(yi~MaxDHW*TimeLag+MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma001,cover.1yr.rma002) # Compare these two models. Better
cover.1yr.rma003 <- rma.mv(yi~MaxDHW*TimeLag+MaxDHW*SSTmean+MaxDHW*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma002,cover.1yr.rma003) # Compare these two models. Worse
cover.1yr.rma004 <- rma.mv(yi~MaxDHW*TimeLag+MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTvar+SSTmean, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma002,cover.1yr.rma004) # Compare these two models. Worse
cover.1yr.rma005 <- rma.mv(yi~TimeLag+MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma002,cover.1yr.rma005) # Compare these two models. Better
cover.1yr.rma006 <- rma.mv(yi~TimeLag+SSTmean+MaxDHW*SSTvar+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma005,cover.1yr.rma006) # Compare these two models. Worse
cover.1yr.rma007 <- rma.mv(yi~TimeLag+MaxDHW*SSTmean+SSTvar+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma005,cover.1yr.rma007) # Compare these two models. Worse
cover.1yr.rma008 <- rma.mv(yi~TimeLag+MaxDHW*SSTmean+MaxDHW*SSTvar+SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma005,cover.1yr.rma008) # Compare these two models. Better
cover.1yr.rma009 <- rma.mv(yi~TimeLag+MaxDHW*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma008,cover.1yr.rma009) # Compare these two models. Better
cover.1yr.rma010 <- rma.mv(yi~MaxDHW*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma009,cover.1yr.rma010) # Compare these two models. Better
cover.1yr.rma011 <- rma.mv(yi~MaxDHW+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma010,cover.1yr.rma011) # Compare these two models. Better
cover.1yr.rma012 <- rma.mv(yi~MaxDHW+SSTmean+SSTvar, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma011,cover.1yr.rma012) # Compare these two models. Better
cover.1yr.rma013 <- rma.mv(yi~MaxDHW+SSTmean, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma012,cover.1yr.rma013) # Compare these two models. Better. BEST MODEL
cover.1yr.rma014 <- rma.mv(yi~MaxDHW, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma013,cover.1yr.rma014) # Compare these two models. Worse
cover.1yr.rma015 <- rma.mv(yi~SSTmean, vi, random=~1|ACC/id, method="ML", data=cover.hedges.1yr)
anova(cover.1yr.rma013,cover.1yr.rma015) # Compare these two models. Worse


# Calculate overall effect size using random effects model
cover.rma <- rma.mv(yi, vi, random=~1|ACC/id, method="REML", data=cover.hedges) # All cover data (up to 2 years after heat stress)
cover.rma.1yr <- rma.mv(yi, vi, random=~1|ACC/id, method="REML", data=cover.hedges.1yr) # Up to 1 year after heat stress

################################################

# %-Bleaching
# Create column called id to include in the random factor (creates unique id)
bleaching.ao.sim.hedges$id <- 1:nrow(bleaching.ao.sim.hedges)
bleaching.ao.sim.plus.hedges$id <- 1:nrow(bleaching.ao.sim.plus.hedges)
# BACKWARDS STEPWISE ANOVAS: Start with the full model, including MaxDHW (maximum degree heat week), TimeLag (time lag between maximum DHW and sampling timepoint), fam.comm (community data vs. taxa-specific data) and Lat.abs (absolute latitude), and all interactions between these moderators. Note that no categorical moderators (e.g. El.Nino.simple, Ocean.Region.simple, etc) are included in the full model because we do not have enough data in order to properly include these moderators.
bleaching.rma001 <- rma.mv(yi~MaxDHW*TimeLag+MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
bleaching.rma002 <- rma.mv(yi~MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma001,bleaching.rma002) # Compare these two models. better
bleaching.rma003 <- rma.mv(yi~MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma002,bleaching.rma003) # Compare these two models. better
bleaching.rma004 <- rma.mv(yi~MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma003,bleaching.rma004) # Compare these two models. worse
bleaching.rma005 <- rma.mv(yi~MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma003,bleaching.rma005) # Compare these two models. better
bleaching.rma006 <- rma.mv(yi~SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma005,bleaching.rma006) # Compare these two models. better
bleaching.rma007 <- rma.mv(yi~MaxDHW+SSTvar+TimeLag*SSTmean, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma006,bleaching.rma007) # Compare these two models. better
bleaching.rma008 <- rma.mv(yi~SSTvar+TimeLag*SSTmean, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma007,bleaching.rma008) # Compare these two models. better 
bleaching.rma009 <- rma.mv(yi~SSTvar+TimeLag+SSTmean, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma008,bleaching.rma009) # Compare these two models. worse
bleaching.rma010 <- rma.mv(yi~TimeLag*SSTmean, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma008,bleaching.rma010) # Compare these two models. worse
bleaching.rma011 <- rma.mv(yi~TimeLag+SSTmean, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma009,bleaching.rma011) # Compare these two models. Better Timelag not sig.
bleaching.rma012 <- rma.mv(yi~SSTmean, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.hedges)
anova(bleaching.rma011,bleaching.rma012) # Compare these two models. worse

# Bleaching plus
bleachingplus.rma001 <- rma.mv(yi~MaxDHW*TimeLag+MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges) # Full model
bleachingplus.rma002 <- rma.mv(yi~MaxDHW*SSTmean+MaxDHW*SSTvar+TimeLag*SSTmean+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges)
anova(bleachingplus.rma001,bleachingplus.rma002) # Compare these two models. better
bleachingplus.rma003 <- rma.mv(yi~MaxDHW*SSTvar+TimeLag*SSTmean+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges)
anova(bleachingplus.rma002,bleachingplus.rma003) # Compare these two models. better
bleachingplus.rma004 <- rma.mv(yi~MaxDHW*SSTvar+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges)
anova(bleachingplus.rma003,bleachingplus.rma004) # Compare these two models. Better
bleachingplus.rma005 <- rma.mv(yi~TimeLag*SSTmean+TimeLag*SSTvar+SSTmean*SSTvar, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges)
anova(bleachingplus.rma004,bleachingplus.rma005) # Compare these two models. better
bleachingplus.rma006 <- rma.mv(yi~TimeLag*SSTmean+TimeLag*SSTvar, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges)
anova(bleachingplus.rma005,bleachingplus.rma006) # Compare these two models. better
bleachingplus.rma007 <- rma.mv(yi~SSTvar+SSTmean*TimeLag, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges)
anova(bleachingplus.rma006,bleachingplus.rma007) # Compare these two models. better
bleachingplus.rma008 <- rma.mv(yi~SSTvar+SSTmean+TimeLag, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges)
anova(bleachingplus.rma007,bleachingplus.rma008) # Compare these two models. better
bleachingplus.rma009 <- rma.mv(yi~SSTmean+TimeLag, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges)
anova(bleachingplus.rma008,bleachingplus.rma009) # Compare these two models. better 
bleachingplus.rma010 <- rma.mv(yi~TimeLag, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges)
anova(bleachingplus.rma009,bleachingplus.rma010) # Compare these two models. better BEST MODEL
bleachingplus.rma011 <- rma.mv(yi~SSTmean, vi, random=~1|ACC/id, method="ML", data=bleaching.ao.sim.plus.hedges)
anova(bleachingplus.rma009,bleachingplus.rma011) # Compare these two models. Worse

# Calculate overall effect size using random effects model
bleaching.rma <- rma.mv(yi, vi, random=~1|ACC/id, method="REML", data=bleaching.ao.sim.hedges)
bleaching.plus.rma <- rma.mv(yi, vi, random=~1|ACC/id, method="REML", data=bleaching.ao.sim.plus.hedges)

################################################

# glmulti Model Selection
#Define rma.glmulti function (from Metafor website, http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti)
rma.glmulti <- function(formula, data,exclude=c(),random=list(""),...) {
rma.mv(as.formula(paste(deparse(formula))),random=random, vi,
data=data, method="ML", ...) # Need ML to be able to compare among models (as opposed to REML)
}

cov <- glmulti(yi~MaxDHW+TimeLag+SSTmean+SSTvar, random=list("~1|ACC/id"), # First try the full model with all 4 potential moderators
data=cover.hedges, fitfunction = rma.glmulti,
level=2, crit="aicc") # Level 2 looks at the interactions between all moderators. AICc allows for correction of a small N
print(cov) # View summary of cov model
tmp <- weightable(cov) # Save tmp weightable for cov
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 10000,]
tmp # print AICc and weights for all models
cov.summ <- summary(cov@objects[[1]]) # create summary object for data extraction
cov.lb <-cov.summ$ci.lb # Extract confidence interval - lower bound
cov.ub <-cov.summ$ci.ub # Extract confidence interval - upper bound
cov.par <- row.names(cov.summ$b)
cov.lb <- rbind(cov.par, cov.lb)
cov.ub <- rbind(cov.par, cov.ub)
cov.est <-cov.summ$b # Extract estimate

# Test for subset of studies only conducted within one year of maximum heat stress event
cov.1yr <- glmulti(yi~MaxDHW+TimeLag+SSTmean+SSTvar, random=list("~1|ACC/id"), # First try the full model with all 4 potential moderators
data=cover.hedges.1yr, fitfunction = rma.glmulti,
level=2, crit="aicc") # Level 2 looks at the interactions between all moderators. AICc allows for correction of a small N
print(cov.1yr) # View summary of cov.1yr model
tmp <- weightable(cov.1yr) # Save tmp weightable
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 10000,]
tmp # print AICc and weights for all models
cov.1yr.summ <- summary(cov.1yr@objects[[1]]) # create summary object for data extraction
cov.1yr.lb <-cov.1yr.summ$ci.lb # Extract confidence interval - lower bound
cov.1yr.ub <-cov.1yr.summ$ci.ub # Extract confidence interval - upper bound
cov.1yr.par <- row.names(cov.1yr.summ$b)
cov.1yr.lb <- rbind(cov.1yr.par, cov.1yr.lb)
cov.1yr.ub <- rbind(cov.1yr.par, cov.1yr.ub)
cov.1yr.est <-cov.1yr.summ$b # Extract estimate

bl <- glmulti(yi~MaxDHW+TimeLag+SSTmean+SSTvar, random=list("~1|ACC/id"), # First try the full model with all 3 potential moderators
data=bleaching.ao.sim.hedges, fitfunction = rma.glmulti,
level=2, crit="aicc") # Level 2 looks at the interactions between all moderators. AICc allows for correction of a small N
print(bl) # View summary of bl model
tmp <- weightable(bl) # Save tmp weightable for bl
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 10000,]
tmp # print AICc and weights for all models
bl.summ <- summary(bl@objects[[1]]) # create summary object for data extraction
bl.lb <- bl.summ$ci.lb # Extract confidence interval - lower bound
bl.ub <- bl.summ$ci.ub # Extract confidence interval - upper bound
bl.par <- row.names(bl.summ$b)
bl.lb <- rbind(bl.par, bl.lb)
bl.ub <- rbind(bl.par, bl.ub)
bl.est <- bl.summ$b # Extract estimate

blplus <- glmulti(yi~MaxDHW+TimeLag+SSTmean+SSTvar, random=list("~1|ACC/id"), # First try the full model with all 3 potential moderators
data=bleaching.ao.sim.plus.hedges, fitfunction = rma.glmulti,
level=2, crit="aicc") # Level 2 looks at the interactions between all moderators. AICc allows for correction of a small N
print(blplus) # View summary of blplus model
tmp <- weightable(blplus) # Save tmp weightable for cov2
tmp <- tmp[tmp$aicc <= min(tmp$aicc) + 10000,]
tmp # print AICc and weights for all models
bl.plus.summ <- summary(blplus@objects[[1]]) # create summary object for data extraction
bl.plus.lb <- bl.plus.summ$ci.lb # Extract confidence interval - lower bound
bl.plus.ub <- bl.plus.summ$ci.ub # Extract confidence interval - upper bound
bl.plus.par <- row.names(bl.plus.summ$b)
bl.plus.lb <- rbind(bl.plus.par, bl.plus.lb)
bl.plus.ub <- rbind(bl.plus.par, bl.plus.ub)
bl.plus.est <- bl.plus.summ$b # Extract estimate

################################################

# Extract glmulti model output for table
cov.1yr.QM <- round(as.numeric(cov.1yr.summ$QM),2) #Extract QM and round to 2 sig figs
cov.1yr.QMp <- round(as.numeric(cov.1yr.summ$QMp),3) #Extract QM p-value and round to 3 sig figs
if (cov.1yr.QMp == "0"){ # If the rounded value of QMp is "0"
cov.1yr.QMp <- "***" # Report it as ***
} else if (cov.1yr.QMp < "0.01"){ # Otherwise, if 0 < cov.1yr.QMp < 0.01
cov.1yr.QMp <- "**" # Report it as **
} else if (cov.1yr.QMp < "0.05"){ # Otherwise, if 0.01 < cov.1yr.QMp < 0.05
cov.1yr.QMp <- "*" # Report it as *
} else { # Otherwise
cov.1yr.QMp <- "ERROR" # Return "ERROR"
}
cov.1yr.QM <- paste(cov.1yr.QM,cov.1yr.QMp,sep="")
cov.1yr.QMdf=data.frame(2) # create a dataframe with df 
colnames(cov.1yr.QMdf) = "df" # Name column in dataframe
cov.1yr.QE <- round(as.numeric(cov.1yr.summ$QE),0) # Extract QE and round to whole number
cov.1yr.QEp <- round(as.numeric(cov.1yr.summ$QEp),3) #Extract QE p-value and round to 3 sig figs
if (cov.1yr.QEp == "0"){ # If the rounded value of QMp is "0"
cov.1yr.QEp <- "***" # Report it as ***
} else if (cov.1yr.QEp < "0.01"){ # Otherwise, if 0 < cov.1yr.QEp < 0.01
cov.1yr.QEp <- "**" # Report it as **
} else if (cov.1yr.QEp < "0.05"){ # Otherwise, if 0.01 < cov.1yr.QEp < 0.05
cov.1yr.QEp <- "*" # Report it as *
} else { # Otherwise
cov.1yr.QEp <- "ERROR" # Return "ERROR"
}
cov.1yr.QE <- paste(cov.1yr.QE,cov.1yr.QEp,sep="")
cov.1yr.QEdf=data.frame(153) # create a dataframe with df
colnames(cov.1yr.QEdf) = "df" # Name column in dataframe
cov.1yr.table <- data.frame("Cover","MaxDHW, Mean Temp",round(cov.1yr.summ$fit.stats$ML[5],2), cov.1yr.QM,cov.1yr.QMdf,cov.1yr.QE,cov.1yr.QEdf) # create a dataframe 
colnames(cov.1yr.table) <- c("Model","Moderators (Top Model)","AICc","QM","df","QE","df") # Name columns in dataframe

bl.QM <- round(as.numeric(bl.summ$QM),0) # Extract QM and round to 2 sig figs
bl.QMp <- round(as.numeric(bl.summ$QMp),3) #Extract QM p-value and round to 3 sig figs
if (bl.QMp == "0"){ # If the rounded value of QMp is "0"
bl.QMp <- "***" # Report it as ***
} else if (bl.QMp < "0.01"){ # Otherwise, if 0 < bl.QMp < 0.01
bl.QMp <- "**" # Report it as **
} else if (bl.QMp < "0.05"){ # Otherwise, if 0.01 < bl.QMp < 0.05
bl.QMp <- "*" # Report it as *
} else { # Otherwise
bl.QMp <- "ERROR" # Return "ERROR"
}
bl.QM <- paste(bl.QM,bl.QMp,sep="")
QMdf=data.frame(2) # create a dataframe with df
colnames(QMdf) = "df" # Name column in dataframe
bl.QE <- round(as.numeric(bl.summ$QE),0) # Extract QE and round to whole number
bl.QEp <- round(as.numeric(bl.summ$QEp),3) #Extract QE p-value and round to 3 sig figs
if (bl.QEp == "0"){ # If the rounded value of QMp is "0"
bl.QEp <- "***" # Report it as ***
} else if (bl.QEp < "0.01"){ # Otherwise, if 0 < bl.QEp < 0.01
bl.QEp <- "**" # Report it as **
} else if (bl.QEp < "0.05"){ # Otherwise, if 0.01 < bl.QEp < 0.05
bl.QEp <- "*" # Report it as *
} else { # Otherwise
bl.QEp <- "ERROR" # Return "ERROR"
}
bl.QE <- paste(bl.QE,bl.QEp,sep="")
QEdf=data.frame(140) # create a dataframe with df
colnames(QEdf) = "df" # Name column in dataframe
bl.table <- data.frame("Bleaching","Time Lag, Mean Temp",round(bl.summ$fit.stats$ML[5],2), bl.QM,QMdf,bl.QE,QEdf) # create dataframe
colnames(bl.table) <- c("Model","Moderators (Top Model)","AICc","QM","df","QE","df") # name columns in dataframe

full.table <- rbind(bl.table,cov.1yr.table) # create a full table of all results
full.table # print full table

################################################

# Proportional reduction in each variance component, i.e. R^2 for the study level and the estimate level. 
# "Total" r^2
cover.r2 <- (sum(cover.rma.1yr$sigma2) - sum(cov.1yr.summ$sigma2)) / sum(cover.rma.1yr$sigma2)*100

bleaching.r2 <- (sum(bleaching.rma$sigma2) - sum(bl.summ$sigma2)) / sum(bleaching.rma$sigma2)*100

################################################

### Make profile plots

tiff(filename = "ms/PLoS/cover_profileplots.tiff", # Open tiff file
     width = 6, height= 6, units = "in", # Set tiff size
     compression = "lzw", res = 300) # Set compression and resolution
par(mfrow=c(2,1), mar=c(4,4.5,3,1)) # Set graphical parameters
cover.sigma <- profile(cov.1yr.summ) # Create profile plots
dev.off() # Close tiff file to complete

tiff(filename = "ms/PLoS/bleaching_profileplots.tiff", # Open tiff file
     width = 6, height= 6, units = "in", # Set tiff size
     compression = "lzw", res = 300) # Set compression and resolution
par(mfrow=c(2,1), mar=c(4,4.5,3,1))# Set graphical paramters
bleaching.sigma <- profile(bl.summ) # Create profile plots
dev.off() # Close tiff file to complete

################################################

# Save current workspace
save.image(file="data/ElNino.RData")
