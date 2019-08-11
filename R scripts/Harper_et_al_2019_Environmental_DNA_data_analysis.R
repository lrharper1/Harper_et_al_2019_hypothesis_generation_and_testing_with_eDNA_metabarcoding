#' ---
#' Title: "Generating and testing ecological hypotheses at the pondscape with environmental DNA metabarcoding: a case study on a threatened amphibian"
#' Author: "Lynsey Rebecca Harper"
#' Date: "1st August 2019"
#' ---
#' 
#' 
#' A total of 532 samples previously analysed by qPCR for great crested newt 
#' (GCN) were re-analysed by eDNA metabarcoding to evaluate this method for 
#' detection of vertebrate communities in freshwater ponds.
#' 
#' Environmental metadata was also obtained for these samples from 
#' Natural England and applied to the species inventories obtained by
#' eDNA metabarcoding.
#' 
#' This analysis aims to generate and test ecological hypotheses using eDNA 
#' metabarcoding, with the GCN as a focal species due to extensive literature 
#' on its ecology. We aim to test biotic and abiotic factors associated with 
#' GCN detection, and test the umbrella species status of GCN.
#' 
#' 
#' ## Prepare working environment
#' 
#' Clear R memory, set working directory and load required packages. Then,
#' load functions for calculating kappa coefficient, plotting model 
#' residuals and testing model fit.
#' 

## Clear memory
rm(list=ls())

## set working directory to the location of the script
#install.packages("rstudioapi") # first time for each computer
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Check working directory
getwd()

## Load required packages
p <- c("ggplot2","grid","gridExtra","forcats","lme4","glmmADMB","MASS","car",
       "scales","AICcmodavg","gtools","reshape2","plyr","dplyr","arm",
       "RVAideMemoire","ResourceSelection","bbmle","RColorBrewer", 
       "MuMIn","rpart","mgcv","ncf","glmmML","digest",
       "LMERConvenienceFunctions")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, 
                                          repos="https://cran.cnr.berkeley.edu/",
                                          dependencies=TRUE)
lapply(p, require, character.only = TRUE)

## Load custom functions
f <- c("CheckResidsFunction.R","OverdispersalFunction.R",
       "CheckConvergenceFunction.R", "HighstatLibV6.R",
       "MyLibrary.R", "glmmML_pres_function.R",
       "ggplotLegendFunction.R", "grid_arrange_shared_legend_Function.R")
lapply(f, source)

#'
#' To ensure reproducibility, print details about the version of R being
#' used for analysis.
#' 

sessionInfo()


#' ---
#' 
#' ## 1) Raw data processing
#' 
#' Examine the raw data output from metaBEAT.
#' 

## Original BLAST 
ass.raw <- read.csv("../Data/GCN_metabarcoding_assigned_raw.csv", header=TRUE)
summary(ass.raw[,c(1:5,length(ass.raw))])
head(ass.raw)
names(ass.raw)
str(ass.raw)

## Unassigned BLAST
unass.raw <- read.csv("../Data/GCN_metabarcoding_unassigned_raw.csv", header=TRUE)
summary(unass.raw[,c(1:5,length(unass.raw))])
head(unass.raw)
names(unass.raw)
str(unass.raw)

## Make first row header names
names(ass.raw) <- lapply(ass.raw[1,], as.character)
names(unass.raw) <- lapply(unass.raw[1,], as.character)
ass.raw <- ass.raw[-1,]
unass.raw <- unass.raw[-1,]

## Reset row names
rownames(ass.raw) <- NULL
rownames(unass.raw) <- NULL

## Remove '-nc.blast'
colnames(ass.raw) <- gsub("_-nc.blast", "", colnames(ass.raw))
colnames(unass.raw) <- gsub("_-nc.blast.blast", "", colnames(unass.raw))

## Rename first column
colnames(ass.raw)[1] <- "Assignment"
colnames(unass.raw)[1] <- "Assignment"

## Correct samples known to have wrong tags placed on them for sequencing 
## causing environmental samples to be labelled controls and vice versa 
## according to sample sheet
## P4-3 = F82
## P4-4 = F84
## P4-5 = F86
## P4-6 = F88
## P6-3 = F132
## P6-4 = F133
## P6-5 = F135
## P6-6 = F137

colnames(ass.raw) <- as.character(colnames(ass.raw))
colnames(ass.raw)[495] <- "P4-3"
colnames(ass.raw)[497] <- "P4-4"
colnames(ass.raw)[499] <- "P4-5"
colnames(ass.raw)[501] <- "P4-6"
colnames(ass.raw)[736] <- "F82"
colnames(ass.raw)[737] <- "F84"
colnames(ass.raw)[738] <- "F86"
colnames(ass.raw)[739] <- "F88"
colnames(ass.raw)[42] <- "P6-3"
colnames(ass.raw)[43] <- "P6-4"
colnames(ass.raw)[45] <- "P6-5"
colnames(ass.raw)[47] <- "P6-6"
colnames(ass.raw)[748] <- "F132"
colnames(ass.raw)[749] <- "F133"
colnames(ass.raw)[750] <- "F135"
colnames(ass.raw)[751] <- "F137"

colnames(unass.raw) <- as.character(colnames(unass.raw))
colnames(unass.raw)[495] <- "P4-3"
colnames(unass.raw)[497] <- "P4-4"
colnames(unass.raw)[499] <- "P4-5"
colnames(unass.raw)[501] <- "P4-6"
colnames(unass.raw)[736] <- "F82"
colnames(unass.raw)[737] <- "F84"
colnames(unass.raw)[738] <- "F86"
colnames(unass.raw)[739] <- "F88"
colnames(unass.raw)[42] <- "P6-3"
colnames(unass.raw)[43] <- "P6-4"
colnames(unass.raw)[45] <- "P6-5"
colnames(unass.raw)[47] <- "P6-6"
colnames(unass.raw)[748] <- "F132"
colnames(unass.raw)[749] <- "F133"
colnames(unass.raw)[750] <- "F135"
colnames(unass.raw)[751] <- "F137"

## Reorder dataframes by column names
ass.raw <- ass.raw[,order(names(ass.raw))]
unass.raw <- unass.raw[,order(names(unass.raw))]

## Make columns factors again
colnames(ass.raw) <- as.factor(colnames(ass.raw))
colnames(unass.raw) <- as.factor(colnames(unass.raw))

## Remove any taxonomic assignments that aren't vertebrate from each dataframe
## using the taxonomy column created during processing with metaBEAT
ass <- ass.raw[which(grepl("Chordata|unassigned", ass.raw$taxomomy)),]
unass <- unass.raw[which(grepl("Chordata|unassigned", unass.raw$taxomomy)),]

## Remove last column containing taxonomy
ass <- ass[-770]
unass <- unass[-770]

## Bind data frames
merged.df <- rbind(ass, unass)

## Remove samples included from other projects
merged.df <- merged.df[,which(!grepl("B|JL", colnames(merged.df)))]

## Remove taxonomic assignments with no read counts across any samples
merged.df <- merged.df[!apply(merged.df == 0, 1, all),]

## Merge read counts by taxonomic assignment
## Anything with the same name will be merged and anything unique will
## be retained
merged.df[,2:761] <- lapply(merged.df[,2:761], function(x) as.numeric(as.character(x)))
merged.df <- ddply(merged.df, .(Assignment), numcolwise(sum))

## Export as .csv file
write.csv(merged.df, "../Data/GCN_metabarcoding_merged.csv", row.names=FALSE)


#' --- 
#' 
#' ## 2) Refine dataset
#' 
#' Now, we need to further refine the metabarcoding dataset.
#' 
#' 1. Any spurious species must be removed. Use NBN atlas to check species
#'    occurrence records and ensure they match with sampling locations. This
#'    is also a good source for checking current taxonomy.
#' 2. Any Genus or Family assignments containing only one species in 
#'    the UK must be changed to that species.
#' 3. Likewise, for any species assignment which is the only species in 
#'    the UK and also has genus/family reads, read counts from all
#'    assignments must be merged.
#'  
#' A record of these changes will be kept as these changes are made.
#'      

## Inspect new dataframe
summary(merged.df[,1:6])
names(merged.df)
str(merged.df)

#' First, remove spurious assignments as follows:
#' 
#' - Gadidae, row 98
#' - Gadiformes (marine cod), row 99
#' - Loricariidae (catfish), row 101
#' - *Ovis canadensis* (Bighorn sheep), row 105
#' - *Pan*, row 106
#' - *Saurida* (Lizardfishes), row 110
#' 

true.assign <- merged.df[-c(98:99,101,105:106,110),]

## Reset row names of data frame for further indexing
rownames(true.assign) <- NULL

## Remove underscore from species names
true.assign$Assignment <- gsub("_", " ", true.assign$Assignment)


#'
#' Now, correct species names:
#' 
#' - *Anas carolinensis* = *Anas*, row 3
#' - *Canis lupus* = *Canis lupus familiaris*, row 17
#' - *Felis silvestris* = *Felis catus*, row 32
#' - *Fulica* = *Fulica atra*, row 34
#' - *Sus scrofa* = *Sus scrofa domesticus*, row 84
#' - *Bos* = *Bos taurus*, row 91
#' - *Canis* = *Canis lupus familiaris*, row 93
#' - Cichlidae = *Rhamphochromis esox*, row 95
#' - *Pungitius* = *Pungitius pungitius*, row 104
#' - *Strix* = *Strix aluco*, row 105
#' 

true.assign$Assignment <- as.character(true.assign$Assignment)
true.assign[3, "Assignment"] <- "Anas"
true.assign[c(17,93), "Assignment"] <- "Canis lupus familiaris"
true.assign[32, "Assignment"] <- "Felis catus"
true.assign[34, "Assignment"] <- "Fulica atra"
true.assign[84, "Assignment"] <- "Sus scrofa domesticus"
true.assign[91, "Assignment"] <- "Bos taurus"
true.assign[95, "Assignment"] <- "Rhamphochromis esox"
true.assign[104, "Assignment"] <- "Pungitius pungitius"
true.assign[105, "Assignment"] <- "Strix aluco"

## Merge read counts again by taxonomic assignment
## Anything with the same name in the data frame will be merged
true.assign <- ddply(true.assign, .(Assignment), numcolwise(sum))


#' ---
#' 
#' ## 3) Clean up dataset
#' 
#' Now the data is in a form that can be manipulated easily, it must be 
#' filtered to remove potential contaminants and false positives. 
#' 
#' There are several ways of doing this:
#' 1. Identify highest level of cichlid DNA contamination across all eDNA
#'    samples.
#' 3. Identify the highest level of contamination in positive (non-cichlid
#'    DNA) controls.
#' 4. Identify taxon-specific thresholds using positive controls, i.e. 
#'    the frequency required to remove a given taxa from the positive 
#'    control as only cichlid DNA should be present.
#' 
#' Arguably, taxon-specific thresholds based on positive controls are 
#' more effective as negative controls have no template DNA for 
#' contaminant DNA to compete with for amplification, thus contaminant
#' DNA amplifies exponentially. However, only 16 positive controls were
#' included on this MiSeq run which may render this approach ineffective.
#'

####################################################
# OPTION 1: highest level of cichlid contamination #
####################################################

## Create copy of dataframe without positive or negative controls for 
## threshold determination
cichlid <- true.assign[,which(grepl("Assignment|F|NU", colnames(true.assign)))]

## Make Assignment column row names
rownames(cichlid) <- cichlid$Assignment
cichlid <- cichlid[,-1]

## Calculate total number of reads in samples
cichlid <- rbind(cichlid, colSums(cichlid))

## Create new dataframe containing the frequency of reads in each sample
cichlid.freq  <- cichlid/c(cichlid[99,])
cichlid.freq[is.na(cichlid.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency of
## DNA from each taxon across samples
cichlid.freq$Threshold <- apply(cichlid.freq, 1, max)

## Check Threshold column has been created properly
head(cichlid.freq[,530:533])

## Combine read frequencies with taxonomic assignment
species <- data.frame(true.assign$Assignment)
colnames(species) <- "Assignment"
total <- data.frame("Total")
colnames(total) <- "Assignment"
species <- rbind(species, total)
cichlid.freq <- cbind(species, cichlid.freq)
rownames(cichlid.freq) <- NULL

## Print contamination threshold based on max. level of cichlid
## DNA contamination
max(cichlid.freq[82,534])   # 1

## In this scenario, any assignments <100% total reads in biological 
## samples would be considered contamination. This threshold cannot be
## applied to the data as it is too stringent and would effectively remove
## all taxonomic assignments. However, we will plot for comparison with
## other thresholds.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
cichlid.test <- cichlid.freq
cichlid.test[cichlid.test <= 1] <- 0

## Now convert back into read counts.
## Remove last row containing frequencies and add the total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
cichlid.test <- cichlid.test[-99,-534]
rownames(cichlid.test) <- NULL

total.counts <- data.frame(colSums(cichlid[-99,]))
total.counts <- data.frame(t(total.counts))
colnames(total.counts) <- gsub('[.]', '-', colnames(total.counts))
total.counts <- cbind(total, total.counts)

## Now convert frequencies back to read counts
cichlid.conversion <- smartbind(cichlid.test, total.counts)
rownames(cichlid.conversion) <- cichlid.conversion$Assignment
cichlid.conversion <- cichlid.conversion[,-1]
cichlid.FP <- cichlid.conversion*c(cichlid.conversion[99,])

## Remove total row, reset row names, recreate Assignment column
cichlid.FP <- cichlid.FP[-99,]
cichlid.FP$Assignment <- rownames(cichlid.FP)
cichlid.FP <- cichlid.FP[,c(533,1:532)]
rownames(cichlid.FP) <- NULL

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- cichlid[-99,]
temp$Assignment <- cichlid.FP$Assignment
temp <- temp[,c(533,1:532)]
temp[,2:533][temp[,2:533] > 0] <- 1
cichlid.FP[,2:533][cichlid.FP[,2:533] > 0] <- 1
cichlid1 <- data.frame(colSums(temp[,2:533]))
cichlid2 <- data.frame(colSums(cichlid.FP[,2:533]))
cichlid.compare <- cbind(cichlid1, cichlid2)

## Calculate proportion of information lost
cichlid.compare$proportion <- cichlid.compare[,2]/cichlid.compare[,1]*100
cichlid.compare[is.na(cichlid.compare)] <- 0

## Examine mean and range of proportion of taxa retained
range(cichlid.compare$proportion)
mean(cichlid.compare$proportion[cichlid.compare$proportion!=0])

## This would result in up to 100% taxa being removed.
## On average, 0% species detections are retained.


########################################################
# OPTION 2: highest level of contamination in controls #
########################################################

## Store PCR negative and PCR positive controls in new dataframes
neg.controls <- true.assign[,which(!grepl("F|NU|P", colnames(true.assign)))]
head(neg.controls)

pos.controls <- true.assign[,which(!grepl("F|NU|N", colnames(true.assign)))]
head(pos.controls)

## Check positive controls for highest level of contamination (any DNA
## except cichlid). Remove column with taxonomic assignment.
pos <- pos.controls[,-1]

## Calculate total number of reads in samples
pos <- rbind(pos, colSums(pos))

## Create new dataframe containing the frequency of reads in each sample
pos.freq <- pos/c(pos[99,])

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
pos.freq$Threshold <- apply(pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(pos.freq[,110:115])

## Combine read frequencies with taxonomic assignment
species <- data.frame(pos.controls$Assignment)
colnames(species) <- "Assignment"
total <- data.frame("Total")
colnames(total) <- "Assignment"
species <- rbind(species, total)
pos.freq <- cbind(species, pos.freq)
rownames(pos.freq) <- NULL

## Print contamination threshold based on max. level of non-cichlid
## DNA contamination
## Exclude unassigned from threshold determination due to large number of
## reads, which cannot be distinguished as poorly amplified cichlid DNA 
## or other taxa
rownames(pos.freq) <- pos.freq$Assignment
pos.freq <- pos.freq[,-1]
max(pos.freq[-c(82,97,99),115])   # 0.8396

## In this scenario, any assignments <83.96% total reads in biological 
## samples would be considered contamination. This threshold cannot be
## applied to the data as it is too stringent and would effectively remove
## all taxonomic assignments. However, we will plot for comparison with other
## thresholds.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
pos.test <- cichlid.freq[,-534]
pos.test[pos.test <= 0.8396] <- 0

## Now convert back into read counts.
## Remove last row and column containing total frequency and thresholds.
## Add the total read counts to convert assignment frequencies back to 
## read counts for all samples.
pos.test <- pos.test[-99,]
rownames(pos.test) <- NULL

total.counts <- data.frame(colSums(cichlid[-99,]))
total.counts <- data.frame(t(total.counts))
colnames(total.counts) <- gsub('[.]', '-', colnames(total.counts))
total.counts <- cbind(total, total.counts)

## Now convert frequencies back to read counts
pos.conversion <- smartbind(pos.test, total.counts)
rownames(pos.conversion) <- pos.conversion$Assignment
pos.conversion <- pos.conversion[,-1]
pos.FP <- pos.conversion*c(pos.conversion[99,])

## Remove total row, reset row names and recreate Assignment column
pos.FP <- pos.FP[-99,]
pos.FP$Assignment <- rownames(pos.FP)
pos.FP <- pos.FP[,c(533,1:532)]
rownames(pos.FP) <- NULL

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- cichlid[-99,]
temp$Assignment <- pos.FP$Assignment
temp <- temp[,c(533,1:532)]
temp[,2:533][temp[,2:533] > 0] <- 1
pos.FP[,2:533][pos.FP[,2:533] > 0] <- 1
pos1 <- data.frame(colSums(temp[,2:533]))
pos2 <- data.frame(colSums(pos.FP[,2:533]))
pos.compare <- cbind(pos1, pos2)

## Calculate proportion of information lost
pos.compare$proportion <- pos.compare[,2]/pos.compare[,1]*100
pos.compare[is.na(pos.compare)] <- 0

## Examine mean and range of proportion of taxa retained
range(pos.compare$proportion)
mean(pos.compare$proportion)

## This would result in up to 100% taxa being removed.
## On average, 3.44% species detections are retained.


#######################################
# OPTION 3: Taxon-specific thresholds #
#######################################

## Check positive controls for highest level of contamination (any DNA
## except cichlid). Remove column with taxonomic assignment.
pos <- pos.controls[,-1]

## Calculate total number of reads in samples
pos <- rbind(pos, colSums(pos))

## Create new dataframe containing the frequency of reads in each sample
pos.freq <- pos/c(pos[99,])

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
pos.freq$Threshold <- apply(pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(pos.freq)

## Combine read frequencies with taxonomic assignment
species <- data.frame(pos.controls$Assignment)
colnames(species) <- "Assignment"
total <- data.frame("Total")
colnames(total) <- "Assignment"
species <- rbind(species, total)
pos.freq <- cbind(species, pos.freq)
rownames(pos.freq) <- NULL

## Manually change the species threshold value for cichlid to 0 as 
## this will be a true contaminant in the dataset and no false positive 
## threshold required.
pos.freq$Threshold[82] <- 0

## Check this manual edit has occurred
head(pos.freq[82,116])

## In this scenario, any assignments less than threshold in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Apply thresholds
SS.test <- cichlid.freq[,-534]
SS.test[SS.test <= pos.freq$Threshold] <- 0

## Now convert back into read counts.
## Remove last row and column containing total frequency and thresholds.
## Add the total read counts to convert assignment frequencies back to 
## read counts for all samples.
SS.test <- SS.test[-99,]
rownames(SS.test) <- NULL

total.counts <- data.frame(colSums(cichlid[-99,]))
total.counts <- data.frame(t(total.counts))
total.counts <- cbind(total, total.counts)

## Now convert frequencies back to read counts
SS.conversion <- smartbind(SS.test, total.counts)
rownames(SS.conversion) <- SS.conversion$Assignment
SS.conversion <- SS.conversion[,-1]
SS.FP <- SS.conversion*c(SS.conversion[99,])

## Remove total row, reset row names and recreate Assignment column
SS.FP <- SS.FP[-99,]
SS.FP$Assignment <- rownames(SS.FP)
SS.FP <- SS.FP[,c(533,1:532)]
rownames(SS.FP) <- NULL

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- cichlid[-99,]
temp$Assignment <- SS.FP$Assignment
temp <- temp[,c(533,1:532)]
temp[,2:533][temp[,2:533] > 0] <- 1
SS.FP[,2:533][SS.FP[,2:533] > 0] <- 1
SS1 <- data.frame(colSums(temp[,2:533]))
SS2 <- data.frame(colSums(SS.FP[,2:533]))
SS.compare <- cbind(SS1,SS2)

## Calculate proportion of information lost
SS.compare$proportion <- SS.compare[,2]/SS.compare[,1]*100
SS.compare[is.na(SS.compare)] <- 0

## Examine mean and range of proportion of taxa retained
range(SS.compare$proportion)
mean(SS.compare$proportion)

## This would result in up to 100% taxa being removed.
## On average, 61.45% species detections are retained.


###########
# SUMMARY # 
###########

## Tidy dataframes
cichlid.compare$Sample <- rownames(cichlid.compare)
cichlid.compare <- cichlid.compare[,c(4,1:3)]
rownames(cichlid.compare) <- NULL
colnames(cichlid.compare)[2:4] <- c("sp_richness_NT",
                                    "sp_richness_TA",
                                    "prop_TA")
cichlid.compare$type <- "Cichlid"

pos.compare$Sample <- rownames(pos.compare)
pos.compare <- pos.compare[,c(4,1:3)]
rownames(pos.compare) <- NULL
colnames(pos.compare)[2:4] <- c("sp_richness_NT",
                                "sp_richness_TA",
                                "prop_TA")
pos.compare$type <- "Positive"

SS.compare$Sample <- rownames(SS.compare)
SS.compare <- SS.compare[,c(4,1:3)]
rownames(SS.compare) <- NULL
colnames(SS.compare)[2:4] <- c("sp_richness_NT",
                               "sp_richness_TA",
                               "prop_TA")
SS.compare$type <- "Taxon-specific"

## Combine dataframes
threshold.df <- rbind(cichlid.compare, pos.compare, SS.compare)

## Plot number of species retained after thresholds applied
p1 <- ggplot(dat=threshold.df, aes(x=Sample, y=sp_richness_TA, fill=type))
p1 <- p1 + geom_bar(stat="identity", position=position_dodge())
p1 <- p1 + scale_y_continuous(limits=c(0,30))
p1 <- p1 + labs(x="Sample", y="Detections remaining after threshold application")
p1 <- p1 + theme_bw()
p1 <- p1 + theme(panel.background = element_rect(fill = 'white'),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.y = element_text(colour="black"),
                 legend.position = "none",
                 text = element_text(size=20))
p1 <- p1 + facet_grid(type ~ .)
p1

## Plot proportion of species retained after thresholds applied
p1 <- ggplot(dat=threshold.df, aes(x=Sample, y=prop_TA, fill=type))
p1 <- p1 + geom_bar(stat="identity", position=position_dodge())
p1 <- p1 + scale_y_continuous(limits=c(0,100))
p1 <- p1 + labs(x="Sample", 
                y="Detections remaining after threshold application (%)")
p1 <- p1 + theme_bw()
p1 <- p1 + theme(panel.background = element_rect(fill = 'white'),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 #axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.y = element_text(colour="black"),
                 legend.position = "none",
                 text = element_text(size=20))
p1 <- p1 + facet_grid(type ~ .)
p1

## Going forward the taxon-specific thresholds will be used as these
## retain the majority of biological information.
## Recreate dataframe with threshold applied.
SS.FP <- SS.conversion*c(SS.conversion[99,])
SS.FP <- SS.FP[-99,]
SS.FP$Assignment <- rownames(SS.FP)
SS.FP <- SS.FP[,c(533,1:532)]
rownames(SS.FP) <- NULL

## Remove taxonomic assignments with no read counts across any samples
SS.FP <- SS.FP[apply(SS.FP[,-1], 1, function(x) !all(x==0)),]


#################
# CONTAMINATION #
#################

## Examine how much contamination occured in negative controls (field blanks,
## filtration blanks, extraction blanks, and PCR negative controls)
contamination <- neg.controls

## Make Assignment column row names
rownames(contamination) <- contamination$Assignment
contamination <- contamination[,-1]

## Calculate total number of reads in each control
contamination <- rbind(contamination, colSums(contamination))

## Calculate frequency of reads in each control
contamination <- contamination/c(contamination[99,])

## Replace NA values with 0
contamination[is.na(contamination)] <- 0

## Remove total frequency row from dataframe
contamination <- contamination[-99,]

## Make row names first column in dataframe
contamination$Assignment <- rownames(contamination)
contamination <- contamination[,c(115,1:114)]
rownames(contamination) <- NULL

## Transpose dataframe for negative controls
contaminants <- setNames(data.frame(t(contamination[,-1])), contamination[,1])

## Make row names first column in dataframe
contaminants$ID <- rownames(contaminants)
contaminants <- contaminants[,c(99,1:98)]
rownames(contaminants) <- NULL

## Melt dataframe for plotting
contaminants <- melt(contaminants, id="ID")

## Rename columns
colnames(contaminants)[2:3] <- c("Assignment", "Frequency")

## Plot contamination found in negative controls
hm0 <- ggplot(contaminants, aes(x=ID, 
                                y=fct_rev(as_factor(Assignment)), 
                                fill=Frequency))
hm0 <- hm0 + geom_tile(colour="grey80")
hm0 <- hm0 + scale_fill_gradientn(name="Proportional\nread counts", 
                                  limits=c(0,1),
                                  breaks=c(0,0.25,0.50,0.75,1),
                                  labels=scales::number_format(accuracy=0.01,
                                                               decimal.mark="."),
                                  colours=c("white","red","black"), 
                                  values=c(0,0.1,1))
hm0 <- hm0 + labs(x="PCR negative controls", y="Taxonomic assignment")
hm0 <- hm0 + theme_bw()
hm0 <- hm0 + theme(panel.grid.major = element_line(colour="white"),
                   panel.grid.minor = element_line(colour="white"), 
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_text(colour="black"),
                   text = element_text(size=20),
                   legend.key.size = unit(1, 'lines'))
hm0



#' ---
#'
#' ## 4) Taxonomic resolution
#'

## Remove any assignments above species level, excluding the genus Anas
vert.dat <- SS.FP[which(grepl(" |Anas", SS.FP$Assignment)),]

## Change Anas to Anas spp.
vert.dat$Assignment <- gsub("Anas", "Anas spp.", vert.dat$Assignment)

## Reset row names of data frame for further indexing
rownames(vert.dat) <- NULL

## Remove PCR positive control.
## Remove human and domestic species as these cannot be ruled out as 
## environmental/laboratory contamination.
## Remove species likely to have resulted from DNA transport by waterfowl and 
## humans from water bodies in the wider catchment.
vert.dat <- vert.dat[which(!grepl("Bos|Canis|Equus|Felis|Homo|Meleagris|Numida|Ovis|Rhamphochromis|Sus",
                                  vert.dat$Assignment)),]

## Reset row names of data frame for further indexing
rownames(vert.dat) <- NULL

## Examine final dataset for downstream analyses
summary(vert.dat)
head(vert.dat)
names(vert.dat)
str(vert.dat)



#' ---
#' 
#' ## 5) Data analysis
#' 
#'
#' Examine the effect of vertebrate group species richness on probability of
#' GCN detection.
#' 
#' First, I need to create columns containing the number of species
#' in each vertebrate group found in each eDNA sample.
#' 
#' Convert the read count data to binary presence-absence data.
#' 

## Transpose dataframe
vert.binary <- setNames(data.frame(t(vert.dat[,-1])), vert.dat[,1])

## Make row names a column in data frame
vert.binary$Sample <- rownames(vert.binary)

## Move to start of dataframe
vert.binary <- vert.binary[,c(53,1:52)]

## Reset row names
rownames(vert.binary) <- NULL

## Convert sequence read counts to presence-absence
vert.binary[,2:53][vert.binary[,2:53] > 0] <- 1

## Replace space with underscore in species names
colnames(vert.binary) <- gsub(" ", "_", colnames(vert.binary))

## Total number of fish species detected in each sample
vert.binary$Fish <- rowSums(vert.binary[,c(3,6:7,10,13:14,16,20:21,34,38,
                                           42,45,49)])
head(vert.binary[,50:54])

## Total number of amphibian species excluding GCN detected in 
## each sample
vert.binary$Amphibian <- rowSums(vert.binary[,c(8,25:26,36,43)])
head(vert.binary[,50:55])

## Total number of waterfowl species detected in each sample
vert.binary$Waterfowl <- rowSums(vert.binary[,c(2,4,17:18,22)])
head(vert.binary[,50:56])

## Total number of terrestrial bird species detected in each sample
vert.binary$Terrestrial.bird <- rowSums(vert.binary[,c(9,11,15,19,23,37,
                                                       39,41,47,50:51)])
head(vert.binary[,55:57])

## Total up number of mammal species detected in each sample
vert.binary$Mammal <- rowSums(vert.binary[,c(5,12,24,27:33,35,40,44,46,
                                             48,53)])
head(vert.binary[,55:58])

## Plot number of vertebrate species present in ponds with and without GCN
p1 <- ggplot(vert.binary, 
             aes(x=factor(Amphibian),
                 fill=factor(Triturus_cristatus)))
p1 <- p1 + geom_bar(stat="count", 
                    alpha=0.7, 
                    position="fill", 
                    colour="black")
p1 <- p1 + geom_text(stat="count", 
                     aes(label=..count..),
                     size=5, position=position_fill(), vjust=1.5)
p1 <- p1 + scale_y_continuous(labels = percent_format(), 
                              expand = c(0,0))
p1 <- p1 + scale_fill_manual(name=expression(italic("T. cristatus")),
                             values=c("grey45","orange"),
                             breaks=c("0","1"),
                             labels=c("Negative",
                                      "Positive"))
p1 <- p1 + labs(title="(a)", 
                x="Number of other amphibian species", 
                y="Ponds")
p1 <- p1 + theme(panel.background = element_rect(fill = "white"),
                 axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                 axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=20))
p1


p2 <- ggplot(vert.binary, 
             aes(x=factor(Fish),
                 fill=factor(Triturus_cristatus)))
p2 <- p2 + geom_bar(stat="count", 
                    alpha=0.7, 
                    position="fill", 
                    colour="black")
p2 <- p2 + geom_text(stat="count", 
                     aes(label=..count..),
                     size=5, position=position_fill(), vjust=1.5)
p2 <- p2 + scale_y_continuous(labels = percent_format(), 
                              expand = c(0,0))
p2 <- p2 + scale_fill_manual(name=expression(italic("T. cristatus")),
                             values=c("grey45","orange"),
                             breaks=c("0","1"),
                             labels=c("Negative",
                                      "Positive"))
p2 <- p2 + labs(title="(b)", 
                x="Number of fish species", 
                y="Ponds")
p2 <- p2 + theme(panel.background = element_rect(fill = "white"),
                 axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                 axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=20))
p2


p3 <- ggplot(vert.binary, 
             aes(x=factor(Waterfowl),
                 fill=factor(Triturus_cristatus)))
p3 <- p3 + geom_bar(stat="count", 
                    alpha=0.7, 
                    position="fill", 
                    colour="black")
p3 <- p3 + geom_text(stat="count", 
                     aes(label=..count..),
                     size=5, position=position_fill(), vjust=1.5)
p3 <- p3 + scale_y_continuous(labels = percent_format(), 
                              expand = c(0,0))
p3 <- p3 + scale_fill_manual(name=expression(italic("T. cristatus")),
                             values=c("grey45","orange"),
                             breaks=c("0","1"),
                             labels=c("Negative",
                                      "Positive"))
p3 <- p3 + labs(title="(c)", 
                x="Number of waterfowl species", 
                y="Ponds")
p3 <- p3 + theme(panel.background = element_rect(fill = "white"),
                 axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                 axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=20))
p3


p4 <- ggplot(vert.binary, 
             aes(x=factor(Terrestrial.bird),
                 fill=factor(Triturus_cristatus)))
p4 <- p4 + geom_bar(stat="count", 
                    alpha=0.7, 
                    position="fill", 
                    colour="black")
p4 <- p4 + geom_text(stat="count", 
                     aes(label=..count..),
                     size=5, position=position_fill(), vjust=1.5)
p4 <- p4 + scale_y_continuous(labels = percent_format(), 
                              expand = c(0,0))
p4 <- p4 + scale_fill_manual(name=expression(italic("T. cristatus")),
                             values=c("grey45","orange"),
                             breaks=c("0","1"),
                             labels=c("Negative",
                                      "Positive"))
p4 <- p4 + labs(title="(d)", 
                x="Number of terrestrial bird species", 
                y="Ponds")
p4 <- p4 + theme(panel.background = element_rect(fill = "white"),
                 axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                 axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=20))
p4


p5 <- ggplot(vert.binary, 
             aes(x=factor(Mammal),
                 fill=factor(Triturus_cristatus)))
p5 <- p5 + geom_bar(stat="count", 
                    alpha=0.7, 
                    position="fill", 
                    colour="black")
p5 <- p5 + geom_text(stat="count", 
                     aes(label=..count..),
                     size=5, position=position_fill(), vjust=1.5)
p5 <- p5 + scale_y_continuous(labels = percent_format(), 
                              expand = c(0,0))
p5 <- p5 + scale_fill_manual(name=expression(italic("T. cristatus")),
                             values=c("grey45","orange"),
                             breaks=c("0","1"),
                             labels=c("Negative",
                                      "Positive"))
p5 <- p5 + labs(title="(e)", 
                x="Number of mammal species", 
                y="Ponds")
p5 <- p5 + theme(panel.background = element_rect(fill = "white"),
                 axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                 axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "bottom",
                 text = element_text(size=20))
p5


## Show all plots:
grid.arrange(arrangeGrob(p1,p2,p3,p4,p5, nrow=5, ncol=1))

## Import habitat metadata for eDNA samples and extract region column
metadata <- read.csv("../Data/HabitatMetadata.csv", header=TRUE)
region.dat <- metadata[,c(1,5)]

## Add County to presence/absence dataset and move to start of dataframe
vert.binary2 <- merge(vert.binary, region.dat, all=TRUE)
vert.binary2 <- vert.binary2[,c(1,59,2:58)]

## Add county information for ADAS ponds surveyed by environmental consultants
vert.binary2$County <- as.character(vert.binary2$County)
vert.binary2$County[c(140,203,321,393,505:532)] <- "Undisclosed"

## Make Sample and County factor columns
vert.binary2$Sample <- as.factor(vert.binary2$Sample)
vert.binary2$County <- as.factor(vert.binary2$County)

## Change column name from Sample to Pond and reset row names
colnames(vert.binary2)[1] <- "Pond"

## Statistically test relationships between vertebrate group richness and 
## GCN occupancy. Perform model with binary response variable:
richnessGCN <- glmer(Triturus_cristatus ~ (1|County) 
                     + Fish + Amphibian + Waterfowl + Terrestrial.bird 
                     + Mammal,
                     family = binomial,
                     data = vert.binary2)

## Null model:
nullGCN <- glmer(Triturus_cristatus ~ (1|County) + 1,
                 family = binomial,
                 data = vert.binary2)

## Compare richness model to null model:
anova(richnessGCN, nullGCN)

## LRT p-value significant and AIC of null model much higher thus richness of 
## vertebrate groups are appropriate explanatory variables

## Examine model:
summary(richnessGCN)
anova(richnessGCN)
drop1(richnessGCN, test = "Chi")  # for binomial/integer y models, 
                                  # statistics for significance of each term
display(richnessGCN)
se.ranef(richnessGCN)   # levels of random factor centred around 0

## summary() gives inflated value for model residual deviance so usual
## methods of calculating residual deviance are unreliable for GLMMs.
## Use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(richnessGCN)

## Model not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion

## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(richnessGCN, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(richnessGCN)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(richnessGCN)

## Plot the residuals against each x vairable
plot(sresid ~ vert.binary2$Fish)
plot(sresid ~ vert.binary2$Amphibian)
plot(sresid ~ vert.binary2$Waterfowl)
plot(sresid ~ vert.binary2$Terrestrial.bird)
plot(sresid ~ vert.binary2$Mammal)

## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(richnessGCN)
## 12.74% explained


#'
#' An additional test of model fit is the Hosmer and Lemeshow Goodness 
#' of Fit Test, which is similar to a chi-square test and indicates the 
#' extent to which the model provides better fit than a null model with 
#' no predictors.
#' 
#' If chi-square goodness of fit is not significant, then the model has 
#' adequate fit. Similarly, if the test is significant, the model does not 
#' adequately fit the data.
#' 

hoslem.test(vert.binary2$Triturus_cristatus, fitted(richnessGCN))

## Model does not fit that well as p-value is highly significant 
## i.e. significant difference between the model and the observed 
## data

## Make dataframe of p-values to calculate Bonferroni and Benjamini-Hochberg
## corrections to account for Type I error (i.e. false positives)
GR.results <- data.frame(Variable = c("Fish", "Amphibian", "Waterfowl",
                                      "Terrestrial bird", "Mammal"),
                         P = c(0.1799168, 0.0001648, 0.0002348, 0.8838733,
                               0.0035487))

## Calculate Bonferroni and Benjamini-Hochberg corrections
GR.results$Bonferroni <- p.adjust(GR.results$P, 
                                  method="bonferroni",
                                  n=length(GR.results$P))

GR.results$Benjamini_Hochberg <- p.adjust(GR.results$P, 
                                          method="BH", 
                                          n=length(GR.results$P))


## Obtain predicted values for the full data set.
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(richnessGCN, newdata=vert.binary2, 
                            se.fit = TRUE, print.matrix=T))

## Create new data set with fitted values and original data for plots 1 and 2
box.dat <- cbind(vert.binary2, fit) 

## Plot showing relationship between GCN presence and vertebrate group richness
## using predicted values from model:
p1m <- ggplot(box.dat) 
p1m <- p1m + geom_jitter(aes(x=factor(Amphibian), 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         height=0.3, width=0.3, cex=1)
p1m <- p1m + scale_colour_gradient(low="gray45", high="orange")
p1m <- p1m + geom_boxplot(aes(x=factor(Amphibian), 
                              y=fit), 
                          outlier.shape=NA, 
                          alpha=0.7)
p1m <- p1m + labs(title="", 
                  y = expression(paste("Probability of ", italic("T. cristatus"), " detection")), 
                  x = "Number of other amphibian species") 
p1m <- p1m + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=20))
p1m


p2m <- ggplot(box.dat) 
p2m <- p2m + geom_jitter(aes(x=factor(Fish), 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         height=0.3, width=0.3, cex=1)
p2m <- p2m + scale_colour_gradient(low="gray45", high="orange")
p2m <- p2m + geom_boxplot(aes(x=factor(Fish), 
                              y=fit), 
                          outlier.shape=NA, 
                          alpha=0.7)
p2m <- p2m + labs(title="", 
                  y = expression(paste("Probability of ", italic("T. cristatus"), " detection")), 
                  x = "Number of fish species") 
p2m <- p2m + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=20))
p2m


p3m <- ggplot(box.dat) 
p3m <- p3m + geom_jitter(aes(x=factor(Waterfowl), 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         height=0.3, width=0.3, cex=1)
p3m <- p3m + scale_colour_gradient(low="gray45", high="orange")
p3m <- p3m + geom_boxplot(aes(x=factor(Waterfowl), 
                              y=fit), 
                          outlier.shape=NA, 
                          alpha=0.7)
p3m <- p3m + labs(title="", 
                  y = expression(paste("Probability of ", italic("T. cristatus"), " detection")), 
                  x = "Number of waterfowl species") 
p3m <- p3m + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=20))
p3m


p4m <- ggplot(box.dat) 
p4m <- p4m + geom_jitter(aes(x=factor(Terrestrial.bird), 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         height=0.3, width=0.3, cex=1)
p4m <- p4m + scale_colour_gradient(low="gray45", high="orange")
p4m <- p4m + geom_boxplot(aes(x=factor(Terrestrial.bird), 
                              y=fit), 
                          outlier.shape=NA, 
                          alpha=0.7)
p4m <- p4m + labs(title="", 
                  y = expression(paste("Probability of ", italic("T. cristatus"), " detection")), 
                  x = "Number of terrestrial bird species") 
p4m <- p4m + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=20))
p4m


p5m <- ggplot(box.dat) 
p5m <- p5m + geom_jitter(aes(x=factor(Mammal), 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         height=0.3, width=0.3, cex=1)
p5m <- p5m + scale_colour_gradient(low="gray45", high="orange")
p5m <- p5m + geom_boxplot(aes(x=factor(Mammal), 
                              y=fit), 
                          outlier.shape=NA, 
                          alpha=0.7)
p5m <- p5m + labs(title="", 
                  y = expression(paste("Probability of ", italic("T. cristatus"), " detection")),
                  x = "Number of mammal species") 
p5m <- p5m + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   text = element_text(size=20),
                   legend.position="none",
                   legend.key.size = unit(1, "cm"))
p5m

## extract common legend
mylegend <- g_legend(p5)

## Show all plots with common legend
p1_5 <- grid.arrange(arrangeGrob(p1,p1m,
                                 p2,p2m,
                                 p3,p3m,
                                 p4,p4m,
                                 p5 + theme(legend.position="none"),p5m, 
                                 nrow=5, ncol=2),
                     mylegend, nrow=2, heights=c(10, 1))


#'
#' Now, examine individual species in GCN and non-GCN ponds. 
#' 
#' First, need to create a new dataframe.
#' 

## Total number of times an assignment occurred in GCN positive ponds
GCNpos <- subset(vert.binary, vert.binary$Triturus_cristatus == 1)
GCNneg <- subset(vert.binary, vert.binary$Triturus_cristatus == 0)

## Create dataframe containing species excluding GCN in GCN positive ponds
GCNpos.spp <- data.frame(GCNpos[,c(2:51,53)])

## Calculate number of GCN positive ponds each species found in
GCNpos.spp <- rbind(GCNpos.spp, colSums(GCNpos.spp))

## Create dataframe containing species excluding GCN in GCN negative ponds
GCNneg.spp <- data.frame(GCNneg[,c(2:51,53)])

## Calculate number of GCN negative ponds each species found in
GCNneg.spp <- rbind(GCNneg.spp, colSums(GCNneg.spp))

## Bind rows containing pond totals from each dataframe
totals <- data.frame(rbind(GCNpos.spp[149,], GCNneg.spp[385,]))
rownames(totals) <- c("P", "N")

## Transpose data frame
totals <- t(totals)

## Melt dataframe to create new factor column containing whether ponds are
## positive or negative for GCN
ponds.spp <- melt(totals, id=c("P","N"))
colnames(ponds.spp) <- c("Species", "Triturus_cristatus", "No_ponds")

## To complete the dataframe, we need to create another factor column
## indicating what vertebrate groups species belong to
ponds.spp$Group <- factor(ifelse(ponds.spp$Species %in% ponds.spp[c(2,5:6,9,12:13,15,19:20,33,37,41,44,48),1], "Fish",
                          ifelse(ponds.spp$Species %in% ponds.spp[c(7,24:25,35,42),1], "Amphibian",
                          ifelse(ponds.spp$Species %in% ponds.spp[c(1,3,8,10,14,16:18,21:22,36,38,40,46,49:50),1], "Bird",
                          ifelse(ponds.spp$Species %in% ponds.spp[c(4,10:11,23,26:32,34,39,43,45,47,51),1], "Mammal",
                                 "NA")))))

## Remove underscores from species names
ponds.spp$Species <- gsub("_", " ", ponds.spp$Species)

## Subset data by vertebrate groups
Fish <- subset(ponds.spp, Group == "Fish")
Amphibian <- subset(ponds.spp, Group == "Amphibian")
Bird <- subset(ponds.spp, Group == "Bird")
Mammal <- subset(ponds.spp, Group == "Mammal")

## Plot proportion of amphibian species in GCN positive and negative ponds
p6 <- ggplot(Amphibian[order(Amphibian$Triturus_cristatus, decreasing=T),],
             aes(x=Species, 
                 y=No_ponds, 
                 fill=factor(Triturus_cristatus, levels=c("N","P"))))
p6 <- p6 + geom_bar(stat="identity", 
                    alpha=0.7, 
                    position="fill", 
                    colour="black")
p6 <- p6 + geom_text(data=subset(Amphibian, No_ponds > 0), 
                     aes(label=No_ponds), 
                     size=5, 
                     position=position_fill(), 
                     vjust=1.5)
p6 <- p6 + scale_y_continuous(labels = percent_format(), 
                              expand = c(0,0))
p6 <- p6 + scale_fill_manual(name=expression(italic("T. cristatus")),
                             values=c("grey45","orange"),
                             breaks=c("P","N"),
                             labels=c("Positive","Negative"))
p6 <- p6 + labs(title="(a)", 
                x="Other amphibian species", 
                y="Ponds")
p6 <- p6 + theme(panel.background = element_rect(fill = "white"),
                 axis.line.x = element_line(colour = "black", 
                                            size=0.5, 
                                            linetype="solid"),
                 axis.line.y = element_line(colour = "black", 
                                            size=0.5, 
                                            linetype="solid"),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black", 
                                            angle = 45, 
                                            vjust = 1, 
                                            hjust=1, 
                                            face="italic"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", 
                                           hjust=0),
                 text = element_text(size=24),
                 legend.position = "bottom")
p6


## Plot proportion of fish species in GCN positive and negative ponds
p7 <- ggplot(Fish[order(Fish$Triturus_cristatus, decreasing=T),],
             aes(x=Species, 
                 y=No_ponds, 
                 fill=factor(Triturus_cristatus, levels=c("N","P"))))
p7 <- p7 + geom_bar(stat="identity", 
                    alpha=0.7, 
                    position="fill", 
                    colour="black")
p7 <- p7 + geom_text(data=subset(Fish, No_ponds > 0), 
                     aes(label=No_ponds), 
                     size=5, 
                     position=position_fill(), 
                     vjust=1.5)
p7 <- p7 + scale_y_continuous(labels = percent_format(), 
                              expand = c(0,0))
p7 <- p7 + scale_fill_manual(name="Great crested newt",
                             values=c("grey45","orange"),
                             breaks=c("P","N"),
                             labels=c("Positive","Negative"))
p7 <- p7 + labs(title="(b)", 
                x="Fish species", 
                y="Ponds")
p7 <- p7 + theme(panel.background = element_rect(fill = "white"),
                 axis.line.x = element_line(colour = "black", 
                                            size=0.5, 
                                            linetype="solid"),
                 axis.line.y = element_line(colour = "black", 
                                            size=0.5, 
                                            linetype="solid"),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black", 
                                            angle = 45, 
                                            vjust = 1, 
                                            hjust=1, 
                                            face="italic"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", 
                                           hjust=0),
                 text = element_text(size=24),
                 legend.position = "none")
p7


## Plot proportion of bird species in GCN positive and negative ponds
p8 <- ggplot(Bird[order(Bird$Triturus_cristatus, decreasing=T),],
             aes(x=Species, 
                 y=No_ponds, 
                 fill=factor(Triturus_cristatus, levels=c("N","P"))))
p8 <- p8 + geom_bar(stat="identity", 
                    alpha=0.7, 
                    position="fill", 
                    colour="black")
p8 <- p8 + geom_text(data=subset(Bird, No_ponds > 0), 
                     aes(label=No_ponds), 
                     size=5, 
                     position=position_fill(), 
                     vjust=1.5)
p8 <- p8 + scale_y_continuous(labels = percent_format(), 
                              expand = c(0,0))
p8 <- p8 + scale_fill_manual(name="Great crested newt",
                             values=c("grey45","orange"),
                             breaks=c("P","N"),
                             labels=c("Positive","Negative"))
p8 <- p8 + labs(title="(c)", 
                x="Bird species", 
                y="Ponds")
p8 <- p8 + theme(panel.background = element_rect(fill = "white"),
                 axis.line.x = element_line(colour = "black", 
                                            size=0.5, 
                                            linetype="solid"),
                 axis.line.y = element_line(colour = "black", 
                                            size=0.5, 
                                            linetype="solid"),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black", 
                                            angle = 45, 
                                            vjust = 1, 
                                            hjust=1, 
                                            face="italic"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", 
                                           hjust=0),
                 text = element_text(size=24),
                 legend.position = "none")
p8


## Plot mammal species in GCN and non-GCN ponds
p9 <- ggplot(Mammal[order(Mammal$Triturus_cristatus, decreasing=T),],
             aes(x=Species, 
                 y=No_ponds, 
                 fill=factor(Triturus_cristatus, levels=c("N","P"))))
p9 <- p9 + geom_bar(stat="identity", 
                    alpha=0.7, 
                    position="fill", 
                    colour="black")
p9 <- p9 + geom_text(data=subset(Mammal, No_ponds > 0), 
                     aes(label=No_ponds), 
                     size=5, 
                     position=position_fill(), 
                     vjust=1.5)
p9 <- p9 + scale_y_continuous(labels = percent_format(), 
                              expand = c(0,0))
p9 <- p9 + scale_fill_manual(name="Great crested newt",
                             values=c("grey45","orange"),
                             breaks=c("P","N"),
                             labels=c("Positive","Negative"))
p9 <- p9 + labs(title="(d)", 
                x="Mammal species", 
                y="Ponds")
p9 <- p9 + theme(panel.background = element_rect(fill = "white"),
                 axis.line.x = element_line(colour = "black", 
                                            size=0.5, 
                                            linetype="solid"),
                 axis.line.y = element_line(colour = "black", 
                                            size=0.5, 
                                            linetype="solid"),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black", 
                                            angle = 45, 
                                            vjust = 1, 
                                            hjust=1, 
                                            face="italic"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", 
                                           hjust=0),
                 text = element_text(size=24),
                 legend.position = "none")
p9

## Show all plots:
grid_arrange_shared_legend(p6,p7,p8,p9, nrow=2, ncol=2)



#' ---
#'
#' ## Are these interactions significant?
#' 
#' The probabilistic model of species co-occurrence (Veech 2013) 
#' measures co-occurrence in the most straightforward way as the 
#' number of sampling sites where two species co-occur. Observed 
#' co-occurrence can be compared to the expected co-occurrence where 
#' the latter is the product of the two species probability of 
#' occurrence multiplied by the number of sampling sites.
#' 
#' This probabilistic model was developed in R to increase 
#' the availability of the model as an easy-to-use method for conducting 
#' pairwise co-occurrence analyses.
#' 
#' The package calculates the probability of selecting a site (or
#' sample) that has species #1 given that it already has species #2.
#' Subsequently, coocur() uses a hypergeometric approach
#' 

## Load package cooccur for probabalistic species coocurrence analysis 
## in R
library("cooccur")

## Adapt presence-absence dataframe to contain only species, no metadata
ponds <- vert.binary[,-c(54:58)]

## Transpose data frame so that species are first column
ponds <- setNames(data.frame(t(ponds[,-1])), ponds[,1])

## Remove underscore from species names
rownames(ponds) <- gsub("_", " ", rownames(ponds))

#'
#' The function accepts community data (e.g., species by site matrix 
#' or vice-versa) in the form of a data frame or matrix and returns 
#' a list containing pairwise species co-occurrence results.
#' 
#' The aim of the sample analysis is to determine the degree to 
#' which communities contain species that are positively, negatively, 
#' and randomly associated with one another, investigate the 
#' contribution of individual species to these patterns, and to 
#' quantify the strength of the positive and negative associations 
#' between species pairs.
#' 

## Run analysis
cooccur.ponds <- cooccur(mat = ponds,        # the dataframe
                         type = "spp_site",  # data organized with species as rows and sites as column
                         thresh = TRUE,      # summarize the most important species associations
                         spp_names = TRUE)   # use species names in data

## NB: thresh=TRUE removes species that simply do not have sufficient 
## occurrence data (i.e. they are expected to share less than one site)

## The cooccur() function produces an output object of class cooccur 
## containing all of the results from the co-occurrence analysis
class(cooccur.ponds)

## Calling summary() produces a readout of the total positive, 
## negative, and random species pairs classified by the algorithm.
summary(cooccur.ponds)

## 31 positive pairs
## 17 negative pairs
## 195 random pairs

#'
#' Of 1326 species pair combinations, 1083 pairs (81.67 %) were removed from 
#' the analysis because expected co-occurrence was < 1 and 243 pairs were 
#' analysed.
#' 
#' There is sufficient statistical power to analyse all pairs as none 
#' were unclassifiable.
#'  

## To obtain the complete set of species pairs analyzed, use prob.table()
prob.table(cooccur.ponds)

## sp1 = Order of species 1 in matrix
## sp2 = Order of species 2 in matric
## sp1_inc = Number of sites (or samples) that have species 1
## sp2_inc = Number of sites that have species 2
## obs_cooccur = Observed number of sites having both species
## prob_cooccur = Probability that both species occur at a site
## exp_cooccur = Expected number of sites having both species
## p_lt = Probability that two species co-occur at a frequency less 
## than observed number of co-occurrence sites if the two species 
## were distributed randomly (independently) of one another
## p_gt = Probability of co-occurrence at a frequency greater than 
## the observed frequency
## sp1_name = supplied name of sp1
## sp2_name = supplied name of sp2

## A list of only significant species combinations can be obtained 
## using the print method
print(cooccur.ponds)

#'
#' For a given species pair, p_lt and p_gt represent the probabilities 
#' that those species could co-occur less than or greater than what 
#' is observed in the data, respectively.
#' 
#' They can be interpreted as p-values, thus indicating significance 
#' levels for negative and positive co-occurrence patterns.
#' 
#' In the ponds data, there was 17 significant negative co-occurrence
#' pattern and 48 significant positive co-occurrence patterns.
#' 
#' However, as there are so many species pair combinations, there 
#' is risk of Type 1 error (i.e. false positives). To control for this,
#' we will adjust the p-values. The most typical correction for Type 1
#' error is Bonferroni, but it is vulnerable to Type 2 error (i.e. false
#' negatives). Therefore, we will also calculate the Benjamini-Hochberg 
#' correction as it is a sequential modified Bonferroni correction and
#' is less stringent.
#' 
#' 

## Make prob.table a dataframe
cooccur.results <- data.frame(prob.table(cooccur.ponds))

## Remove redundant columns and reorder dataframe
cooccur.results <- cooccur.results[,-c(1:7)]
cooccur.results <- cooccur.results[,c(3:4,1:2)]

## Calculate Bonferroni and Benjamini-Hochberg corrections for negative
## and positive species association p-values
cooccur.results$neg_Bonferroni <- p.adjust(cooccur.results$p_lt, 
                                           method="bonferroni", 
                                           n=length(cooccur.results$p_lt))

cooccur.results$neg_Benjamini_Hochberg <- p.adjust(cooccur.results$p_lt, 
                                                   method="BH",
                                                   n=length(cooccur.results$p_lt))

cooccur.results$pos_Bonferroni <- p.adjust(cooccur.results$p_gt, 
                                           method="bonferroni", 
                                           n=length(cooccur.results$p_gt))

cooccur.results$pos_Benjamini_Hochberg <- p.adjust(cooccur.results$p_gt, 
                                                   method="BH",
                                                   n=length(cooccur.results$p_gt))

## Reorder dataframe
cooccur.results <- cooccur.results[,c(1:3,5:6,4,7:8)]

## Subset data for species associations involving GCN
cooccur.GCN <- subset(cooccur.results, 
                      sp1_name == "Triturus cristatus" | sp2_name == "Triturus cristatus")

## Keep only the associations that were identified as significant by coocccur
## package
cooccur.GCN <- subset(cooccur.GCN,
                      p_lt < 0.05 | p_gt < 0.05)

#'
#' Use the plot() method on the results object. This will produce a 
#' visualization of all of the pairwise combinations of species and 
#' their co-occurrence signs (positive or negative) using a ggplot2 
#' heatmap.
#' 
#' The plot trims out any species that do not have any significant 
#' negative or positive associations and orders the remaining species 
#' starting from those with the most negative interactions to those 
#' with the most positive interactions.
#' 

plot(cooccur.ponds)

#'
#' The pair() function can be used to specifiy a specific species, 
#' by name or number, to inspect, by default only significant results 
#' are shown, but if all = TRUE then all results will be shown
#' 

## Inspect great crested newt
pair(mod = cooccur.ponds, "Triturus cristatus")
pair(mod = cooccur.ponds, "Triturus cristatus", all = TRUE)

#'
#' To understand each species individual contribution to the 
#' positive and negative species associations, need to create a 
#' pairing profile.
#' 
#' pair.attributes() produces a table of the percentage of each 
#' species total pairings that were classified as positive, negative,
#' and random (columns with prefix num are counts).
#' 

pair.attributes(cooccur.ponds)

#'
#' The primary goal is to compare the numbers of positive versus 
#' negative associations (unclassifiable pairings treated as random).
#' Same results can be visualized across all species by using the
#' function pair.profile() to create a box plot of these percentages.
#' This plot will show the percent of species pairs that were positive, 
#' negative, and random for all species.
#' 

p10 <- pair.profile(cooccur.ponds)
p10

#'
#' This communicates whether or not species tend to have mostly negative
#' or mostly positive interactions and suggests whether interactions 
#' are evenly distributed among species as opposed to being clustered
#' in a few species.
#' 
#' Effect sizes can also be calculated from co-occurrence analyses; 
#' they allow for comparisons among studies and methods as well as 
#' providing a quantitative measurement of co-occurrence for use in 
#' downstream analyses.
#' 
#' Effect size = difference between expected and observed frequency 
#' of co-occurrence
#' 
#' Values are standardized by dividing these differences by the 
#' number of sampling sites in the dataset. In standardized form, 
#' these values are bounded from -1 to 1, with positive values 
#' indicating positive associations and negative values indicating 
#' negative associations.
#' 

effects.table <- cooccur(mat = ponds, 
                         type = "spp_site", 
                         thresh = FALSE,
                         spp_names = TRUE, 
                         only_effects = TRUE, 
                         eff_standard = TRUE, 
                         eff_matrix = TRUE)

## Convert dist object produced by cooccur() to a matrix
effects.table <- as.matrix(dist(effects.table)) 

## Write matrix of effect sizes as csv file to working directory
write.csv(effects.table, "../Data/EffectSize_new.csv")

#'
#' Cohen provided rules of thumb for interpreting these effect sizes, 
#' suggesting that an r of .1 represents a 'small' effect size, .3 
#' represents a 'medium' effect size and .5 represents a 'large' 
#' effect size. 
#' 
#' To inspect the degree to pond species pairs deviate from their 
#' expected co-occurrence levels, plot the observed values against 
#' the expected value as a visual diagnostic. 
#' 

p11 <- obs.v.exp(cooccur.ponds)
p11

#' 
#' The probability calculations are based on the number of sites and 
#' the individual frequencies of occurrence and co-occurrence for 
#' each species pair. Therefore the conditions determining 
#' statistical power change with sample size and it is valuable to 
#' examine effect sizes for species pairs regardless of statistical 
#' significance. A detailed discussion of power, Type I and II error 
#' rates, and a comparison with other methods can be found in Veech 
#' (2013).
#' 
#' 
#' There are several species pairs that exhibit fewer and greater 
#' expected co-occurrences. These pairs are largely clustered towards 
#' having low or high expected co-occurrences in the first place.
#' 

## Plot pairwise profile and observed-expected plot alongside each 
## other
grid.arrange(p10, p11, ncol = 2,
             top = textGrob("Vertebrate cooccurrence at the UK pondscape", 
                            gp = gpar(cex = 2), 
                            vjust = 0.75))



#' ---
#' 
#' ## Biotic and abiotic determinants of GCN occupancy
#' 

## Adapt presence-absence dataframe and remove columns created for previous 
## analysis
vert.spp <- vert.binary[,1:53]

## Create column containing the total number of species detected in
## each sample i.e. species richness
vert.spp$Sp.richness <- rowSums(vert.spp[,2:53])

## Merge habitat metadata with the metabarcoding data using the Sample column 
habitat.dat <- merge(metadata, vert.spp, all=TRUE)

## Check merge was successful
head(habitat.dat[,1:6])
head(habitat.dat[,85:89])

## Remove samples which did not have associated pond metadata
## (ADAS UK Ltd samples and 4 samples from FERA)
habitat.dat <- habitat.dat[-c(505:532),]

## Change column name from Sample to Pond
colnames(habitat.dat)[1] <- "Pond"

## Now, check dataframe structure
summary(habitat.dat)
names(habitat.dat)
head(habitat.dat)
str(habitat.dat)

## corvif function and pairplot with the Pearson correlation 
## coefficients were imported earlier

## Matrix of Spearman's rank correlations for dataset
## Use Spearman's rank as no assumptions are made about linearity in a
## relationship between two variables
cor(habitat.dat[,10:19], method = "spearman")

## Booth et al. (1994) and Zuur et al. (2009) suggest that correlations
## between pairs of variables with magnitudes ± 0.5 indicate high
## collinearity

## View pairwise plots of variables
pairs(habitat.dat[,10:19])

## Below produces a pairwise plot but the strength of colinearity 
## strength is indicated by the size of the text
plot(habitat.dat[,10:19], 
     lower.panel = panel.cor, 
     upper.panel = panel.smooth2)


#' 
#' There appears to be collinearity between pond circumference, pond
#' length, pond width, and pond area.
#' 
#' Pond area encompasses length and width thus accounting for the same variance 
#' in the data as these variables. Therefore, it is logical to remove pond 
#' circumference, pond length and pond width to counter high collinearity 
#' between explanatory variables.
#' 
#' Shading (% total pond margin shaded) and terrestrial overhang (% pond
#' overhung by trees and shrubs) also appear to be collinear. Shading is a 
#' known driver of biodiversity in ponds (Sayer et al. 2012), thus this will 
#' be retained as an explanatory variable.
#' 
#' Similarly, there is no need to include overall habitat score and HSI
#' band as these are encompassed by HSI score.
#' 
#' 
#' Now check the variance inflation factors (VIFs) of variables to 
#' assess the extent of remaining collinearity.
#' 
#' This is done using a binomial GLM and logit link function containing
#' all explanatory variables to the GCN presence/absence data, then 
#' calculating the VIFs for each variable from the resulting model.
#' 

glmGCN <- glm(Triturus_cristatus ~ Max_Depth + Pond_Area 
              + Pond_Density + Shading + Macrophytes + HSI 
              + Bufo_bufo + Lissotriton_vulgaris 
              + Cyprinus_carpio + Pungitius_pungitius + Gasterosteus_aculeatus 
              + Fulica_atra + Gallinula_chloropus 
              + Permanance + Quality + Pond_Substrate + Inflow 
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals,
              family = binomial,
              data = habitat.dat)

## Use vif() function in car package to calculate VIFs
vif(glmGCN)

## All VIF values below 5, except HSI score which is likely higher as it
## is calculated from the other variables.
## Check possible relationships between response variable and nominal
## variables (factors) in a design plot.
factors <- data.frame(Triturus_cristatus = habitat.dat$Triturus_cristatus,
                      Permanance = habitat.dat$Permanance,
                      Quality = habitat.dat$Quality,
                      Pond_Substrate = habitat.dat$Pond_Substrate,
                      Inflow = habitat.dat$Inflow,
                      Outflow = habitat.dat$Outflow,
                      Pollution = habitat.dat$Pollution,
                      Amphibians = habitat.dat$Amphibians,
                      Waterfowl = habitat.dat$Waterfowl,
                      Fish = habitat.dat$Fish,
                      Woodland = habitat.dat$Woodland,
                      Rough_Grass = habitat.dat$Rough_Grass,
                      Scrub_Hedge = habitat.dat$Scrub_Hedge,
                      Ruderals = habitat.dat$Ruderals)

plot.design(Triturus_cristatus ~ Permanance + Quality + Pond_Substrate 
            + Inflow + Outflow + Pollution + Amphibians + Waterfowl 
            + Fish + Woodland + Rough_Grass + Scrub_Hedge + Ruderals, 
            data = factors, axes = T, xtick = T)

#' 
#' The highest mean values were observed with pond substrate,
#' amphibian presence and permanance.
#' 

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables
f1 <- formula(Triturus_cristatus ~ Max_Depth + Pond_Area 
              + Pond_Density + Shading + Macrophytes 
              + Bufo_bufo + Lissotriton_vulgaris 
              + Cyprinus_carpio + Pungitius_pungitius + Gasterosteus_aculeatus 
              + Fulica_atra + Gallinula_chloropus 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals)

GCN_tree <- rpart(f1, 
                  data = habitat.dat, 
                  method = "class",
                  minsplit=5, 
                  cp = 0.001)

par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(GCN_tree,uniform=TRUE, margin=0.1)
text(GCN_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - L. vulgaris presence
#' - fish presence
#' - toad presence
#' - amphibian presence
#' - pond area
#' - common moorhen presence
#' - pond substrate
#' - water quality
#' - pond density
#' - woodland
#' - permanance
#' - max depth
#' - outflow
#' - inflow
#' - scrub/hedge
#' - macrophytes
#' - shading
#' - ruderals
#' - waterfowl presence
#'  

## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(GCN_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 4 is optimal i.e. 4 explanatory variables. 
#' The best tree will be subjective.
#' 
#' L. vulgaris and T. cristatus known to have shared ecology including
#' foraging and habitat preference.
#' 
#' Large predatory fish species feed on newt larvae. Others increase turbidity
#' and remove macrophytes as well as invertebrate prey.
#' 
#' Presence of other amphibians could play positive or negative role
#' in GCN occupancy as prey availability will be increased but 
#' predation and competition for resources will also increase.
#' 
#' Larger ponds are likely to have more fish species but may also 
#' support greater GCN populations.
#' 
#' Waterfowl species may have shared habitat/prey preference as GCN.
#' 
#' Pond substrate and water quality (pollution, run-off) will influence 
#' macrophyte composition and macroinvertebrate prey availability.
#' 
#' Pond density can influence terrestrial migration and population viability.
#' 
#' GCN are known to spend majority of their life on land so terrestrial
#' habitat will influence pond occupancy, e.g. shading and grassland.
#' 
#' Permanant ponds will be more suitable for GCN breeding but may
#' contain more fish species. 
#' 
#' Depth may play a role in GCN occupancy as more fish species are 
#' likely to inhabit deep ponds and less breeding substrate available
#' for GCN.
#' 
#' Inflow will facilitate fish movement and agricultural runoff whereas
#' outflow may indicate whether pollutants and excess nutriets are 
#' flushed out.
#' 
#' Macrophyte cover is a known driver of biodiversity in ponds and
#' influence GCN occupancy.
#'  
#' Shading is an indication of terrestrial plant life encroaching 
#' on pond.
#' 
#' Great crested newt are dependent upon rough grass habitat for terrestrial
#' migration as well as scrub/hedge and ruderal habitat for hibernacula.
#'
#'  
#' Variables not included in the classification tree were:
#' 
#' - presence of other species that have interaction with GCN
#' - pollution
#' 
#' 
#' Therefore, variables to be taken forward in model selection are:
#' 
#' - L. vulgaris presence
#' - fish presence
#' - toad presence
#' - amphibian presence
#' - pond area
#' - common moorhen presence
#' - pond substrate
#' - water quality
#' - pond density
#' - woodland
#' - permanance
#' - max depth
#' - outflow
#' - inflow
#' - scrub/hedge
#' - macrophytes
#' - shading
#' - ruderals
#' - waterfowl presence
#' 
#' Although not identified by the classification tree, we will include 
#' presence of common carp, three-spined stickleback and ninespine 
#' stickleback as these fish predate GCN. We will also include common coot 
#' presence as co-occur identified a significant positive association.
#'
#' Several variables occur more than once in the tree which indicates 
#' weak non-linear relationships with the response variable thus a 
#' GAM may be more appropriate to model the data.
#' 
#' A GAM with binomial distribution and logistic link function will be
#' used to relate GCN presence/absence with the explanatory variables.
#' It is actually a Bernoulli distribution thus overdispersion cannot 
#' occur which is a major problem in ecological distribution studies.
#' 
#' Forward selection will be applied to find the optimal set of 
#' explanatory variables based on model AIC values.

## Use package mgcv for gam modelling
GCN.gam1 <- gam(Triturus_cristatus ~ Lissotriton_vulgaris, data=habitat.dat, family=binomial)
GCN.gam2 <- gam(Triturus_cristatus ~ Fish, data=habitat.dat, family=binomial)
GCN.gam3 <- gam(Triturus_cristatus ~ Bufo_bufo, data=habitat.dat, family=binomial)
GCN.gam4 <- gam(Triturus_cristatus ~ Amphibians, data=habitat.dat, family=binomial)
GCN.gam5 <- gam(Triturus_cristatus ~ s(Pond_Area), data=habitat.dat, family=binomial)
GCN.gam6 <- gam(Triturus_cristatus ~ Gallinula_chloropus, data=habitat.dat, family=binomial)
GCN.gam7 <- gam(Triturus_cristatus ~ Pond_Substrate, data=habitat.dat, family=binomial)
GCN.gam8 <- gam(Triturus_cristatus ~ Quality, data=habitat.dat, family=binomial)
GCN.gam9 <- gam(Triturus_cristatus ~ s(Pond_Density), data=habitat.dat, family=binomial)
GCN.gam10 <- gam(Triturus_cristatus ~ Woodland, data=habitat.dat, family=binomial)
GCN.gam11 <- gam(Triturus_cristatus ~ Permanance, data=habitat.dat, family=binomial)
GCN.gam12 <- gam(Triturus_cristatus ~ s(Max_Depth), data=habitat.dat, family=binomial)
GCN.gam13 <- gam(Triturus_cristatus ~ Outflow, data=habitat.dat, family=binomial)
GCN.gam14 <- gam(Triturus_cristatus ~ Inflow, data=habitat.dat, family=binomial)
GCN.gam15 <- gam(Triturus_cristatus ~ Scrub_Hedge, data=habitat.dat, family=binomial)
GCN.gam16 <- gam(Triturus_cristatus ~ s(Macrophytes), data=habitat.dat, family=binomial)
GCN.gam17 <- gam(Triturus_cristatus ~ s(Shading), data=habitat.dat, family=binomial)
GCN.gam18 <- gam(Triturus_cristatus ~ Ruderals, data=habitat.dat, family=binomial)
GCN.gam19 <- gam(Triturus_cristatus ~ Waterfowl, data=habitat.dat, family=binomial)
GCN.gam20 <- gam(Triturus_cristatus ~ Cyprinus_carpio, data=habitat.dat, family=binomial)
GCN.gam21 <- gam(Triturus_cristatus ~ Gasterosteus_aculeatus, data=habitat.dat, family=binomial)
GCN.gam22 <- gam(Triturus_cristatus ~ Pungitius_pungitius, data=habitat.dat, family=binomial)
GCN.gam23 <- gam(Triturus_cristatus ~ Fulica_atra, data=habitat.dat, family=binomial)

AIC(GCN.gam1,GCN.gam2,GCN.gam3,GCN.gam4,GCN.gam5,GCN.gam6, GCN.gam7,
    GCN.gam8,GCN.gam9,GCN.gam10,GCN.gam11,GCN.gam12,GCN.gam13,
    GCN.gam14,GCN.gam15,GCN.gam16,GCN.gam17,GCN.gam18,GCN.gam19,
    GCN.gam20,GCN.gam21,GCN.gam22,GCN.gam23)

## Smooth newt identified as best explanatory variable and most optimal model
summary(GCN.gam1)
summary(GCN.gam2)
summary(GCN.gam3)
summary(GCN.gam4)
summary(GCN.gam5)
summary(GCN.gam6)
summary(GCN.gam7)
summary(GCN.gam8)
summary(GCN.gam9)
summary(GCN.gam10)
summary(GCN.gam11)
summary(GCN.gam12)
summary(GCN.gam13)
summary(GCN.gam14)
summary(GCN.gam15)
summary(GCN.gam16)
summary(GCN.gam17)
summary(GCN.gam18)
summary(GCN.gam19)
summary(GCN.gam20)
summary(GCN.gam21)
summary(GCN.gam22)
summary(GCN.gam23)


#' 
#' Other variables are biologically important for GCN and several have an 
#' estimated 1 degree of freedom for the smoothers which is equivalent to a 
#' linear relationship. Therefore, GLM should be applied instead of GAM. 
#' GLM is also a parametric method. 
#' 
#' Spatial autocorrrelation is also common in ecological studies of species 
#' presence/absence as sites that are located within an animal's ranging 
#' capability are likely to be inhabited. In the case of GCN, individuals 
#' can migrate distances of 2km when they have finished breeding. Therefore, 
#' occurrence of GCN is likely in ponds that are closely located to one 
#' another in a given area. Furthermore, if pond distances are substantially 
#' smaller than the spatial extent of the study area, this can also lead to 
#' spatial auto-correlation. The underlying spatial patterns of habitat 
#' surrounding ponds could also be a source of spatial autocorrelation but the 
#' explanatory variables which describe the habitat should account for this.
#' 
#' One way to assess spatial autocorrelation is through correlograms 
#' of the data. These are graphical representations of spatial 
#' correlation between locations at a range of lag distances. Positive 
#' spatial correlation indicates spatial autocorrelation between data 
#' points. Negative spatial correlation may also indicate problems but 
#' is very uncommon.
#' 
#' Spline correlograms are smoothed correlograms using a spline 
#' function (Zuur *et al.* 2009).
#' 

## The ncf package will be used to produce spline correlograms
## Location data is required 
Correlog <- spline.correlog(x=habitat.dat$Easting,
                            y=habitat.dat$Northing,
                            z=habitat.dat$Triturus_cristatus, 
                            xmax=10000)
plot(Correlog)

#' 
#' A spline correlogram with 95% pointwise bootstrap CIs and maximum 
#' lag distance of 10km has been produced.
#' 
#' Some positive spatial auto-correlation is present but only
#' at short lag distances of less than or around 1km. This suggests
#' that spatial auto-correlation may be a problem for sites that are
#' located near to one another.
#' 
#' Check for spatial auto-correlation in Pearson residuals of logistic
#' regression model containing all explanatory variables fitted to the
#' presence/absence data.
#' 

glmGCN <- glm(Triturus_cristatus ~ Lissotriton_vulgaris + Fish + Bufo_bufo
              + Amphibians + Pond_Area + Gallinula_chloropus + Pond_Substrate
              + Quality + Pond_Density + Woodland + Permanance + Max_Depth
              + Outflow + Inflow + Scrub_Hedge + Macrophytes + Shading
              + Ruderals + Waterfowl + Cyprinus_carpio + Gasterosteus_aculeatus
              + Pungitius_pungitius + Fulica_atra,
              family = binomial,
              data = habitat.dat)

Correlog_glmGCN <- spline.correlog(x=habitat.dat$Easting,
                                   y=habitat.dat$Northing,
                                   z=residuals(glmGCN, type="pearson"),
                                   xmax=10000)
plot(Correlog_glmGCN)


#' 
#' A GLMM can account for dependencies within sites (i.e. counties) and is 
#' appropriate where nested model selection will be used and the data are 
#' structured hierarchically. Ponds are nested within counties thus a mixed 
#' model is necessary to account for spatial dependencies and county must be 
#' treated as a random effect.
#' 
#' Apply GLMM to variables most likely to influence GCN presence/absence
#' based on exploratory analysis. Check that a GLMM will adequately account 
#' for spatial auto-correlation that is present by fitting a model including 
#' all explanatory variables and examining a spline correlogram of the Pearson 
#' residuals.
#' 
#' glmmML estimates model parameters by maximum likelihood and allows 
#' AIC values to be calculated. However, lmer in lme4 can also be 
#' used.
#' 
#' The cluster argument indicates the grouping level for the random 
#' effect.
#' 

modelML <- glmmML(formula=Triturus_cristatus ~ Lissotriton_vulgaris + Fish + Bufo_bufo
                  + Amphibians + Pond_Area + Gallinula_chloropus + Pond_Substrate
                  + Quality + Pond_Density + Woodland + Permanance + Max_Depth
                  + Outflow + Inflow + Scrub_Hedge + Macrophytes + Shading
                  + Ruderals + Waterfowl + Cyprinus_carpio + 
                  + Gasterosteus_aculeatus + Pungitius_pungitius 
                  + Fulica_atra,
                  cluster = County,
                  family = binomial,
                  data = habitat.dat)

## Check for spatial autocorrelation of the Pearson residuals
Correlog.modelML <- spline.correlog(x=habitat.dat$Easting,
                                    y=habitat.dat$Northing,
                                    z=pres.glmmML(model = modelML,
                                                  data = habitat.dat), 
                                    na.rm=TRUE, 
                                    xmax=10000)
plot(Correlog.modelML)

#' 
#' Spatial autocorrelation has been removed at short lag distances
#' thus both models are successfully accomodating some of the 
#' spatial autocorrelation within sites. However, there is no real
#' difference between a GLM and GLMM therefore mixed effect modeling
#' is not necessarily required to remove spatial autocorrelation. 
#' To examine GCN detection at the pondscape and remove any pontential
#' effects of county though, we have to model county as a random effect.
#' 
#' Now a suitable set of explanatory variables and modelling framework
#' has been identified, we now need to identify which variables are 
#' important determinants of GCN detection and choose a suitable,
#' parsimonious approximating model to make predictions. This is
#' done using an information-theoretic approach using AIC criteria
#' to evaluate model fit. Low AIC models are more parsimonious than
#' high AIC models.
#' 
#' First, group explanatory variables into different functional groups:
#' 
#' - Biotic (smooth newt, toad, amphibian, fish, waterfowl, common moorhen,
#'   common coot, common carp, three-spined stickleback, ninespine stickleback)
#'    
#' - Abiotic (pond substrate, max depth, pond area, water quality, outflow, 
#'   inflow, pond density, permanence, macrophytes, shading, scrub/hedge, 
#'   ruderals, woodland)
#'   
#' We will apply step-wise down nested model selection to each group using 
#' glmer. The drop1() function will be used to select the order that variables 
#' should be dropped from the model. The function computes all the single terms 
#' in the scope argument that can be added to or dropped from the model, fits 
#' those models and computes a table of the changes in fit.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled
GCNdist <- data.frame(habitat.dat[,c(1,5,10,14:15,17:18,20:24,26:29,31:32,
                                     43,49,52:53,55,61,77,87)])


################
# Biotic model #
################

biotic.m0 <- glmer(Triturus_cristatus ~ (1|County) 
                   + Lissotriton_vulgaris + Bufo_bufo + Amphibians 
                   + Fish + Waterfowl + Cyprinus_carpio 
                   + Gasterosteus_aculeatus + Pungitius_pungitius 
                   + Gallinula_chloropus + Fulica_atra,
                   family = binomial(link="logit"),
                   data = GCNdist)

drop1(biotic.m0)

## Remove amphibian presence as this variable has the lowest AIC value
biotic.m1 <- glmer(Triturus_cristatus ~ (1|County)
                   + Lissotriton_vulgaris + Bufo_bufo
                   + Fish + Waterfowl + Cyprinus_carpio 
                   + Gasterosteus_aculeatus + Pungitius_pungitius 
                   + Gallinula_chloropus + Fulica_atra,
                   family = binomial(link="logit"),
                   data = GCNdist)

anova(biotic.m1, biotic.m0)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibian presence has resulted in better model fit

drop1(biotic.m1)

## Remove common coot presence as this variable has the lowest AIC value
biotic.m2 <- glmer(Triturus_cristatus ~ (1|County)
                   + Lissotriton_vulgaris + Bufo_bufo
                   + Fish + Waterfowl + Cyprinus_carpio 
                   + Gasterosteus_aculeatus + Pungitius_pungitius 
                   + Gallinula_chloropus,
                   family = binomial(link="logit"),
                   data = GCNdist)

anova(biotic.m2, biotic.m1)

## p-value from LRT not significant and reduction in AIC value
## so removal of common coot presence has resulted in better model fit

drop1(biotic.m2)

## Remove fish presence as this variable has the lowest AIC value
biotic.m3 <- glmer(Triturus_cristatus ~ (1|County)
                   + Lissotriton_vulgaris + Bufo_bufo
                   + Waterfowl + Cyprinus_carpio 
                   + Gasterosteus_aculeatus + Pungitius_pungitius 
                   + Gallinula_chloropus,
                   family = binomial(link="logit"),
                   data = GCNdist)

anova(biotic.m3, biotic.m2)

## p-value from LRT not significant and reduction in AIC value
## so removal of fish presence has resulted in better model fit

drop1(biotic.m3)

## Remove ninespine stickleback presence as this variable has the lowest 
## AIC value
biotic.m4 <- glmer(Triturus_cristatus ~ (1|County)
                   + Lissotriton_vulgaris + Bufo_bufo + Waterfowl 
                   + Cyprinus_carpio + Gasterosteus_aculeatus
                   + Gallinula_chloropus,
                   family = binomial(link="logit"),
                   data = GCNdist)

anova(biotic.m4, biotic.m3)

## p-value from LRT not significant and reduction in AIC value
## so removal of ninespine stickleback presence has resulted in better model 
## fit

drop1(biotic.m4)

## Remove waterfowl presence as this variable has the lowest AIC value
biotic.m5 <- glmer(Triturus_cristatus ~ (1|County)
                   + Lissotriton_vulgaris + Bufo_bufo + Cyprinus_carpio 
                   + Gasterosteus_aculeatus + Gallinula_chloropus,
                   family = binomial(link="logit"),
                   data = GCNdist)

anova(biotic.m5, biotic.m4)

## p-value from LRT significant and increase in AIC value so removal of 
## waterfowl presence has reduced model fit

drop1(biotic.m4)

## Remove three-spined stickleback presence as this variable has the lowest 
## AIC value after waterfowl presence to see if final model has been reached
biotic.m6 <- glmer(Triturus_cristatus ~ (1|County)
                   + Lissotriton_vulgaris + Bufo_bufo + Waterfowl 
                   + Cyprinus_carpio + Gallinula_chloropus,
                   family = binomial(link="logit"),
                   data = GCNdist)

anova(biotic.m6, biotic.m4)

## p-value from LRT significant and increase in AIC value so removal of 
## three-spined stickleback presence has reduced model fit


#'
#' Variation in GCN occupancy resulting from biotic factors is best explained 
#' by smooth newt presence, common toad presence, common carp presence,
#' three-spined stickleback presence, common moorhen presence, and
#' waterfowl presence.
#' 

# Final biotic model: 
biotic <- glmer(Triturus_cristatus ~ (1|County)
                + Lissotriton_vulgaris + Bufo_bufo + Waterfowl 
                + Cyprinus_carpio + Gasterosteus_aculeatus
                + Gallinula_chloropus,
                family = binomial(link="logit"),
                data = GCNdist)

summary(biotic)
anova(biotic)
drop1(biotic, test = "Chi")   # for binomial/integer y models, 
                              # statistics for significance of each term
display(biotic)
se.ranef(biotic)              # levels of random factor centred around 0

## Calculate R-squared of model
r.squaredGLMM(biotic)
## marginal R-squared = 27.77% (proportion of variance explained by fixed effects)
## conditional R-squared = 27.77% (proportion of variance explained by fixed + random effects)

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs.
## Use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(biotic)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion

## Plot the fitted data against the observed data
plot(GCNdist$Triturus_cristatus ~ fitted(biotic))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(GCNdist$Triturus_cristatus, fitted(biotic))

## Model does not fit well as p-value is significant 
## i.e. significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(biotic, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P < 2.2e-16

## Some deviation from normality as residuals are not normally distributed
## therefore model may not be that reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to identify 
## source of heterogeneity i.e. independent variable that is non-linearly 
## associated with y
plot(sresid ~ GCNdist$Lissotriton_vulgaris)   
plot(sresid ~ GCNdist$Bufo_bufo)
plot(sresid ~ GCNdist$Cyprinus_carpio)
plot(sresid ~ GCNdist$Gasterosteus_aculeatus)
plot(sresid ~ GCNdist$Gallinula_chloropus)
plot(sresid ~ GCNdist$Waterfowl)

## Assumption 3: no collinearity
## All collinearity was removed during preliminary examination of
## data.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by influential observations
## Data point influential if hi<2p/n
leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(GCNdist$Lissotriton_vulgaris)
hi <- (2*2)/10
plot(leverage(GCNdist$Lissotriton_vulgaris), type = "h")
abline(0.4, 0, lty = 2)
points(GCNdist$Lissotriton_vulgaris)

## Make dataframe of p-values to calculate Bonferroni and Benjamini-Hochberg
## corrections to account for Type I error (i.e. false positives)
biotic.results <- data.frame(Variable = c("Lissotriton vulgaris",
                                          "Bufo bufo", "Waterfowl",
                                          "Cyprinus carpio", 
                                          "Gasterosteus aculeatus",
                                          "Gallinula chloropus"),
                             P = c(1.519e-12, 0.0008291, 0.0417882, 0.0050708,
                                   0.0120750, 7.411e-05))

## Calculate Bonferroni and Benjamini-Hochberg corrections
biotic.results$Bonferroni <- round(p.adjust(biotic.results$P, 
                                            method="bonferroni", 
                                            n=length(biotic.results$P)), 3)

biotic.results$Benjamini_Hochberg <- round(p.adjust(biotic.results$P, 
                                                    method="BH", 
                                                    n=length(biotic.results$P)), 3)

## PLOT MODEL
## Obtain predicted values for the full data set: GCNdist
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(biotic, newdata=GCNdist, 
                            se.fit=TRUE, print.matrix=T))

# Create new data set with fitted values and original data for plots 1 and 2
box.dat <- cbind(GCNdist, fit) 

## Plot smooth newt presence-absence and GCN presence-absence 
p12 <- ggplot(box.dat) + ggtitle("(a)")
p12 <- p12 + geom_jitter(aes(x=factor(Lissotriton_vulgaris), 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         width=0.3, height=0.3, cex=1)
p12 <- p12 + geom_boxplot(aes(x=factor(Lissotriton_vulgaris), 
                              y=fit), 
                          alpha=0.7, outlier.shape=NA)
p12 <- p12 + scale_colour_gradient(name=expression(paste(italic("T. cristatus"), " occupancy")),
                                   low="grey45", 
                                   high="orange")
p12 <- p12 + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Absent", "Present"))
p12 <- p12 + labs(x = expression(paste(italic("L. vulgaris"))),
                  y = expression(paste("Probability of ", italic("T. cristatus"))))
p12 <- p12 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p12

## Plot common toad presence-absence and GCN presence-absence 
p13 <- ggplot(box.dat) + ggtitle("(b)")
p13 <- p13 + geom_jitter(aes(x=factor(Bufo_bufo), 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         width=0.3, height=0.3, cex=1)
p13 <- p13 + geom_boxplot(aes(x=factor(Bufo_bufo), 
                              y=fit), 
                          alpha=0.7, outlier.shape=NA)
p13 <- p13 + scale_colour_gradient(name=expression(paste(italic("T. cristatus"), " occupancy")),
                                   low="grey45", 
                                   high="orange")
p13 <- p13 + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Absent", "Present"))
p13 <- p13 + labs(x = expression(paste(italic("B. bufo"))),
                  y = "")
p13 <- p13 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p13

## Plot common carp presence-absence and GCN presence-absence 
p14 <- ggplot(box.dat) + ggtitle("(c)")
p14 <- p14 + geom_jitter(aes(x=factor(Cyprinus_carpio), 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         width=0.3, height=0.3, cex=1)
p14 <- p14 + geom_boxplot(aes(x=factor(Cyprinus_carpio), 
                              y=fit), 
                          alpha=0.7, outlier.shape=NA)
p14 <- p14 + scale_colour_gradient(name=expression(paste(italic("T. cristatus"), " occupancy")),
                                   low="grey45", 
                                   high="orange")
p14 <- p14 + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Absent", "Present"))
p14 <- p14 + labs(x = expression(paste(italic("C. carpio"))),
                  y = expression(paste("Probability of ", italic("T. cristatus"))))
p14 <- p14 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p14

## Plot three-spined stickleback presence-absence and GCN presence-absence 
p15 <- ggplot(box.dat) + ggtitle("(d)")
p15 <- p15 + geom_jitter(aes(x=factor(Gasterosteus_aculeatus), 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         width=0.3, height=0.3, cex=1)
p15 <- p15 + geom_boxplot(aes(x=factor(Gasterosteus_aculeatus), 
                              y=fit), 
                          alpha=0.7, outlier.shape=NA)
p15 <- p15 + scale_colour_gradient(name=expression(paste(italic("T. cristatus"), " occupancy")),
                                   low="grey45", 
                                   high="orange")
p15 <- p15 + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Absent", "Present"))
p15 <- p15 + labs(x = expression(paste(italic("G. aculeatus"))),
                  y = "")
p15 <- p15 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p15

## Plot common moorhen presence-absence and GCN presence-absence 
p16 <- ggplot(box.dat) + ggtitle("(e)")
p16 <- p16 + geom_jitter(aes(x=factor(Gallinula_chloropus), 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         width=0.3, height=0.3, cex=1)
p16 <- p16 + geom_boxplot(aes(x=factor(Gallinula_chloropus), 
                              y=fit), 
                          alpha=0.7, outlier.shape=NA)
p16 <- p16 + scale_colour_gradient(name=expression(paste(italic("T. cristatus"), " occupancy")),
                                   low="grey45", 
                                   high="orange")
p16 <- p16 + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Absent", "Present"))
p16 <- p16 + labs(x = expression(paste(italic("G. chloropus"))),
                  y = expression(paste("Probability of ", italic("T. cristatus"))))
p16 <- p16 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p16

## Plot waterfowl presence-absence and GCN presence-absence 
f.waterfowl <- ordered(GCNdist$Waterfowl, 
                       levels = c("Absent", "Minor", "Major"))
f.waterfowl

p17 <- ggplot(box.dat) + ggtitle("(f)")
p17 <- p17 + geom_jitter(aes(x=f.waterfowl, 
                             y=Triturus_cristatus, 
                             colour=Triturus_cristatus), 
                         width=0.3, height=0.3, cex=1)
p17 <- p17 + geom_boxplot(aes(x=f.waterfowl, 
                              y=fit), 
                          alpha=0.7, outlier.shape=NA)
p17 <- p17 + scale_colour_gradient(name=expression(paste(italic("T. cristatus"), " occupancy")),
                                   low="grey45", 
                                   high="orange")
p17 <- p17 + labs(x = "Waterfowl",
                  y = "")
p17 <- p17 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p17


#################
# Abiotic model #
#################

abiotic.m0 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Max_Depth + Pond_Density + Inflow
                    + Pond_Area + Pond_Substrate + Outflow
                    + Macrophytes + Quality + Permanance
                    + Shading + Ruderals + Scrub_Hedge + Woodland,
                    family = binomial(link="logit"),
                    data = GCNdist)

drop1(abiotic.m0)

## Remove pond substrate from model as this has the lowest AIC value
abiotic.m1 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Max_Depth + Pond_Density + Inflow + Pond_Area  
                    + Outflow + Macrophytes + Quality + Permanance
                    + Shading + Ruderals + Scrub_Hedge + Woodland,
                    family = binomial(link="logit"),
                    data = GCNdist)

anova(abiotic.m1, abiotic.m0)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond susbtrate has resulted in better model fit.

drop1(abiotic.m1)

## Remove scrub/hedge from model as this has the lowest AIC value
abiotic.m2 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Max_Depth + Pond_Density + Inflow + Pond_Area  
                    + Outflow + Macrophytes + Quality + Permanance
                    + Shading + Ruderals + Woodland,
                    family = binomial(link="logit"),
                    data = GCNdist)

anova(abiotic.m2, abiotic.m1)

## p-value from LRT not significant and reduction in AIC value
## so removal of scrub/hedge has resulted in better model fit.

drop1(abiotic.m2)

## Remove permanance from model as this has the lowest AIC value
abiotic.m3 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Max_Depth + Pond_Density + Inflow + Pond_Area  
                    + Outflow + Macrophytes + Quality
                    + Shading + Ruderals + Woodland,
                    family = binomial(link="logit"),
                    data = GCNdist)

anova(abiotic.m3, abiotic.m2)

## p-value from LRT not significant and reduction in AIC value
## so removal of permanance has resulted in better model fit.

drop1(abiotic.m3)

## Remove woodland from model as this has the lowest AIC value
abiotic.m4 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Max_Depth + Pond_Density + Inflow + Pond_Area  
                    + Outflow + Macrophytes + Quality
                    + Shading + Ruderals,
                    family = binomial(link="logit"),
                    data = GCNdist)

anova(abiotic.m4, abiotic.m3)

## p-value from LRT not significant and reduction in AIC value
## so removal of woodland has resulted in better model fit.

drop1(abiotic.m4)

## Remove macrophyte cover from model as this has the lowest AIC value
abiotic.m5 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Max_Depth + Pond_Density + Inflow + Pond_Area  
                    + Outflow + Quality + Shading + Ruderals,
                    family = binomial(link="logit"),
                    data = GCNdist)

anova(abiotic.m5, abiotic.m4)

## p-value from LRT not significant and reduction in AIC value
## so removal of macrophyte cover has resulted in better model fit.

drop1(abiotic.m5)

## Remove outflow from model as this has the lowest AIC value
abiotic.m6 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Max_Depth + Pond_Density + Inflow + Pond_Area  
                    + Quality + Shading + Ruderals,
                    family = binomial(link="logit"),
                    data = GCNdist)

anova(abiotic.m6, abiotic.m5)

## p-value from LRT not significant and reduction in AIC value
## so removal of outflow has resulted in better model fit.

drop1(abiotic.m6)

## Remove pond density from model as this has the lowest AIC value
abiotic.m7 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Max_Depth + Inflow + Pond_Area + Quality 
                    + Shading + Ruderals,
                    family = binomial(link="logit"),
                    data = GCNdist)

anova(abiotic.m7, abiotic.m6)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond density has resulted in better model fit.

drop1(abiotic.m7)

## Remove ruderals from model as this has the lowest AIC value
abiotic.m8 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Max_Depth + Inflow + Pond_Area + Quality + Shading,
                    family = binomial(link="logit"),
                    data = GCNdist)

anova(abiotic.m8, abiotic.m7)

## p-value from LRT not significant and reduction in AIC value
## so removal of ruderals has resulted in better model fit.

drop1(abiotic.m8)

## Remove maximum depth from model as this has the lowest AIC value
abiotic.m9 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Inflow + Pond_Area + Quality + Shading,
                    family = binomial(link="logit"),
                    data = GCNdist)

anova(abiotic.m9, abiotic.m8)

## p-value from LRT not significant and increase in AIC value (<2)
## so removal of maximum depth has resulted in better model fit.

drop1(abiotic.m9)

## Remove water quality from model as this has the lowest AIC value
abiotic.m10 <- glmer(Triturus_cristatus ~ (1|County) 
                    + Inflow + Pond_Area + Shading,
                    family = binomial(link="logit"),
                    data = GCNdist)

anova(abiotic.m10, abiotic.m9)

## p-value from LRT not significant and increase in AIC value (<2)
## so removal of water quality has resulted in better model fit.

drop1(abiotic.m10)

## Remove pond area from model as this has the lowest AIC value
abiotic.m11 <- glmer(Triturus_cristatus ~ (1|County) 
                     + Inflow + Shading,
                     family = binomial(link="logit"),
                     data = GCNdist)

anova(abiotic.m11, abiotic.m10)

## p-value from LRT significant and increase in AIC value (>2)
## so removal of pond area has reduced model fit.

drop1(abiotic.m10)

## Remove inflow from model as this has the next lowest AIC value
## after pond area to see if final model has been reached
abiotic.m12 <- glmer(Triturus_cristatus ~ (1|County) 
                     + Inflow + Pond_Area,
                     family = binomial(link="logit"),
                     data = GCNdist)

anova(abiotic.m12, abiotic.m10)

## p-value from LRT significant and increase in AIC value (>2)
## so removal of inflow has reduced model fit.


## Final abiotic model:
abiotic <- glmer(Triturus_cristatus ~ (1|County)
                 + Inflow + Pond_Area + Shading,
                 family = binomial(link="logit"),
                 data = GCNdist)

summary(abiotic)
anova(abiotic)
drop1(abiotic, test = "Chi")   # for binomial/integer y models, 
                               # statistics for significance of each term
display(abiotic)
se.ranef(abiotic)              # levels of random factor centred around 0

## Calculate R-squared of model
r.squaredGLMM(abiotic)
## marginal R-squared = 9.27% (proportion of variance explained by fixed effects)
## conditional R-squared = 10.00% (proportion of variance explained by fixed + random effects)

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs.
## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(abiotic)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion

## Plot the fitted data against the observed data
plot(GCNdist$Triturus_cristatus ~ fitted(abiotic))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(GCNdist$Triturus_cristatus, fitted(abiotic))

## Model does not fit well as p-value is significant 
## i.e. significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(abiotic, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P < 2.2e-16

## Some deviation from normality as residuals are not normally distributed
## therefore model may not be that reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to identify 
## source of heterogeneity i.e. independent variable that is non-linearly 
## associated with y
plot(sresid ~ GCNdist$Inflow)   
plot(sresid ~ GCNdist$Pond_Area)
plot(sresid ~ GCNdist$Shading)

## Assumption 3: no collinearity
## All collinearity was removed during preliminary examination of
## data.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by influential observations
## Data point influential if hi<2p/n
leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(GCNdist$Lissotriton_vulgaris)
hi <- (2*2)/10
plot(leverage(GCNdist$Lissotriton_vulgaris), type = "h")
abline(0.4, 0, lty = 2)
points(GCNdist$Lissotriton_vulgaris)

## Make dataframe of p-values to calculate Bonferroni and Benjamini-Hochberg
## corrections to account for Type I error (i.e. false positives)
abiotic.results <- data.frame(Variable = c("Inflow", "Pond Area", "Shading"),
                              P = c(0.0007199, 0.0420322, 0.0004612))

## Calculate Bonferroni and Benjamini-Hochberg corrections
abiotic.results$Bonferroni <- round(p.adjust(abiotic.results$P, 
                                             method="bonferroni", 
                                             n=length(abiotic.results$P)), 3)

abiotic.results$Benjamini_Hochberg <- round(p.adjust(abiotic.results$P, 
                                                     method="BH", 
                                                     n=length(abiotic.results$P)), 3)


## PLOT MODEL
## Obtain predicted values for the full data set: GCNdist
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
abiotic.fit <- data.frame(predictSE(abiotic, newdata=GCNdist, 
                                    se.fit=TRUE, print.matrix=T))

# Create new data set with fitted values and original data for plots 1 and 2
abiotic.dat <- cbind(GCNdist, abiotic.fit) 

## Plot inflow and GCN presence-absence
p18 <- ggplot(abiotic.dat) + ggtitle("(g)")
p18 <- p18 + geom_jitter(aes(x=Inflow, 
                             y=Triturus_cristatus, 
                             colour=GCNdist$Triturus_cristatus), 
                         width=0.3, height=0.3, cex=1)
p18 <- p18 + geom_boxplot(aes(x=Inflow, 
                              y=fit), 
                          alpha=0.7, outlier.shape=NA)
p18 <- p18 + scale_colour_gradient(low="grey45", high="orange")
p18 <- p18 + labs(x = "Inflow", 
                  y = expression(paste("Probability of ", italic("T. cristatus"))))  
p18 <- p18 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p18

## Plot pond area and GCN presence-absence
## Create a range of pond area values which increase by 18.6012
## to predict GCN occupancy as pond area increases
range <- seq(from=min(GCNdist$Pond_Area), 
             to=max(GCNdist$Pond_Area), 
             by=18.6012)

## Create new data frame where only pond area changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset
d1 <- data.frame(Sample=rep("F1", length(range)),
                 Pond_Area=range,
                 Inflow=rep("Present", length(range)),
                 Shading=rep(40.95238, length(range)))

## Get predictions for this new dataset
fit <- data.frame(predictSE(abiotic, newdata=d1, 
                            se.fit=TRUE, print.matrix=T)) 

## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

## New dataframe with upper and lower 95% CIs
dat.PA <- cbind(d1, fit, ciu, cil)

## Plot showing relationship between pond area and GCN presence-absence 
## with predicted values from model
p19 <- ggplot() + ggtitle("(h)")
p19 <- p19 + geom_jitter(aes(x=GCNdist$Pond_Area,
                             y=GCNdist$Triturus_cristatus, 
                             colour=GCNdist$Triturus_cristatus), 
                         height=0.3, cex=1)
p19 <- p19 + geom_line(aes(x=dat.PA$Pond_Area, 
                           y=dat.PA$fit), 
                       size = 1)
p19 <- p19 + geom_ribbon(aes(x=dat.PA$Pond_Area, 
                             ymin = dat.PA$cil, 
                             ymax = dat.PA$ciu), 
                         alpha = 0.25)
p19 <- p19 + scale_colour_gradient(low="grey45", high="orange")
p19 <- p19 + coord_cartesian(xlim=c(0,10000))
p19 <- p19 + labs(x = "Pond area ("~m^2~")", 
                  y = "")
p19 <- p19 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p19


## Plot percentage of shading and GCN presence-absence
## Create a range of shading values which increase by 0.1984127
## to predict GCN occupancy as shading increases
range <- seq(from=min(GCNdist$Shading), 
             to=max(GCNdist$Shading), 
             by=0.1984127)

## Create new data frame where only percentage of shading changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset
d2 <- data.frame(Sample=rep("F1", length(range)),
                 Pond_Area=rep(599.9762, length(range)),
                 Inflow=rep("Present", length(range)),
                 Shading=range)

## Get predictions for this new dataset
fit <- data.frame(predictSE(abiotic, newdata=d2, 
                            se.fit=TRUE, print.matrix=T)) 

## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

## New dataframe with upper and lower 95% CIs
dat.S <- cbind(d2, fit, ciu, cil)

## Plot showing relationship between pond area and GCN presence-absence 
## with predicted values from model
p20 <- ggplot() + ggtitle("(i)")
p20 <- p20 + geom_jitter(aes(x=GCNdist$Shading,
                             y=GCNdist$Triturus_cristatus, 
                             colour=GCNdist$Triturus_cristatus), 
                         height=0.3, cex=1)
p20 <- p20 + geom_line(aes(x=dat.S$Shading, 
                           y=dat.S$fit), 
                       size = 1)
p20 <- p20 + geom_ribbon(aes(x=dat.S$Shading, 
                             ymin = dat.S$cil, 
                             ymax = dat.S$ciu), 
                         alpha = 0.25)
p20 <- p20 + scale_colour_gradient(low="grey45", high="orange")
p20 <- p20 + coord_cartesian(xlim=c(0,100))
p20 <- p20 + labs(x = "Shading (%)", 
                  y = expression(paste("Probability of ", italic("T. cristatus"))))
p20 <- p20 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p20



#' ---
#' 
#' Now, test HSI score for correlation with GCN detection probability.
#' 

## GLMM:
HSIm <- glmer(Triturus_cristatus ~ (1|County) + HSI,
              family = binomial(link="logit"),
              data = habitat.dat)

## Compare to null model:
HSInull <- glmer(Triturus_cristatus ~ (1|County) + 1,
                 family = binomial(link="logit"),
                 data = habitat.dat)

anova(HSIm, HSInull)

## Model with HSI score is better fit to the data than the null model.

## Model output:
summary(HSIm)
anova(HSIm)
drop1(HSIm, test="Chi")  # for binomial/integer y models, 
                         # statistics for significance of each term
display(HSIm)
se.ranef(HSIm)           # levels of random factor centred around 0

## Calculate R-squared of model
r.squaredGLMM(HSIm)
## marginal R-squared = 5.14% (proportion of variance explained by fixed effects)
## conditional R-squared = 6.66% (proportion of variance explained by fixed + random effects)

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs.
## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSIm)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion

## Plot the fitted data against the observed data
plot(GCNdist$Triturus_cristatus ~ fitted(HSIm))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(GCNdist$Triturus_cristatus, fitted(HSIm))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(HSIm, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P < 2.2e-16

## Some deviation from normality as residuals are not normally distributed
## therefore model may not be that reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to identify 
## source of heterogeneity i.e. independent variable that is non-linearly 
## associated with y
plot(sresid ~ habitat.dat$HSI)   

## Assumption 3: no collinearity
## Only one variable being modelled thus no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by influential observations
## Data point influential if hi<2p/n
leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(habitat.dat$HSI)
hi <- (2*2)/10
plot(leverage(habitat.dat$HSI), type = "h")
abline(0.4, 0, lty = 2)
points(habitat.dat$HSI)


## PLOT MODEL
## Obtain predicted values for the full data set: habitat.dat
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
HSI.fit <- data.frame(predictSE(HSIm, newdata=habitat.dat, 
                                se.fit=TRUE, print.matrix=T))

# Create new data set with fitted values and original data for plots 1 and 2
HSI.dat <- cbind(habitat.dat, HSI.fit) 


## Plot HSI score and GCN presence-absence
## Create a range of HSI score values which increase by 0.00142
## to predict GCN occupancy as HSI score increases
range <- seq(from=min(habitat.dat$HSI), 
             to=max(habitat.dat$HSI), 
             by=0.00142)

## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset
d3 <- data.frame(sample=rep("F1", length(range)),
                 HSI=range)

## Get predictions for this new dataset
fit <- data.frame(predictSE(HSIm, newdata=d3, 
                            se.fit = TRUE, print.matrix=T)) 

## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

## New dataframe containing upper and lower 95% CIs
dat.HSI <- cbind(d3, fit, ciu, cil) # This is now the data frame

# Plot showing relationship between HSI score and GCN presence-absence 
# with predicted values from model
p21 <- ggplot() + ggtitle("(j)")
p21 <- p21 + geom_jitter(aes(x=habitat.dat$HSI, 
                             y=habitat.dat$Triturus_cristatus, 
                             colour=habitat.dat$Triturus_cristatus), 
                         height=0.3, cex=1)
p21 <- p21 + geom_line(aes(x=dat.HSI$HSI, 
                           y=dat.HSI$fit), 
                       size = 1)
p21 <- p21 + geom_ribbon(aes(x=dat.HSI$HSI, 
                             ymin = dat.HSI$cil, 
                             ymax = dat.HSI$ciu), 
                         alpha = 0.25)
p21 <- p21 + scale_colour_gradient(low="grey45", high="orange")
p21 <- p21 + coord_cartesian(xlim=c(0,1))
p21 <- p21 + labs(x = "Habitat Suitability Index (HSI) score",
                  y = "")
p21 <- p21 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p21

## Combine HSI score plot with plots from biotic and abiotic GLMMs
## Load 'ggpubr' package
library(ggpubr)
p12_21 <- ggarrange(arrangeGrob(p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,
                                ncol=2, nrow=5),
                    mylegend, nrow=2, heights=c(15,1))
p12_21



#' ---
#' 
#' Does the GCN deserve its status as a umbrella species for pond biodiversity?
#' 
#' If the GCN is an indictor species for vertebrate biodiversity, then 
#' vertebrate species richness should be positively associated with GCN 
#' detection and the GCN habitat suitability index score (i.e. habitat 
#' suitable for GCN is suitable for other species).
#' 
#' Model vertebrate species richness against GCN detection and HSI score 
#' in a GLMM with Poisson distribution.
#' 

SRm <- glmer(Sp.richness ~ (1|County) + Triturus_cristatus + HSI,
             family=poisson,
             data=habitat.dat)

## Compare to null model
SRnull <- glmer(Sp.richness ~ (1|County) + 1,
                family=poisson,
                data=habitat.dat)

anova(SRm, SRnull)

## LRT p-value significant and AIC value of null model greater.
## Therefore, GCN presence-absence and HSI score are appropriate 
## explanatory variables.

## Model output:
summary(SRm)
anova(SRm)
drop1(SRm, test="Chi")  # for binomial/integer y models, 
                        # statistics for significance of each term
display(SRm)
se.ranef(SRm)           # levels of random factor centred around 0

## Calculate R-squared of model
r.squaredGLMM(SRm)
## marginal R-squared = 8.70% (proportion of variance explained by fixed effects)
## conditional R-squared = 9.14% (proportion of variance explained by fixed + random effects)

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs.
## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(SRm)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion

## Plot the fitted data against the observed data
plot(habitat.dat$Sp.richness ~ fitted(SRm))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(habitat.dat$Sp.richness, fitted(SRm))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(SRm, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 1.407e-12

## Some deviation from normality as residuals are not normally distributed
## therefore model may not be that reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to identify 
## source of heterogeneity i.e. independent variable that is non-linearly 
## associated with y
plot(sresid ~ habitat.dat$HSI) 
plot(sresid ~ habitat.dat$Triturus_cristatus)

## Assumption 3: no collinearity
vif(SRm)

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by influential observations
## Data point influential if hi<2p/n
leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(habitat.dat$HSI)
hi <- (2*2)/10
plot(leverage(habitat.dat$HSI), type = "h")
abline(0.4, 0, lty = 2)
points(habitat.dat$HSI)

leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(habitat.dat$Triturus_cristatus)
hi <- (2*2)/10
plot(leverage(habitat.dat$Triturus_cristatus), type = "h")
abline(0.4, 0, lty = 2)
points(habitat.dat$Triturus_cristatus)

## Make dataframe of p-values to calculate Bonferroni and Benjamini-Hochberg
## corrections to account for Type I error (i.e. false positives)
SR.results <- data.frame(Variable = c("Triturus cristatus", "HSI"),
                         P = c(1.033e-11, 0.268))

## Calculate Bonferroni and Benjamini-Hochberg corrections
SR.results$Bonferroni <- round(p.adjust(SR.results$P, 
                                        method="bonferroni", 
                                        n=length(SR.results$P)), 3)

SR.results$Benjamini_Hochberg <- round(p.adjust(SR.results$P, 
                                                method="BH", 
                                                n=length(SR.results$P)), 3)


## PLOT MODEL
## Obtain predicted values for the full data set: GCNdist
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
SR.fit <- data.frame(predictSE(SRm, newdata=habitat.dat,
                               se.fit=TRUE, print.matrix=T))

## Create new data set with fitted values and original data for plots 1 and 2
SR.dat <- cbind(habitat.dat, SR.fit) 

## Plot relationship between vertebrate species richness and GCN 
## presence-absence
p22 <- ggplot(SR.dat) + ggtitle("(a)")
p22 <- p22 + geom_jitter(aes(x=factor(Triturus_cristatus), 
                             y=Sp.richness), 
                         colour="black", 
                         width=0.3, 
                         height=0.3, 
                         cex=1)
p22 <- p22 + geom_boxplot(aes(x=factor(Triturus_cristatus), 
                              y=fit), 
                          alpha=0.7, 
                          outlier.shape=NA)
p22 <- p22 + scale_y_continuous(limits=c(0,10), breaks=seq(0,10,2))
p22 <- p22 + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Negative", "Positive"))
p22 <- p22 + labs(x = expression(italic("T. cristatus")), 
                  y = "Vertebrate species richness") 
p22 <- p22 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p22


## Plot HSI score and vertebrate species richness
## Create a range of HSI score values which increase by 0.00142
## to predict GCN occupancy as HSI score increases
range <- seq(from=min(habitat.dat$HSI), 
             to=max(habitat.dat$HSI), 
             by=0.00142)

## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset
d4 <- data.frame(Sample=rep("F1", length(range)),
                 Triturus_cristatus=rep(1, length(range)),
                 HSI=range)

## Get predictions for this new dataset
fit <- data.frame(predictSE(SRm, newdata=d4, 
                            se.fit = TRUE, print.matrix=T)) 

## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

## New dataframe containing upper and lower 95% CIs
dat.SR <- cbind(d4, fit, ciu, cil) # This is now the data frame

## Plot showing relationship between HSI score and species richness
## with predicted values from model
p23 <- ggplot() + ggtitle("(b)")
p23 <- p23 + geom_jitter(aes(x=habitat.dat$HSI, 
                             y=habitat.dat$Sp.richness), 
                         colour="black", 
                         height=0.7, 
                         cex=1)
p23 <- p23 + geom_line(aes(x=dat.SR$HSI, 
                           y=dat.SR$fit), 
                       size = 1)
p23 <- p23 + geom_ribbon(aes(x=dat.SR$HSI, 
                             ymin = dat.SR$cil, 
                             ymax = dat.SR$ciu), 
                         alpha = 0.25)
p23 <- p23 + scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.1))
p23 <- p23 + scale_y_continuous(limits=c(0,10), breaks=seq(0,10,2))
p23 <- p23 + labs(y = "", 
                  x = "Habitat Suitability Index (HSI) score")
p23 <- p23 + theme(panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p23

## SHOW ALL PLOTS
ggarrange(p22, p23, ncol=2, nrow=1)
                      
