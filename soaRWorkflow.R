rm(list=ls())

options(java.parameters = c("-Xmx30000m"))
gc()
#library(devtools)

devtools::install_github("https://github.com/robboyd/soaR")

library(BRCmap)
library(randomForest)
library(ggplot2)
library(colorRamps)
library(raster)
library(dismo)
library(rgdal)
library(maptools)
library(sp)
library(sf)
library(maxnet)
library(zoon)
library(PresenceAbsence)
library(soaR)

## source functions to save outputs for validation

source("G:/SDMValidation/R/extractData.R")

taxa <- "Ants"

file <- list.files("G:/TSDA/SDMs/Data/BRCRecords/",
                   full.names = T,
                   pattern = taxa)

load(file)


## format BRC records for use with soaR

dat <- formatData(occ = taxa_data,
                  taxa = taxa,
                  minYear = 2000,
                  maxYear = 2015,
                  degrade = TRUE)

#save(dat, file = paste0("G:/TSDA/SDMs/Data/occurrence/", taxa, ".rdata"))

load(paste0("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/TSDA/SDMs/Data/occurrence/2022/", taxa, ".rdata"))

## check for duplicates

any(duplicated(dat))

if(any(duplicated(dat))) dat <- dat[-which(duplicated(dat)), ]

## get species names in group

spp <- getSpNames(inPath = "W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/TSDA/SDMs/Data/occurrence/2022/",
                  taxa = taxa)

## load climate and land cover covariates

lcVars <- raster::stack(list.files("D:/TSDA/SDMs/Data/targetLandCover2007/selectedCovs/", full.names = T))

climVars <- raster::stack(list.files("D:/TSDA/SDMs/Data/UKBioClimPCA/",
                             full.names = T,
                             pattern = ".asc"))

elev <- raster::raster("D:/TSDA/SDMs/Data/elevation/UKelv.asc")

climVars <- raster::crop(climVars, raster::extent(0, 7e+5, 0, 1250000))

lcVars <- raster::crop(lcVars, raster::extent(0, 7e+5, 0, 1250000))

elev <- raster::crop(elev, raster::extent(0, 7e+5, 0, 1250000))

covs <- raster::stack(lcVars, climVars, elev)

screenRast <- raster::stack(elev, climVars[[1]], lcVars[[1]])


## create pseudo absences and combine with presence data

for (i in c("lrReg", "else")) {

        if (i == "lrReg") {

                matchPres <- FALSE

                nAbs <- 10000

        } else {

                matchPres <- TRUE

                nAbs <- NULL

        }

        spDat <- lapply(X = spp,
                        FUN = createPresAb,
                        dat = dat,
                        taxon = taxa,
                        matchPres = matchPres,
                        nAbs = nAbs,
                        minYear = 2000,
                        maxYear = 2015,
                        recThresh = 10,
                        screenRast = screenRast)

        names(spDat) <- spp

        outPath <- ifelse(i == "lrReg", paste0("D:/TSDA/SDMs/Data/presAb/", taxa, "_", 10000, ".rdata"), paste0("D:/TSDA/SDMs/Data/presAb/", taxa, ".rdata"))

        save(spDat, file = outPath)

        rm(spDat) # clean up

}

## fit SDMs

model <- "lrReg" # "lrReg" or "else"

## load 10,000 pseudo absences if model == lrReg and equal number to presences otherwise

file <- ifelse(model == "lrReg", paste0("D:/TSDA/SDMs/Data/presAb/", taxa, "_", 10000, ".rdata"), paste0("D:/TSDA/SDMs/Data/presAb/", taxa, ".rdata"))

load(file)

nrecs <- lapply(X =spp,
                FUN = function(x) {
                        ind <- which(names(spDat) == x)
                        data.frame(species = x,
                                nrecs = ifelse(is.null(nrow(spDat[[ind]]$Presence)), 0, nrow(spDat[[ind]]$Presence)))
                })

nrecs <- do.call("rbind", nrecs)

spp <- nrecs$species[nrecs$nrecs >= 10]

#alreadyGot <- list.files("G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/Moths/", pattern = "_rf.rdata")
#alreadyGot <- gsub("_rf.rdata", "", alreadyGot)
#alreadyGot

#spp <- spp[!spp %in% alreadyGot]

lapply(X = spp,
       FUN = fitSDM,
       model = model,
       envDat = covs,
       spDat = spDat,
       k = 5,
       write = TRUE,
       plot = TRUE,
       predict = TRUE,
       outPath = paste0("D:/TSDA/SDMs/Data/SDMOutputs_April_2022/", taxa, "/"))

## get model AUC scores

skill <- purrr::map_df(.x = taxa,
                       .f = soaR::getSkill,
                       inPath = "G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/",
                       write = TRUE,
                       outPath = "G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/Skill/",
                       models = c("lrReg", "rf", "max"))


write.csv(skill, paste0("G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/Skill/", taxa, ".csv"),
          row.names = F)

#skill <- read.csv(paste0("G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/Skill/", taxa, ".csv"))

## create AUC-weighted ensembles


lapply(X = unique(skill$species)[453:length(unique(skill$species))],
       FUN = modelAverage,
       inPath = "G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/",
       outPath = "G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/",
       skillDat =skill,
       plot = TRUE,
       models = c("lrReg", "rf", "max"),
       skillThresh = 0.5)

createMeans <- function(species) {

        if (file.exists(paste0("G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/Moths/", species, "_rf.rdata"))) {

            load(paste0("G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/Moths/", species, "_rf.rdata"))

                if(length(out$meanPrediction) < 870000) {

                        print(out$meanPrediction)

                        print(paste("no ensemble for", species))

                        out$meanPrediction <- mean(out$predictions)

                        save(out, file = paste0("G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/Moths/", species, "_rf.rdata"))

                }
        }



}

lapply(unique(skill$species),
       createMeans)

## generate binary presence absence maps

files <- list.files(paste0("G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/", taxa, "/"),
                    pattern = "ensemble.asc")

spNames <- gsub("_ensemble.asc",
                "",
                files)

optStats <- purrr::map_df(.x = spNames,
                          .f = optOccThresh,
                          group = taxa,
                          presAb = spDat,
                          inPath = "G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/",
                          outPath = paste0("G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/", taxa, "/"),
                          map = TRUE)

save(optStats, file = paste0("G:/TSDA/SDMs/Data/skillStatistics/", taxa, ".rdata"))


## save continuous probabilities as a raster stack

files <- list.files(paste0("G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/", taxa, "/"),
                    pattern = "_ensemble",
                    full.names = T)

spp <- list.files(paste0("G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/", taxa, "/"),
                  pattern = "_ensemble")

spp <- gsub("_ensemble.asc", "", spp)

stack <- stack(files)

names(stack) <- spp

rich <- sum(stack)

plot(rich)

writeRaster(rich, filename = "G:/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/Bees/richness.asc",
            format = "ascii")

stack <- readAll(stack)

save(stack, file = paste0("G:/SDMValidation/", taxa, "/continuousPreds.rdata"))

rm(stack)

## plot records

library(BRCmap)

data(UK)

gb <- UK$britain

ir <- UK$ireland

col <- adjustcolor( "blue", alpha.f = 0.5)

load(paste0("G:/TSDA/SDMs/Data/occurrence/", taxa, ".rdata"))

lapply(spp,
       plotRecs,
       outPath = "G:/SDMValidation/")


## plot binary maps
#
#lapply(spp,
#       plotPresAb,
#       outPath = "G:/SDMValidation/")
#
## extract and save prevalence

file <- list.files("G:/TSDA/SDMs/Data/BRCRecords/",
                   full.names = T,
                   pattern = taxa)

load(file)

getMonad <- function(taxa, species) {

        data.frame(group = taxa,
                   charlieName = species,
                   monads = nrow(dat[dat$species == species, ]))

}

monads <- purrr::map_df(.x = spp,
                        .f = getMonad,
                        taxa = taxa)


range <- purrr::map_df(.x = spp,
                      .f = getRecRange,
                      taxa = taxa)

mean(c(2,3,4,5,6,7,8))

range <- lapply(spp,
                getRecRange,
                taxa = taxa)

range <- do.call("rbind", range)

lookup <- read.csv("G:/masterLookup.csv",
                   stringsAsFactors = F)


getLatName <- function(species, lookup) {

        if (species %in% lookup$CONCEPT) {

                lat <- lookup$NAME[which(lookup$CONCEPT == species)]

        } else {

                lat <- species

        }

        out <- data.frame(charlieName = species,
                          latName = lat)


}

spLookup <- lapply(spp,
                   getLatName,
                   lookup = lookup)

spLookup <- do.call("rbind", spLookup)

sppInfo <- merge(monads, range, by = c("charlieName", "group"))
sppInfo <- merge(sppInfo, spLookup, by = "charlieName")
sppInfo

write.csv(sppInfo, file = paste0("G:/SDMValidation/", taxa, "/sppInfo.csv"),
          row.names = F)

load(paste0("G:/SDMValidation/",taxa, "/continuousPreds.rdata"))

par(mfrow = c(1,2))

visSpp <- names(stack)

visSpp <- gsub(" ", ".", spp)

visSpp <- gsub("-", ".", visSpp)

for (i in 1:length(spp)) {

        plot(stack[[visSpp[i]]])

        points(dat$eastings[dat$species == spp[i]], dat$northings[dat$species == spp[i]], cex = 0.2)

        #points(spDat[[spp[i]]]$pseudoAbsence$eastings, spDat[[spp[i]]]$pseudoAbsence$northings, col = "red", cex = 0.2)

        Sys.sleep(7)

}

rm(stack)
getMonad
library(soaR)
