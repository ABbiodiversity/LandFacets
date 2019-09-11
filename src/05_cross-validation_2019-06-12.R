#
# Title: Leave one out cross validation of species models
# Created: June 12th, 2019
# Last Updated: September 4th, 2019
# Author: Brandon Allen
# Objectives: Perform a leave one out cross validation to determine the stability of the two types of models
# Keywords: Facet models, Facet validate, Soil models, Soil validate
#

################
# Facet models # 
################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load relevant libraries, scripts, and lookups
library(mgcv)  # For binomial GAM
library(mapproj)  # For projected maps
library(binom)  # For exact binomial confidence intervals
library(MuMIn)
library(RcmdrMisc)
library(pROC)
library(arm)

source("src/cross-validation-functions_2019-06-17.R")

# Load all relevant datasets
# Prediction grid
pm <- read.csv("data/lookup/prediction-matrix/predictionmatrix_terrain-water-dunes.csv")
rownames(pm) <- colnames(pm)

# Landscape, climate, space
load("data/processed/facet-soil-proportions_2019-08-16.Rdata")
load("data/processed/spatial-climate_2019-08-16.Rdata")

# Add the number of visits to the climate data
site.visit <- table(droplevels(climate.raw$abmi_site))
colnames(climate.raw)[4] <- "visit"
climate.raw$visit <- rep(0, nrow(climate.raw))

for(x in 1:length(site.visit)) {
    
    climate.raw[climate.raw$abmi_site %in% names(site.visit)[x], "visit"] <- as.numeric(site.visit[x])
    
}

rm(site.visit, x)

# Species data, some sites are not found across data sets after accounting for name changes
# 103 sites still do not match up after reclassification
site.reclass <- read.csv("data/lookup/site-reclass_2019-06-14.csv")
species.data <- read.csv("data/base/species/vplants_SiteBinom_2018-05-29.csv")
species.data$site_year <- as.character(species.data$site_year)
species.data$site_year[species.data$site_year %in% site.reclass$species_data] <- as.character(site.reclass$site_data)

# Current landscape, space, and climate 
model.data <- as.data.frame(landscape.summaries$landfacet.site$curr)
model.data["site_year"] <- rownames(model.data)
model.data <- model.data[, -12] # Remove HFor footprint

# Merging of data sets and filter of species with less than 20 occurrences
model.data <- merge.data.frame(model.data, climate.raw, by = "site_year")
model.data <- merge.data.frame(model.data, species.data, by = "site_year")
model.data <- model.data[, c(rep(TRUE, 40), as.logical(colSums(ifelse(model.data[, -c(1:40)] > 0, 1, 0)) > 20))]

# List containing models that will be looped through by the function
habitat.models <- list(as.formula(paste("pcount ~ Cool_slopes + Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Warm_slopes + UrbInd + Crop + TameP + RoughP + SoftLin + HardLin")),
                       as.formula(paste("pcount ~ Cool_slopes + Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Warm_slopes + UrbInd + Crop + I(TameP + RoughP) + SoftLin + HardLin")), 
                       as.formula(paste("pcount ~ Cool_slopes + Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Warm_slopes + UrbInd + I(Crop + TameP + RoughP) + SoftLin + HardLin")), 
                       as.formula(paste("pcount ~ Cool_slopes + I(Dry + Mesic + Wet) + I(Sloped_mesic + Sloped_wet) + Warm_slopes + UrbInd + I(Crop + TameP + RoughP) + SoftLin + HardLin")),
                       as.formula(paste("pcount ~ Cool_slopes + Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Warm_slopes + UrbInd + Crop + TameP + RoughP + SoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Cool_slopes + Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Warm_slopes + UrbInd + Crop + I(TameP + RoughP) + SoftLin + HardLin+ paspen")),
                       as.formula(paste("pcount ~ Cool_slopes + Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Warm_slopes + UrbInd + I(Crop + TameP + RoughP) + SoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Cool_slopes + I(Dry + Mesic + Wet) + I(Sloped_mesic + Sloped_wet) + Warm_slopes + UrbInd + I(Crop + TameP + RoughP) + SoftLin + HardLin + paspen")))

landscape.names <- c("Water", "Cool_slopes", "Dry", 
                     "Mesic", "Sloped_mesic", "Sloped_wet", 
                     "Warm_slopes", "Wet", "Crop", "RoughP", 
                     "TameP", "UrbInd", "SoftLin", "HardLin")

climate.space.names <- c("Intercept", "AHM", "Eref", "FFP", "MAP", "MAT", 
                         "MCMT", "MWMT", "PET", "Lat", "Long", 
                         "Lat2", "Long2", "LatLong", "MAT2", "MWMT2", 
                         "MAPPET", "MAPFFP", "MATAHM")

species.names <- colnames(model.data)[41:350]


species.coef <- list(matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                     matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                     matrix(ncol = length(climate.space.names), nrow = length(species.names), dimnames = list(c(species.names), c(climate.space.names))), 
                     matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))), 
                     matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))))

names(species.coef) <- c("landscape.coef", "landscape.se", "climate.coef", "paspen.coef", "paspen.se")

#for (row.remove in 1:nrow(model.data)) {
for (row.remove in 1:nrow(model.data)) {
        
        # Species modeling
        for (spp in 1:length(species.names)) {
                
                species.coef <- cv_southern_models(data.analysis = model.data[-(row.remove), ], results.store = species.coef, 
                                                landscape.models = habitat.models, prediction.matrix = pm, 
                                                species.ID = species.names[spp])
                
        }
        
        save(species.coef, file = paste("results/boot/facet/land-facet-coef_vplants_removal_", row.remove, "_2019-09-06.Rdata", sep = ""))
        
        # Reset coefficient list
        species.coef <- list(matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                             matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                             matrix(ncol = length(climate.space.names), nrow = length(species.names), dimnames = list(c(species.names), c(climate.space.names))), 
                             matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))), 
                             matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))))
        names(species.coef) <- c("landscape.coef", "landscape.se", "climate.coef", "paspen.coef", "paspen.se")
        
        print(row.remove)
}


##################
# Facet validate # 
##################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load relevant libraries, scripts, and lookups
library(mgcv)  # For binomial GAM
library(mapproj)  # For projected maps
library(binom)  # For exact binomial confidence intervals
library(MuMIn)
library(RcmdrMisc)
library(pROC)
library(arm)

source("src/cross-validation-functions_2019-06-17.R")

# Landscape, climate, space
load("data/processed/facet-soil-proportions_2019-08-16.Rdata")
load("data/processed/spatial-climate_2019-08-16.Rdata")

# Species data, some sites are not found across data sets after accounting for name changes
# 103 sites still do not match up after reclassification
site.reclass <- read.csv("data/lookup/site-reclass_2019-06-14.csv")
species.data <- read.csv("data/base/species/vplants_SiteBinom_2018-05-29.csv")
species.data$site_year <- as.character(species.data$site_year)
species.data$site_year[species.data$site_year %in% site.reclass$species_data] <- as.character(site.reclass$site_data)

# Current landscape, space, and climate 
model.data <- as.data.frame(landscape.summaries$landfacet.site$curr)
model.data["site_year"] <- rownames(model.data)
model.data <- model.data[, -12] # Remove HFor footprint

# Add additional climate variables to the data set
climate.raw["LatLong"] <- climate.raw$Lat * climate.raw$Long
climate.raw["Lat2"] <- (climate.raw$Lat^2)
climate.raw["Long2"] <- (climate.raw$Long^2)
climate.raw["MAT2"] <- climate.raw$MAT * climate.raw$MAT
climate.raw["MWMT2"] <- (climate.raw$MWMT^2)
climate.raw["MAPPET"] <- climate.raw$MAP * climate.raw$PET
climate.raw["MATAHM"] <- climate.raw$MAT * climate.raw$AHM
climate.raw["MAPFFP"] <- climate.raw$MAP * climate.raw$FFP

# Merging of data sets and filter of species with less than 20 occurrences
model.data <- merge.data.frame(model.data, climate.raw, by = "site_year")
model.data <- merge.data.frame(model.data, species.data, by = "site_year")
model.data <- model.data[, c(rep(TRUE, 48), as.logical(colSums(ifelse(model.data[, -c(1:48)] > 0, 1, 0)) > 20))]

rm(climate.raw, species.data, site.reclass, landscape.summaries)

# Make site level predictions and calculate AUC scores for each species

removal.list <- list.files("results/boot/facet/", full.names = TRUE)

for (row.remove in 1:length(removal.list)) {
        
        load(removal.list[row.remove])
        species.list <- rownames(species.coef$landscape.coef)
        
        for (spp in 1:nrow(species.coef$landscape.coef)) {
                
                temp.auc <- cv_site_predictions(model.coef = species.coef,                          
                                             landscape.data = model.data[row.remove, 2:15], 
                                             climate.data = model.data[row.remove, c(19:26, 28, 29, 36:43)], 
                                             paspen.data = model.data[row.remove, "paspen"], 
                                             species.data = model.data[row.remove, c(species.list[spp], "nQuadrant")], 
                                             species.ID = species.list[spp],
                                             landcover.only = TRUE)
                
                if (spp == 1) {
                        
                        site.pred <- temp.auc
                        
                } else {
                        
                        site.pred <- cbind(site.pred, temp.auc)
                        
                }
                
                rm (temp.auc)
        
        }
        
        if (row.remove == 1) {
                
                pred.results <- site.pred
                
        } else {
                
                pred.results <- rbind(pred.results, site.pred)
                
        }
        
        rm(species.coef)
        print(row.remove)
}

save(pred.results, file = "results/boot/pred/land-facet-boot-predictions_vplants_2019-09-09.Rdata")

###############
# Soil models # 
###############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load relevant libraries, scripts, and lookups
library(mgcv)  # For binomial GAM
library(mapproj)  # For projected maps
library(binom)  # For exact binomial confidence intervals
library(MuMIn)
library(RcmdrMisc)
library(pROC)
library(arm)

source("src/cross-validation-functions_2019-06-17.R")

# Load all relevant datasets
# Prediction grid
pm <- read.csv("data/lookup/prediction-matrix/predictionmatrix_soil.csv")
rownames(pm) <- colnames(pm)

# Landscape, climate, space
load("data/processed/facet-soil-proportions_2019-08-16.Rdata")
load("data/processed/spatial-climate_2019-08-16.Rdata")

# Add the number of visits to the climate data
site.visit <- table(droplevels(climate.raw$abmi_site))
colnames(climate.raw)[4] <- "visit"
climate.raw$visit <- rep(0, nrow(climate.raw))

for(x in 1:length(site.visit)) {
    
    climate.raw[climate.raw$abmi_site %in% names(site.visit)[x], "visit"] <- as.numeric(site.visit[x])
    
}

rm(site.visit, x)

# Species data, some sites are not found across data sets after accounting for name changes
# 103 sites still do not match up after reclassification
site.reclass <- read.csv("data/lookup/site-reclass_2019-06-14.csv")
species.data <- read.csv("data/base/species/vplants_SiteBinom_2018-05-29.csv")
species.data$site_year <- as.character(species.data$site_year)
species.data$site_year[species.data$site_year %in% site.reclass$species_data] <- as.character(site.reclass$site_data)

# Current landscape, space, and climate 
model.data <- as.data.frame(landscape.summaries$soil.site$curr)
model.data["site_year"] <- rownames(model.data)
model.data <- model.data[, -9] # Remove HFor footprint

# Merging of data sets and filter of species with less than 20 occurrences
model.data <- merge.data.frame(model.data, climate.raw, by = "site_year")
model.data <- merge.data.frame(model.data, species.data, by = "site_year")
model.data <- model.data[, c(rep(TRUE, 37), as.logical(colSums(ifelse(model.data[, -c(1:37)] > 0, 1, 0)) > 20))]

# List containing models that will be looped through by the function
habitat.models <- list(as.formula(paste("pcount ~ Productive + Clay + Saline + RapidDrain + UrbInd + Crop + TameP + RoughP + SoftLin + HardLin")),
                       as.formula(paste("pcount ~ Productive + Clay + Saline + RapidDrain + UrbInd + Crop + I(TameP + RoughP) + SoftLin + HardLin")), 
                       as.formula(paste("pcount ~ Productive + Clay + Saline + RapidDrain + UrbInd + I(Crop + TameP + RoughP) + SoftLin + HardLin")), 
                       as.formula(paste("pcount ~ Productive + I(Clay + Saline) + RapidDrain + UrbInd + I(Crop + TameP + RoughP) + SoftLin + HardLin")),
                       as.formula(paste("pcount ~ Productive + Clay + Saline + RapidDrain + UrbInd + Crop + TameP + RoughP + SoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Productive + Clay + Saline + RapidDrain + UrbInd + Crop + I(TameP + RoughP) + SoftLin + HardLin+ paspen")),
                       as.formula(paste("pcount ~ Productive + Clay + Saline + RapidDrain + UrbInd + I(Crop + TameP + RoughP) + SoftLin + HardLin + paspen")),
                       as.formula(paste("pcount ~ Productive + I(Clay + Saline) + RapidDrain + UrbInd + I(Crop + TameP + RoughP) + SoftLin + HardLin + paspen")))

landscape.names <- c("Productive", "Clay", "Saline", "RapidDrain", 
                     "Crop", "RoughP", "TameP", "UrbInd", 
                     "SoftLin", "HardLin")


climate.space.names <- c("Intercept", "AHM", "Eref", "FFP", "MAP", "MAT", 
                         "MCMT", "MWMT", "PET", "Lat", "Long", 
                         "Lat2", "Long2", "LatLong", "MAT2", "MWMT2", 
                         "MAPPET", "MAPFFP", "MATAHM")

species.names <- colnames(model.data)[38:347]

species.coef <- list(matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                     matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                     matrix(ncol = length(climate.space.names), nrow = length(species.names), dimnames = list(c(species.names), c(climate.space.names))), 
                     matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))), 
                     matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))))

names(species.coef) <- c("landscape.coef", "landscape.se", "climate.coef", "paspen.coef", "paspen.se")

for (row.remove in 1:nrow(model.data)) {
        
        # Species modeling
        for (spp in 1:length(species.names)) {
                
                species.coef <- cv_southern_models(data.analysis = model.data[-(row.remove), ], results.store = species.coef, 
                                                   landscape.models = habitat.models, prediction.matrix = pm, 
                                                   species.ID = species.names[spp])
                
        }
        
        save(species.coef, file = paste("results/boot/soil/soil-coef_vplants_removal_", row.remove, "_2019-09-06.Rdata", sep = ""))
        
        # Reset coefficient list
        species.coef <- list(matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                             matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                             matrix(ncol = length(climate.space.names), nrow = length(species.names), dimnames = list(c(species.names), c(climate.space.names))), 
                             matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))), 
                             matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))))
        names(species.coef) <- c("landscape.coef", "landscape.se", "climate.coef", "paspen.coef", "paspen.se")
        
        print(row.remove)
        
}


rm(list=ls())
gc()


#################
# Soil validate # 
#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load relevant libraries, scripts, and lookups
library(mgcv)  # For binomial GAM
library(mapproj)  # For projected maps
library(binom)  # For exact binomial confidence intervals
library(MuMIn)
library(RcmdrMisc)
library(pROC)
library(arm)

source("src/cross-validation-functions_2019-06-17.R")

# Landscape, climate, space
load("data/processed/facet-soil-proportions_2019-08-16.Rdata")
load("data/processed/spatial-climate_2019-08-16.Rdata")

# Species data, some sites are not found across data sets after accounting for name changes
# 103 sites still do not match up after reclassification
site.reclass <- read.csv("data/lookup/site-reclass_2019-06-14.csv")
species.data <- read.csv("data/base/species/vplants_SiteBinom_2018-05-29.csv")
species.data$site_year <- as.character(species.data$site_year)
species.data$site_year[species.data$site_year %in% site.reclass$species_data] <- as.character(site.reclass$site_data)

# Current landscape, space, and climate 
model.data <- as.data.frame(landscape.summaries$soil.site$curr)
model.data["site_year"] <- rownames(model.data)
model.data <- model.data[, -9] # Remove HFor footprint

# Add additional climate variables to the data set
climate.raw["LatLong"] <- climate.raw$Lat * climate.raw$Long
climate.raw["Lat2"] <- (climate.raw$Lat^2)
climate.raw["Long2"] <- (climate.raw$Long^2)
climate.raw["MAT2"] <- climate.raw$MAT * climate.raw$MAT
climate.raw["MWMT2"] <- (climate.raw$MWMT^2)
climate.raw["MAPPET"] <- climate.raw$MAP * climate.raw$PET
climate.raw["MATAHM"] <- climate.raw$MAT * climate.raw$AHM
climate.raw["MAPFFP"] <- climate.raw$MAP * climate.raw$FFP

# Merging of data sets and filter of species with less than 20 occurrences
model.data <- merge.data.frame(model.data, climate.raw, by = "site_year")
model.data <- merge.data.frame(model.data, species.data, by = "site_year")
model.data <- model.data[, c(rep(TRUE, 45), as.logical(colSums(ifelse(model.data[, -c(1:45)] > 0, 1, 0)) > 20))]

rm(climate.raw, species.data, site.reclass, landscape.summaries)

# Make site level predictions and calculate AUC scores for each species

removal.list <- list.files("results/boot/soil/", full.names = TRUE)

for (row.remove in 1:nrow(model.data)) {
        
        load(removal.list[row.remove])
        species.list <- rownames(species.coef$landscape.coef)
        
        for (spp in 1:nrow(species.coef$landscape.coef)) {
                
                temp.auc <- cv_site_predictions(model.coef = species.coef,                          
                                                landscape.data = model.data[row.remove, 2:12], 
                                                climate.data = model.data[row.remove, c(16:23, 25, 26, 33:40)], 
                                                paspen.data = model.data[row.remove, "paspen"], 
                                                species.data = model.data[row.remove, c(species.list[spp], "nQuadrant")], 
                                                species.ID = species.list[spp],
                                                landcover.only = TRUE)
                
                if (spp == 1) {
                        
                        site.pred <- temp.auc
                        
                } else {
                        
                        site.pred <- cbind(site.pred, temp.auc)
                        
                }
                
                rm (temp.auc)
                 
                
        }
        
        if (row.remove == 1) {
                
                pred.results <- site.pred
                
        } else {
                
                pred.results <- rbind(pred.results, site.pred)
                
        }
        
        rm(species.coef)
        print(row.remove)
        
}

save(pred.results, file = "results/boot/pred/soil-boot-predictions_vplants_2019-09-09.Rdata")

##############
# Comparison # 
##############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load libraries, scripts, and data sets
library(pROC)

source("src/cross-validation-functions_2019-06-17.R")

load("results/boot/pred/land-facet-boot-predictions_vplants_2019-09-09.Rdata")
facet.predict <- pred.results
load("results/boot/pred/soil-boot-predictions_vplants_2019-09-09.Rdata")
soil.predict <- pred.results
rm(pred.results)

# Calculate AUC from each model set and store in a matrix for comparison

table(colnames(facet.predict) %in% colnames(soil.predict)) # Confirm all colnames match

spp.sequence <- seq(from = 1, to = ncol(facet.predict), by = 2)

for(spp in spp.sequence) {
        
        if(spp == 1) {
                
                species.name <- gsub("_obs", "", colnames(facet.predict)[spp])
                species.auc <- cv_auc(soil.matrix = soil.predict, 
                                      facet.matrix = facet.predict, 
                                      obs = spp, 
                                      pred = (spp + 1), 
                                      species.ID = species.name)
                
        }else {
                
                species.name <- gsub("_obs", "", colnames(facet.predict)[spp])
                temp.auc <- cv_auc(soil.matrix = soil.predict, 
                                   facet.matrix = facet.predict, 
                                   obs = spp, 
                                   pred = (spp + 1), 
                                   species.ID = species.name)
                
                species.auc <- rbind(species.auc, temp.auc)
                rm(temp.auc)
                
                
        }
     
        print(spp)
           
}

save(species.auc, file = "results/coef/soil-facet-only-boot-AUC_2019-09-11.Rdata")

# Assess variation between the different approaches
cor(species.auc[, 1:2]) # They are 83.3% correlated

png(filename = "results/figures/facet-soil-cv-AUC_2019-09-09.png",
    height = 1200,
    width = 1200)

plot(species.auc[, 1] ~ species.auc[, 2], 
     xlab = "Soil", 
     ylab = "Land facet", 
     col = ifelse(species.auc[, 3] < 0.05 & species.auc[, 1] > species.auc[, 2], "red", 
                  ifelse(species.auc[, 3] < 0.05 & species.auc[, 1] < species.auc[, 2], "blue","black")),
     pch = 15,
     xlim = c(0.5,0.9),
     ylim = c(0.5, 0.9),
     cex = 2)
abline(0,1)
dev.off()

# table(species.auc[, 3] < 0.05 & species.auc[, 1] > species.auc[, 2]) # 117 improved
# table(species.auc[, 3] < 0.05 & species.auc[, 1] < species.auc[, 2]) # 36 worsened


