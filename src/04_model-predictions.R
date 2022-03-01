#
# Title: Habitat model predictions
# Created: June 12th, 2019
# Last Updated: August 19th, 2019
# Author: Brandon Allen
# Objectives: Creating predictions for both sets of species habitat models.
# Keywords: Facet site, Facet grid, Soil site, Soil grid
#

##############
# Facet site # 
##############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

source("src/habitat-models-functions_2019-06-12.R")

# Landscape, climate, space
load("data/processed/facet-soil-proportions_2019-08-16.Rdata")
load("data/processed/spatial-climate_2019-08-16.Rdata")
load("results/coef/land-facet-coefficients_vplants_2019-08-16.Rdata")

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
model.data <- model.data[, c(rep(TRUE, 49), as.logical(colSums(ifelse(model.data[, -c(1:49)] > 0, 1, 0)) > 20))]

rm(climate.raw, species.data, site.reclass, landscape.summaries)

# # Make site level predictions and calculate AUC scores for each species
# 
# species.list <- rownames(species.coef$landscape.coef)
# 
# for (spp in 1:nrow(species.coef$landscape.coef)) {
#         
#         temp.auc <- site_predictions(model.coef = species.coef,                          
#                                      landscape.data = model.data[, 2:16], 
#                                      climate.data = model.data[, c(20:27, 29, 30, 37:44)], 
#                                      paspen.data = model.data$paspen, 
#                                      species.data = model.data[, c(species.list[spp], "nQuadrant")], 
#                                      species.ID = species.list[spp])
#         
#         if (spp == 1) {
#                 
#                 facet.auc <- temp.auc
#                 
#         } else {
#                 
#                 facet.auc <- cbind(facet.auc, temp.auc)
#                 
#         }
# 
#         rm (temp.auc)
#         print(spp)   
#         
# }

# Make site level predictions and calculate AUC scores for each species

species.list <- rownames(species.coef$landscape.coef)

for (spp in 1:nrow(species.coef$landscape.coef)) {

        temp.auc <- landscape_only_predictions(model.coef = species.coef,
                                     landscape.data = model.data[, c(2:15)],
                                     species.data = model.data[, c(species.list[spp], "nQuadrant")],
                                     species.ID = species.list[spp])

        if (spp == 1) {

                facet.auc <- temp.auc

        } else {

                facet.auc <- cbind(facet.auc, temp.auc)

        }

        rm (temp.auc)
        print(spp)

}

save(facet.auc, file = "results/coef/land-facet-only-auc-vplants_2019-08-16.Rdata")

##############
# Facet grid # 
##############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load relevant libraries, scripts, and lookups
library(arm)
library(mgcv)  # For binomial GAM
library(mapproj)  # For projected maps
library(binom)  # For exact binomial confidence intervals
library(MuMIn)
library(RcmdrMisc)
library(RColorBrewer)
library(pROC)
library(scales)

source("src/habitat-models-functions_2019-06-12.R")

# Landscape, climate, space
load("data/processed/facet-kgrid_2019-08-16.Rdata")
load("data/base/habitat/kgrid_table_km.Rdata")
load("results/coef/land-facet-coefficients_vplants_2019-08-16.Rdata")

# Current landscape, space, and climate 
colnames(kgrid)[4:6] <- c("LinkID", "Long", "Lat")

# Add additional climate variables to the data set
kgrid["LatLong"] <- kgrid$Lat * kgrid$Long
kgrid["Lat2"] <- (kgrid$Lat^2)
kgrid["Long2"] <- (kgrid$Long^2)
kgrid["MAT2"] <- kgrid$MAT * kgrid$MAT
kgrid["MWMT2"] <- (kgrid$MWMT^2)
kgrid["MAPPET"] <- kgrid$MAP * kgrid$PET
kgrid["MATAHM"] <- kgrid$MAT * kgrid$AHM
kgrid["MAPFFP"] <- kgrid$MAP * kgrid$FFP
kgrid <- kgrid[, c(4:6, 10:17, 36:43, 7, 8)]

# Merging of data sets and filter of species with less than 20 occurrences
km2.info <- merge.data.frame(facet.grid, kgrid, by = "LinkID")
rm(facet.grid, kgrid)

# Make site level predictions and calculate AUC scores for each species
species.list <- rownames(species.coef$landscape.coef)
col.pal <- colorRampPalette(c("Blue", "White", "Orange", "Red"))
col.mask <- !(km2.info$NRNAME %in% c("Grassland", "Parkland") | (km2.info$NSRNAME %in% c("Dry Mixedwood")))

for (spp in 1:nrow(species.coef$landscape.coef)) {
        
        temp.predictions <- grid_predictions(model.coef = species.coef,                          
                                     landscape.data = km2.info[, 2:9], 
                                     climate.data = km2.info[, c(10:18, 20:27)], 
                                     paspen.data = km2.info$pAspen,
                                     species.ID = species.list[spp])
        temp.predictions[col.mask] <- 0
        temp.predictions <- rescale(temp.predictions, to = c(1,100))
        temp.colors <- col.pal(100)[temp.predictions]
        temp.colors[col.mask] <- "Grey"
        
        png(filename = paste("results/figures/facet/reference/", species.list[spp], ".png", sep = ""),
            height = 1200,
            width = 800)
        
        plot(km2.info$Long, km2.info$Lat, pch = 15, cex = 0.25, col = temp.colors)

        dev.off()
        
        rm(temp.predictions, temp.colors)
        print(spp)   
        
}

#############
# Soil site # 
#############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

source("src/habitat-models-functions_2019-06-12.R")

# Landscape, climate, space
load("data/processed/facet-soil-proportions_2019-08-16.Rdata")
load("data/processed/spatial-climate_2019-08-16.Rdata")
load("results/coef/soil-coefficients_vplants_2019-08-16.Rdata")

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

# # Make site level predictions and calculate AUC scores for each species
# 
# species.list <- rownames(species.coef$landscape.coef)
# 
# for (spp in 1:nrow(species.coef$landscape.coef)) {
#         
#         temp.auc <- site_predictions(model.coef = species.coef,                          
#                                      landscape.data = model.data[, 2:12], 
#                                      climate.data = model.data[, c(16:23, 25, 26, 33:40)], 
#                                      paspen.data = model.data$paspen, 
#                                      species.data = model.data[, c(species.list[spp], "nQuadrant")], 
#                                      species.ID = species.list[spp])
#         
#         if (spp == 1) {
#                 
#                 soil.auc <- temp.auc
#                 
#         } else {
#                 
#                 soil.auc <- cbind(soil.auc, temp.auc)
#                 
#         }
#         
#         rm (temp.auc)
#         print(spp)   
#         
# }

# Make site level predictions and calculate AUC scores for each species

species.list <- rownames(species.coef$landscape.coef)

for (spp in 1:nrow(species.coef$landscape.coef)) {
    
    temp.auc <- landscape_only_predictions(model.coef = species.coef,                          
                                 landscape.data = model.data[, c(2:5, 7:12)],
                                 species.data = model.data[, c(species.list[spp], "nQuadrant")], 
                                 species.ID = species.list[spp])
    
    if (spp == 1) {
        
        soil.auc <- temp.auc
        
    } else {
        
        soil.auc <- cbind(soil.auc, temp.auc)
        
    }
    
    rm (temp.auc)
    print(spp)   
    
}

save(soil.auc, file = "results/coef/soil-only-auc-vplants_2019-08-16.Rdata")

#############
# Soil grid # 
#############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load relevant libraries, scripts, and lookups
library(arm)
library(mgcv)  # For binomial GAM
library(mapproj)  # For projected maps
library(binom)  # For exact binomial confidence intervals
library(MuMIn)
library(RcmdrMisc)
library(RColorBrewer)
library(pROC)
library(scales)

source("src/habitat-models-functions_2019-06-12.R")

# Landscape, climate, space
load("data/processed/soil-kgrid_2019-08-16.Rdata")
load("data/base/habitat/kgrid_table_km.Rdata")
load("results/coef/soil-coefficients_vplants_2019-08-16.Rdata")

# Current landscape, space, and climate 
colnames(kgrid)[4:6] <- c("LinkID", "Long", "Lat")

# Add additional climate variables to the data set
kgrid["LatLong"] <- kgrid$Lat * kgrid$Long
kgrid["Lat2"] <- (kgrid$Lat^2)
kgrid["Long2"] <- (kgrid$Long^2)
kgrid["MAT2"] <- kgrid$MAT * kgrid$MAT
kgrid["MWMT2"] <- (kgrid$MWMT^2)
kgrid["MAPPET"] <- kgrid$MAP * kgrid$PET
kgrid["MATAHM"] <- kgrid$MAT * kgrid$AHM
kgrid["MAPFFP"] <- kgrid$MAP * kgrid$FFP
kgrid <- kgrid[, c(4:6, 10:17, 36:43, 7, 8)]

# Merging of data sets and filter of species with less than 20 occurrences
km2.info <- merge.data.frame(soil.grid, kgrid, by = "LinkID")
rm(soil.grid, kgrid)

# Make site level predictions and calculate AUC scores for each species
species.list <- rownames(species.coef$landscape.coef)
col.pal <- colorRampPalette(c("Blue", "White", "Orange", "Red"))
col.mask <- !(km2.info$NRNAME %in% c("Grassland", "Parkland") | (km2.info$NSRNAME %in% c("Dry Mixedwood")))

for (spp in 1:nrow(species.coef$landscape.coef)) {
        
        temp.predictions <- grid_predictions(model.coef = species.coef,                          
                                             landscape.data = km2.info[, c(2:5)], 
                                             climate.data = km2.info[, c(8:16, 18:25)], 
                                             paspen.data = km2.info$pAspen,
                                             species.ID = species.list[spp])
        temp.predictions[col.mask] <- 0
        temp.predictions <- rescale(temp.predictions, to = c(1,100))
        temp.colors <- col.pal(100)[temp.predictions]
        temp.colors[col.mask] <- "Grey"
        
        png(filename = paste("results/figures/soil/reference/", species.list[spp], ".png", sep = ""),
            height = 1200,
            width = 800)
        
        plot(km2.info$Long, km2.info$Lat, pch = 15, cex = 0.25, col = temp.colors)
        
        dev.off()
        
        rm(temp.predictions, temp.colors)
        print(spp)   
        
}

##############
# Comparison # 
##############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load libraries, scripts, and data sets
library(pROC)

source("src/habitat-models-functions_2019-06-12.R")

load("results/coef/soil-only-auc-vplants_2019-08-16.Rdata")
load("results/coef/land-facet-only-auc-vplants_2019-08-16.Rdata")

# Calculate AUC from each model set and store in a matrix for comparison

table(colnames(facet.auc) %in% colnames(soil.auc)) # Confirm all colnames match

spp.sequence <- seq(from = 1, to = ncol(facet.auc), by = 2)

for(spp in spp.sequence) {
        
        if(spp == 1) {
                
                species.name <- gsub("_obs", "", colnames(facet.auc)[spp])
                species.auc <- model_auc(soil.matrix = soil.auc, 
                                      facet.matrix = facet.auc, 
                                      obs = spp, 
                                      pred = (spp + 1), 
                                      species.ID = species.name)
                
        }else {
                
                species.name <- gsub("_obs", "", colnames(facet.auc)[spp])
                temp.auc <- model_auc(soil.matrix = soil.auc, 
                                   facet.matrix = facet.auc, 
                                   obs = spp, 
                                   pred = (spp + 1), 
                                   species.ID = species.name)
                
                species.auc <- rbind(species.auc, temp.auc)
                rm(temp.auc)
                
                
        }
        
        print(spp)
        
}

save(species.auc, file = "results/coef/soil-facet-only-AUC_2019-09-11.Rdata")

# Assess variation between the different approaches
cor(species.auc[, 1:2]) # They are 60.7% correlated

png(filename = "results/figures/facet-soil-only-model-AUC_2019-09-04.png",
    height = 1200,
    width = 1200)

plot(species.auc[, 1] ~ species.auc[, 2], 
     xlab = "Soil", 
     ylab = "Land facet",
     xlim = c(0.4,1),
     ylim = c(0.4,1), 
     col = ifelse(species.auc[, 3] < 0.05 & species.auc[, 1] > species.auc[, 2], "red", 
                  ifelse(species.auc[, 3] < 0.05 & species.auc[, 1] < species.auc[, 2], "blue","black")),
     pch = 15,
     cex = 2)
abline(0,1)
dev.off()

# Cluster analysis Facet
load("results/coef/land-facet-coefficients_vplants_2019-08-16.Rdata")

# Perform PCA on the diversity variability
pca.out1 <- prcomp(species.coef$landscape.coef, scale = TRUE, center = TRUE)
summary(pca.out1) # PCA loadinds are much better not having the life strategy separated. 96.7% in the first four axis
plot(pca.out1$rotation, pch = 15, xlim = c(-0.1, 0.4), ylim = c(-0.5, 0.5))
text(pca.out1$rotation, labels = rownames(pca.out1$rotation), pos = 2)
dev.off()

# Cluster analysis Soil
load("results/coef/soil-coefficients_vplants_2019-08-16.Rdata")

# Perform PCA on the diversity variability
pca.out1 <- prcomp(species.coef$landscape.coef, scale = TRUE, center = TRUE)
summary(pca.out1) 
png(file = "results/figures/soil-pca_2019-09-04.png",
    height = 800,
    width = 800)
plot(pca.out1$rotation, pch = 15, xlim = c(0.1, 0.4), ylim = c(-0.5, 0.5),
     xlab = "PC1 53.5%",
     ylab = "PC2 21.5%")
text(pca.out1$rotation, labels = rownames(pca.out1$rotation), pos = 1)
dev.off()

# Cluster analysis Facet
load("results/coef/land-facet-coefficients_vplants_2019-08-16.Rdata")

# Perform PCA on the diversity variability
pca.out1 <- prcomp(species.coef$landscape.coef, scale = TRUE, center = TRUE)
summary(pca.out1) 
png(file = "results/figures/facet-pca_2019-09-04.png",
    height = 800,
    width = 800)
plot(pca.out1$rotation, pch = 15, xlim = c(-0.1, 0.4), ylim = c(-0.5, 0.5),
     xlab = "PC1 38.6%",
     ylab = "PC2 18.4%")
text(pca.out1$rotation, labels = rownames(pca.out1$rotation), pos = 2)
dev.off()





# Native visualization

vplant.pca <- as.data.frame(pca.out1$x)
vplant.pca["sppid"] <- rownames(vplant.pca)

spp.lookup <- read.csv("data/lookup/vplants-modeled_lookup.csv")

vplant.pca <- merge.data.frame(vplant.pca, spp.lookup, by = "sppid")

plot(x = vplant.pca$PC1, y = vplant.pca$PC2, 
     pch = 15,
     col = ifelse(vplant.pca$Origin == "Exotic", "red", "blue"))


