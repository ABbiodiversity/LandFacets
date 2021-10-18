#
# Title: Species habitat models 
# Created: June 12th, 2019
# Last Updated: October 8th, 2021
# Author: Brandon Allen
# Objectives: Species habitat modeling for both land facet and soil data
# Keywords: Facet models , Soil models 
#

################
# Facet models # MODLES ARE MOSTLY CORRECT, BUT SOME SPECIES AUC DECREASES WITH THE FULL MODEL. INVESTIGATE!!
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

source("src/habitat-models-functions_2019-06-12.R")

# Load all relevant datasets
# Prediction grid
pm <- read.csv("data/lookup/prediction-matrix/predictionmatrix_terrain-water.csv")
rownames(pm) <- pm$VegType
pm <- pm[, -1]

# Landscape, climate, space
load("data/processed/facet-soil-proportions_2021-10-07.Rdata")
load("data/processed/spatial-climate_2021-10-07.Rdata")

# Add the number of visits to the climate data
site.visit <- table(droplevels(climate.raw$abmi_site))
colnames(climate.raw)[4] <- "visit"
climate.raw$visit <- rep(0, nrow(climate.raw))

for(x in 1:length(site.visit)) {
        
        climate.raw[climate.raw$abmi_site %in% names(site.visit)[x], "visit"] <- as.numeric(site.visit[x])
        
}

rm(site.visit, x)

# Current landscape, space, and climate 
site.data <- as.data.frame(landscape.summaries$landfacet.site$curr)
site.data["site_year"] <- rownames(site.data)
site.data <- site.data[, !(colnames(site.data) %in% c("HFor", "HWater"))] # Remove unmodeld HF

# Add the new categories
for (col.id in colnames(pm)) {
        
        # Identify coefficient set
        coef.set <- rownames(pm)[as.logical(pm[, col.id])]
        
        if (length(coef.set) == 1) {
                
                site.data[, col.id] <- site.data[, colnames(site.data) %in% coef.set]
                
        } else {
                
                site.data[, col.id] <- rowSums(site.data[, colnames(site.data) %in% coef.set])
                
        }
        
}

# Define species occurrence files
occurrence.list <- list.files("data/processed/", full.names = TRUE)

for (taxon in c("lichen", "mite", "moss", "vplant")) {
        
        # Load the species data (vascular plants test)
        load(occurrence.list[grep(paste0(taxon, "-site-occurrence"), occurrence.list)])
        
        # Merging of data sets and filter of species with less than 20 occurrences
        model.data <- merge.data.frame(site.data, climate.raw, by = "site_year")
        model.data <- merge.data.frame(model.data, occurrence.in, by = "site_year")
        model.data <- model.data[, c(rep(TRUE, 59), as.logical(colSums(ifelse(model.data[, -c(1:59)] > 0, 1, 0)) > 20))]
        
        # List containing models that will be looped through by the function
        habitat.models <- list(as.formula(paste("pcount ~ Cool_slopes + Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Warm_slopes + UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin")),
                               as.formula(paste("pcount ~ Cool_slopes + Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Warm_slopes + UrbInd + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin")),
                               as.formula(paste("pcount ~ Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Sloped_Dry + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin")),
                               as.formula(paste("pcount ~ Dry + Mesic + Wet + Sloped_Moist + Sloped_Dry + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin")),
                               as.formula(paste("pcount ~ Dry + Moist + Sloped_Dry + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin")),
                               as.formula(paste("pcount ~ Dry + Moist + Sloped_Dry + UrbIndWellsites + Rural + Cult + SoftLin + HardLin")),
                               as.formula(paste("pcount ~ Dry + Moist + Sloped_Dry + Alien + SoftLin")),
                               as.formula(paste("pcount ~ Cool_slopes + Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Warm_slopes + UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Cool_slopes + Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Warm_slopes + UrbInd + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Dry + Mesic + Wet + Sloped_mesic + Sloped_wet + Sloped_Dry + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Dry + Mesic + Wet + Sloped_Moist + Sloped_Dry + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Dry + Moist + Sloped_Dry + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Dry + Moist + Sloped_Dry + UrbIndWellsites + Rural + Cult + SoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Dry + Moist + Sloped_Dry + Alien + SoftLin + paspen")))
        
        landscape.names <- c("Cool_slopes", "Dry", 
                             "Mesic", "Sloped_mesic", "Sloped_wet", 
                             "Warm_slopes", "Wet", "EnSeismic", "EnSoftLin", 
                             "TrSoftLin", "HardLin", "UrbInd", "Rural",
                             "Wellsites", "Crop", "TameP", "RoughP")
        
        climate.space.names <- c("Intercept", "AHM", "Eref", "FFP", "MAP", "MAT", 
                                 "MCMT", "MWMT", "PET", "Lat", "Long", 
                                 "Lat2", "Long2", "LatLong", "MAT2", "MWMT2", 
                                 "MAPPET", "MAPFFP", "MATAHM")
        
        species.names <- colnames(model.data)[60:ncol(model.data)]
        
        species.coef <- list(matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                             matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                             matrix(ncol = length(climate.space.names), nrow = length(species.names), dimnames = list(c(species.names), c(climate.space.names))), 
                             matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))), 
                             matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))),
                             matrix(ncol = 4, nrow = length(species.names), dimnames = list(c(species.names), c("auc_LC", "auc_both", "dect", "survey"))))
        
        names(species.coef) <- c("landscape.coef", "landscape.se", "climate.coef", "paspen.coef", "paspen.se", "fit")
        
        # Species modeling
        for (spp in 1:length(species.names)) {
                
                species.coef <- southern_models(data.analysis = model.data, results.store = species.coef, 
                                                landscape.models = habitat.models, prediction.matrix = pm, 
                                                species.ID = species.names[spp])
                print(paste0(taxon, " ", spp, "/", length(species.names)))
                
        }
        
        save(species.coef, file = paste0("results/coef/land-facet-coefficients_", taxon, "_2021-10-08.Rdata"))
        
}



###############
# Soil models # NEED TO RESOLVE THE SOIL MODELS
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

source("src/habitat-models-functions_2019-06-12.R")

# Load all relevant datasets
# Prediction grid
pm <- read.csv("data/lookup/prediction-matrix/soil-prediction-matrix_2020.csv")
rownames(pm) <- pm$VegType
pm <- pm[, -1]

# Landscape, climate, space
load("data/processed/facet-soil-proportions_2021-10-07.Rdata")
load("data/processed/spatial-climate_2021-10-07.Rdata")

# Add the number of visits to the climate data
site.visit <- table(droplevels(climate.raw$abmi_site))
colnames(climate.raw)[4] <- "visit"
climate.raw$visit <- rep(0, nrow(climate.raw))

for(x in 1:length(site.visit)) {
        
        climate.raw[climate.raw$abmi_site %in% names(site.visit)[x], "visit"] <- as.numeric(site.visit[x])
        
}

rm(site.visit, x)

# Current landscape, space, and climate 
site.data <- as.data.frame(landscape.summaries$soil.site$curr)
site.data["site_year"] <- rownames(site.data)
site.data <- site.data[, !(colnames(site.data) %in% c("HFor", "HWater"))] # Remove unmodeld HF

# Add the new categories
for (col.id in colnames(pm)) {
        
        # Identify coefficient set
        coef.set <- rownames(pm)[as.logical(pm[, col.id])]
        
        if (length(coef.set) == 1) {
                
                site.data[, col.id] <- site.data[, colnames(site.data) %in% coef.set]
                
        } else {
                
                site.data[, col.id] <- rowSums(site.data[, colnames(site.data) %in% coef.set])
                
        }
        
}

# Define species occurrence files
occurrence.list <- list.files("data/processed/", full.names = TRUE)

for (taxon in c("lichen", "mite", "moss", "vplant")) {
        
        # Load the species data (vascular plants test)
        load(occurrence.list[grep(paste0(taxon, "-site-occurrence"), occurrence.list)])
        
        # Merging of data sets and filter of species with less than 20 occurrences
        model.data <- merge.data.frame(site.data, climate.raw, by = "site_year")
        model.data <- merge.data.frame(model.data, occurrence.in, by = "site_year")
        model.data <- model.data[, c(rep(TRUE, 59), as.logical(colSums(ifelse(model.data[, -c(1:59)] > 0, 1, 0)) > 20))]
        
        # List containing models that will be looped through by the function
        habitat.models <- list(as.formula(paste("pcount ~ Blowout + ClaySub + Loamy + RapidDrain + SandyLoam + ThinBreak + Other + UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin")),
                               as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid + Other +  UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin")),
                               as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid + Other +  UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin")),
                               as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin")),
                               as.formula(paste("pcount ~ Blowout + ClaySub + Loamy + RapidDrain + SandyLoam + ThinBreak  + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin")),
                               as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin")),
                               as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin")),
                               as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbIndWellsites + Rural + Cult + SoftLin + HardLin")),
                               as.formula(paste("pcount ~ Productive + Nonproductive + Other + Alien + SoftLin")),
                               as.formula(paste("pcount ~ Blowout + ClaySub + Loamy + RapidDrain + SandyLoam + ThinBreak  + Other + UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid +  Other +  UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid +  Other +  UrbInd + Rural + Wellsites + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbInd + Wellsites + Rural + Crop + TameP + RoughP + EnSoftLin + EnSeismic + TrSoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Blowout + ClaySub + Loamy +RapidDrain + SandyLoam + ThinBreak  + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Blowout + ClaySubThin + Loamy + SandyRapid + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbIndWellsites + Rural + Crop + Pasture + SoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Productive + Nonproductive + Other + UrbIndWellsites + Rural + Cult + SoftLin + HardLin + paspen")),
                               as.formula(paste("pcount ~ Productive + Nonproductive + Other + Alien + SoftLin + paspen")))
        
        landscape.names <- c("Loamy", "SandyLoam", 
                             "ClaySub", "RapidDrain", "Blowout", 
                             "ThinBreak", "Other", "EnSeismic", "EnSoftLin", 
                             "TrSoftLin", "HardLin", "UrbInd", "Rural",
                             "Wellsites", "Crop", "TameP", "RoughP")
        
        climate.space.names <- c("Intercept", "AHM", "Eref", "FFP", "MAP", "MAT", 
                                 "MCMT", "MWMT", "PET", "Lat", "Long", 
                                 "Lat2", "Long2", "LatLong", "MAT2", "MWMT2", 
                                 "MAPPET", "MAPFFP", "MATAHM")
        
        species.names <- colnames(model.data)[60:ncol(model.data)]
        
        species.coef <- list(matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                             matrix(ncol = length(landscape.names), nrow = length(species.names), dimnames = list(c(species.names), c(landscape.names))), 
                             matrix(ncol = length(climate.space.names), nrow = length(species.names), dimnames = list(c(species.names), c(climate.space.names))), 
                             matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))), 
                             matrix(ncol = 1, nrow = length(species.names), dimnames = list(c(species.names), c("paspen"))),
                             matrix(ncol = 4, nrow = length(species.names), dimnames = list(c(species.names), c("auc_LC", "auc_both", "dect", "survey"))))
        
        names(species.coef) <- c("landscape.coef", "landscape.se", "climate.coef", "paspen.coef", "paspen.se", "fit")
        
        # Species modeling
        for (spp in 1:length(species.names)) {
                
                species.coef <- southern_models(data.analysis = model.data, results.store = species.coef, 
                                                landscape.models = habitat.models, prediction.matrix = pm, 
                                                species.ID = species.names[spp])
                print(paste0(taxon, " ", spp, "/", length(species.names)))
                
        }
        
        save(species.coef, file = paste0("results/coef/soil-coefficients_", taxon, "_2021-10-08.Rdata"))
        
}

rm(list=ls())
gc()

##########################
# Exploration of results #
##########################

for(taxon in c("lichen", "mite", "moss", "vplant")) {
        
        load(paste0("results/coef/soil-coefficients_", taxon, "_2021-10-08.Rdata"))
        soil.results <- species.coef$fit
        soil.results <- data.frame(Species = rownames(soil.results),
                                   Soil_AUC = soil.results[,2])
        load(paste0("results/coef/land-facet-coefficients_", taxon, "_2021-10-08.Rdata"))
        facet.results <- species.coef$fit
        facet.results <- data.frame(Species = rownames(facet.results),
                                   Facet_AUC = facet.results[,2])
        
        combined.results <- merge.data.frame(soil.results, facet.results, by = "Species")
        
        # Visualize
        print(ggplot(data = combined.results) +
                      geom_point(aes(x = Soil_AUC, y = Facet_AUC, color = "Area"), show.legend = FALSE) +
                      geom_abline(intercept = 0) +
                      scale_color_manual(values = c("#2D415B")) +
                      xlab("Soil model AUC") +
                      ylab("Facet model AUC") + 
                      ggtitle(paste(taxon, "; proportion improved = ", round(as.numeric(table(combined.results$Facet_AUC > combined.results$Soil_AUC)["TRUE"]/nrow(combined.results)),3))) +
                      theme_light())
        
}
