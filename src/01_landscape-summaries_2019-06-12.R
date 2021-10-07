#
# Title: Site summaries for both soil and land facet layers
# Created: June 12th, 2019
# Last Updated: October 7th, 2021
# Author: Brandon Allen
# Objectives: Create site and quadrant level summaries of the soil and land facet layers
# Keywords: Initialization, Long form, Site proportions, Kgrid proportions
# Notes: 1) Updated the soil and human footprint lookup tables to the 2021 standard
#

##################
# Initialization # 
##################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load relevant libraries, scripts, and lookups
library(DBI)
library(ggplot2)
library(mefa4)
library(rgdal)
library(RSQLite)

source("src/landscape-summary-functions_2019-06-13.R")
soil.facet.hf.lookup <- read.csv("data/lookup/soil-facet-hf-lookup-v61_2021.csv")

# Load and clean land facet data 
f <- file.path("data/base/habitat/20190318_SummaryTables_1ha_TerrestrialSites_Veg61_vHF_LandFacets_SurveyYear_2003_2018.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
landscape.raw <- dbReadTable(db, "Summary_1ha")
dbDisconnect(db)
rm(db, f)

# Remove duplicate rows and sites where there is no soil information available
landscape.raw <- landscape.raw[!duplicated(landscape.raw[, !(colnames(landscape.raw) %in% "year")]),]
landscape.raw <- landscape.raw[!(landscape.raw$UID %in% unique(landscape.raw$UID[is.na(landscape.raw$Soil_Type_1)])), ]
topo.continuous <- landscape.raw[, c("UID", "cti", "solar_local", "solar_lat", "aspect_deg", "slope_deg")]

landscape.raw <- data.frame(quadrant = landscape.raw$Section,
                            site_year = landscape.raw$UID,
                            shape_area = landscape.raw$Shape_Area,
                            nr = landscape.raw$NRNAME,
                            nsr = landscape.raw$NSRNAME,
                            soil_class = factor(landscape.raw$Soil_Type_1),
                            feature_ty = factor(landscape.raw$FEATURE_TY),
                            cti_123 = landscape.raw$cti_123,
                            solar_123 = landscape.raw$solar_123,
                            water_01 = landscape.raw$water_01,
                            snowice_01 = landscape.raw$snowice_01,
                            springs_01 = landscape.raw$springs_01,
                            dunes_01 = landscape.raw$dunes_01)

levels(landscape.raw$feature_ty) <- c(levels(landscape.raw$feature_ty), "NATIVE")
landscape.raw$feature_ty[is.na(landscape.raw$feature_ty)] <- as.factor("NATIVE")

# Load and clean climate and spatial data
load("data/base/habitat/veg-hf_SiteCenter_v6verified.Rdata")
xx <- xx[xx$NATURAL_REGIONS %in% c("Grassland", "Parkland") | xx$NATURAL_SUBREGIONS == "Dry Mixedwood", ]

climate.raw <- data.frame(abmi_site = xx$ABMI_Assigned_Site_ID,
                          year = xx$survey_year,
                          site_year = xx$Site_YEAR,
                          on_off = xx$On_Off,
                          AHM = xx$AHM,
                          Eref = xx$Eref,
                          FFP = xx$FFP,
                          MAP = xx$MAP,
                          MAT = xx$MAT,
                          MCMT = xx$MCMT,
                          MWMT = xx$MWMT,
                          PET = xx$PET,
                          paspen = xx$pAspen,
                          Lat = xx$PUBLIC_LATTITUDE,
                          Long = xx$PUBLIC_LONGITUDE,
                          Elevation = xx$ELEVATION)

rm(dd_150m, dd_1ha, dd_564m, dd_point, dd_qha, dw_1ha, dw_qha, xx)

# Filter the climate and landscape data to include only matching sites
landscape.raw <- landscape.raw[landscape.raw$site_year %in% climate.raw$site_year, ]
climate.raw <- climate.raw[climate.raw$site_year %in% unique(landscape.raw$site_year), ]

# Aggregate the continuous topographic variables and bind to the climate data
name.store <- colnames(topo.continuous)
topo.continuous <- aggregate(topo.continuous[, 2:6], by = list(topo.continuous[, "UID"]), FUN = mean)
colnames(topo.continuous) <- c("site_year", name.store[-1])
topo.continuous$site_year <- as.factor(topo.continuous$site_year)
climate.raw <- merge.data.frame(climate.raw, topo.continuous, by = "site_year")

save(climate.raw, file = "data/processed/spatial-climate_2021-10-07.Rdata")

rm(climate.raw, topo.continuous, name.store)

#############
# Long form # 
#############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Definitions of sub-model cetegories relevant in creating the Terrain x Water x Dunes land facet layer
# cti_123:
# 1 = dry, 2 = mesic, 3 = wet
# solar_123: 
# 1 = cool slopes, 2 = normal, 3 = warm slopes
# springs_01:
# 1 = springs
# snowice_01:
# 1 = snow/ice
# water_01:
# 1 = water
# dunes_01:
# 1 = dunes

# Rule set for creating the Terrain land facet layer
# cti_123 x solar_123 = Terrain
# 1_1 = dry cool slope = Cool slopes 
# 1_2 = dry normal = Dry
# 1_3 = dry warm slope = Warm slopes
# 2_1 = mesic cool slope = Sloped mesic
# 2_2 = mesic normal = Mesic
# 2_3 = mesic warm slope = Sloped mesic
# 3_1 = wet cool sloped = Sloped wet
# 3_2 = wet normal flat = Wet
# 3_3 = wet warm slopes = Sloped wet

# Create Terrain layer
landscape.raw["Terrain"] <- paste(landscape.raw$cti_123, landscape.raw$solar_123, sep = "_")
landscape.raw$Terrain[landscape.raw$Terrain == "1_1"] <- "Cool_slopes"
landscape.raw$Terrain[landscape.raw$Terrain == "1_2"] <- "Dry"
landscape.raw$Terrain[landscape.raw$Terrain == "1_3"] <- "Warm_slopes"
landscape.raw$Terrain[landscape.raw$Terrain == "2_1"] <- "Sloped_mesic"
landscape.raw$Terrain[landscape.raw$Terrain == "2_2"] <- "Mesic"
landscape.raw$Terrain[landscape.raw$Terrain == "2_3"] <- "Sloped_mesic"
landscape.raw$Terrain[landscape.raw$Terrain == "3_1"] <- "Sloped_wet"
landscape.raw$Terrain[landscape.raw$Terrain == "3_2"] <- "Wet"
landscape.raw$Terrain[landscape.raw$Terrain == "3_3"] <- "Sloped_wet"

# Rule set for creating the Terrain x Water land facet layer 
# Terrain x water_01 x springs_01 x snowice_01 = Terrain
# If any of the three water groups = 1, replace value from the Terrain model

# Create land_facet layer
landscape.raw["land_facet"] <- landscape.raw$Terrain
landscape.raw$land_facet[landscape.raw$water_01 == 1] <- "Water"
landscape.raw$land_facet[landscape.raw$snowice_01 == 1] <- "Snowice"
landscape.raw$land_facet[landscape.raw$springs_01 == 1] <- "Springs"

# Factor the land facet and soil categories
landscape.raw$land_facet <- factor(landscape.raw$land_facet)
landscape.raw$soil_class <- factor(landscape.raw$soil_class)

# In this analysis, we are going to exclude Dunes from the layer set as it is more dificult to perform the summaries myself.
# If Eric has time, it would be good to perform the summaries for the 1km2 grid.
# There are only 17 sites with Dunes in the vascular plant data. However, the proportion of the site
# ranges between 0 and 80 so it has a good distribution. 
# # Rule set for creating the Terrain x Water x Dunes land facet layer
# # Terrain x Dunes = land_facet
# # If dunes = 1, and the Terrain model = Cool slopes, Warm slopes, Sloped mesic, Dry, or Mesic, replace value from the Terrain model 
# 
# # Create land_facet layer
# landscape.raw["land_facet"] <- landscape.raw$Terrain_Water
# landscape.raw$land_facet[landscape.raw$dunes_01 == 1 & landscape.raw$land_facet == "Cool_slopes"] <- "Dunes"
# landscape.raw$land_facet[landscape.raw$dunes_01 == 1 & landscape.raw$land_facet == "Warm_slopes"] <- "Dunes"
# landscape.raw$land_facet[landscape.raw$dunes_01 == 1 & landscape.raw$land_facet == "Sloped_mesic"] <- "Dunes"
# landscape.raw$land_facet[landscape.raw$dunes_01 == 1 & landscape.raw$land_facet == "Dry"] <- "Dunes"
# landscape.raw$land_facet[landscape.raw$dunes_01 == 1 & landscape.raw$land_facet == "Mesic"] <- "Dunes"

# Store the list of factor from the FEATURE_TY, Soil_Type_1, and land_facet columns
features.lookup <- list(as.factor(sort(unique(landscape.raw$land_facet))), 
                       droplevels(sort(unique(landscape.raw$soil_class))),
                       droplevels(sort(unique(landscape.raw$feature_ty))))
names(features.lookup) <- c("land_facet", "soil_class", "feature_ty")
save(features.lookup, file = "data/lookup/features-lookup_2021-10-07.Rdata")

# Create and save long forms for the four combinations (land facet quadrant, land facet site, soil quadrant, soil site)

landfacet.long.form <- list(aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "quadrant"], 
                                                landscape.raw[, "site_year"], 
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"], 
                                                landscape.raw[, "land_facet"]), 
                                      FUN = sum),
                            aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "site_year"], 
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"], 
                                                landscape.raw[, "land_facet"]), 
                                      FUN = sum),
                            aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "quadrant"], 
                                                landscape.raw[, "site_year"], 
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"], 
                                                landscape.raw[, "soil_class"]), 
                                      FUN = sum),
                            aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "site_year"], 
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"], 
                                                landscape.raw[, "soil_class"]), 
                                      FUN = sum))

names(landfacet.long.form) <- c("landfacet.quadrant", "landfacet.site", "soil.quadrant", "soil.site")
names(landfacet.long.form$landfacet.quadrant) <- c("quadrant", "site_year", "nr", "nsr", "feature_ty", "land_facet", "area")
names(landfacet.long.form$landfacet.site) <- c("site_year", "nr", "nsr", "feature_ty", "land_facet", "area")
names(landfacet.long.form$soil.quadrant) <- c("quadrant", "site_year", "nr", "nsr", "feature_ty", "soil_class", "area")
names(landfacet.long.form$soil.site) <- c("site_year", "nr", "nsr", "feature_ty", "soil_class", "area")

save(landfacet.long.form, file = "data/processed/facet-soil-longform_2021-10-07.Rdata")
rm(landfacet.long.form, landscape.raw)
gc()

####################
# Site proportions # 
####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load the landcover layer 
load("data/processed/facet-soil-longform_2021-10-07.Rdata")

# Create blank matrix to store the resutls in
# Should have a two lists (facet and soil) with the current and reference summaries at the site and quadrant level

# Facet Quadrant
facet.quad <- matrix_creation(site.list = unique(paste(landfacet.long.form$landfacet.quadrant$site_year, 
                                                 landfacet.long.form$landfacet.quadrant$quadrant, sep = "_")), 
                        feature.lookup = features.lookup$feature_ty[!(features.lookup$feature_ty %in% "NATIVE")], 
                        landscape.lookup = features.lookup$land_facet)
names(facet.quad) <- c("curr", "ref")

# Facet Site
facet.site <- matrix_creation(site.list = unique(paste(landfacet.long.form$landfacet.site$site_year)), 
                              feature.lookup = features.lookup$feature_ty[!(features.lookup$feature_ty %in% "NATIVE")], 
                              landscape.lookup = features.lookup$land_facet)
names(facet.site) <- c("curr", "ref")

# Soil Quadrant
soil.quad <- matrix_creation(site.list = unique(paste(landfacet.long.form$soil.quadrant$site_year, 
                                                       landfacet.long.form$soil.quadrant$quadrant, sep = "_")), 
                              feature.lookup = features.lookup$feature_ty[!(features.lookup$feature_ty %in% "NATIVE")], 
                              landscape.lookup = features.lookup$soil_class)
names(soil.quad) <- c("curr", "ref")

# Soil Site
soil.site <- matrix_creation(site.list = unique(paste(landfacet.long.form$soil.site$site_year)), 
                              feature.lookup = features.lookup$feature_ty[!(features.lookup$feature_ty %in% "NATIVE")], 
                              landscape.lookup = features.lookup$soil_class)
names(soil.site) <- c("curr", "ref")

# Create proportion summaries
facet.quad <- proportion_summary(facet.quad, landfacet.long.form$landfacet.quadrant, soil.facet.hf.lookup, TRUE, TRUE)
facet.site <- proportion_summary(facet.site, landfacet.long.form$landfacet.site, soil.facet.hf.lookup, TRUE, TRUE)
soil.quad <- proportion_summary(soil.quad, landfacet.long.form$soil.quadrant, soil.facet.hf.lookup, FALSE, TRUE)
soil.site <- proportion_summary(soil.site, landfacet.long.form$soil.site, soil.facet.hf.lookup, FALSE, TRUE)

landscape.summaries <- list(facet.quad, facet.site, soil.quad, soil.site)
names(landscape.summaries) <- c("landfacet.quadrant", "landfacet.site", "soil.quadrant", "soil.site")

save(landscape.summaries, file = "data/processed/facet-soil-proportions_2021-10-07.Rdata")

rm(facet.quad, facet.site, soil.quad, soil.site, landfacet.long.form, landscape.summaries)

#####################
# Kgrid proportions # 
#####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

facet.grid <- read.csv("data/base/habitat/land-facet-kgrid.csv")
colnames(facet.grid) <- c("LinkID", "Water_only", "Springs", "SnowIce", 
                          "Wet", "Mesic", "Dry", "Sloped_wet", 
                          "Sloped_mesic", "Cool_slopes", "Warm_slopes")
facet.grid["Water"] <- rowSums(facet.grid[, 2:4]) # Combine water categories
facet.grid <- facet.grid[, c(1, 5:12)]
facet.grid[, 2:9] <- facet.grid[, 2:9] / rowSums(facet.grid[, 2:9])

save(facet.grid, file = "data/processed/facet-kgrid_2021-10-07.Rdata")

load("data/base/habitat/veg-hf_1kmgrid_v6-fixage0.Rdata")
soil.ref <- as.data.frame(as.matrix(dd1km_pred$soil_reference))

rm(dd1km_pred)

# Soil grid
soil.grid <- as.data.frame(landscape_hf_reclass(soil.ref, soil.facet.hf.lookup))
soil.grid <- soil.grid / rowSums(soil.grid)
soil.grid["LinkID"] <- rownames(soil.grid)
soil.grid <- soil.grid[, c(7, 1:6)]

save(soil.grid, file = "data/processed/soil-kgrid_2021-10-07.Rdata")

rm(list=ls())
gc()
