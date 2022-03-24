#
# Title: Site summaries for both soil and land facet layers
# Created: June 12th, 2019
# Last Updated: March 24th, 2022
# Author: Brandon Allen
# Objectives: Create site and quadrant level summaries of the soil and land facet layers
# Keywords: Initialization, Exploring layer, Long form, Site proportions, Kgrid proportions, Occurrences
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

# Load multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        library(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }
        
        if (numPlots==1) {
                print(plots[[1]])
                
        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                        
                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}

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
                            dunes_01 = landscape.raw$dunes_01,
                            land_facet = landscape.raw$land_facet)

levels(landscape.raw$feature_ty) <- c(levels(landscape.raw$feature_ty), "NATIVE")
landscape.raw$feature_ty[is.na(landscape.raw$feature_ty)] <- as.factor("NATIVE")

# Load and clean climate and spatial data
load("data/base/habitat/veg-hf_SITE-all-years-combined.RData")
clim_qha <- clim_qha[clim_qha$NRNAME %in% c("Grassland", "Parkland") | clim_qha$NSRNAME == "Dry Mixedwood", ]


# climate.raw <- data.frame(abmi_site = clim_qha$site_id,
#                           year = clim_qha$survey_year,
#                           site_year = xx$Site_YEAR,
#                           on_off = xx$On_Off,
#                           AHM = xx$AHM,
#                           Eref = xx$Eref,
#                           FFP = xx$FFP,
#                           MAP = xx$MAP,
#                           MAT = xx$MAT,
#                           MCMT = xx$MCMT,
#                           MWMT = xx$MWMT,
#                           PET = xx$PET,
#                           paspen = xx$pAspen,
#                           Lat = xx$PUBLIC_LATTITUDE,
#                           Long = xx$PUBLIC_LONGITUDE,
#                           Elevation = xx$ELEVATION,
#                           Lat2 = xx$PUBLIC_LATTITUDE * xx$PUBLIC_LATTITUDE,
#                           Long2 = xx$PUBLIC_LONGITUDE * xx$PUBLIC_LONGITUDE,
#                           LatLong = xx$PUBLIC_LATTITUDE * xx$PUBLIC_LONGITUDE,
#                           MWMT2 = xx$MWMT * xx$MWMT,
#                           MAT2 = xx$MAT * xx$MAT,
#                           MAPPET = xx$MAP * xx$PET,
#                           MAPFFP = xx$MAP * xx$FFP,
#                           MATAHM = xx$MAT * xx$AHM)

climate.raw <- data.frame(abmi_site = clim_qha$site_id,
                          year = clim_qha$survey_year,
                          site_year = paste(clim_qha$site_id, clim_qha$survey_year, sep = "_"),
                          on_off = clim_qha$offgrid,
                          AHM = clim_qha$AHM,
                          FFP = clim_qha$FFP,
                          MAP = clim_qha$MAP,
                          MAT = clim_qha$MAT,
                          MCMT = clim_qha$MCMT,
                          MWMT = clim_qha$MWMT,
                          PET = clim_qha$PET,
                          paspen = clim_qha$pAspen,
                          Lat = clim_qha$X,
                          Long = clim_qha$Y,
                          Lat2 = clim_qha$X * clim_qha$X,
                          Long2 = clim_qha$Y * clim_qha$Y,
                          LatLong = clim_qha$X * clim_qha$Y,
                          MWMT2 = clim_qha$MWMT * clim_qha$MWMT,
                          MAT2 = clim_qha$MAT * clim_qha$MAT,
                          MAPPET = clim_qha$MAP * clim_qha$PET,
                          MAPFFP = clim_qha$MAP * clim_qha$FFP,
                          MATAHM = clim_qha$MAT * clim_qha$AHM)

rm(dd_1ha, dd_564m, dd_qha, clim_1ha, clim_qha)

# Filter the climate and landscape data to include only matching sites ( there are 40 sites of 968 in the climate data not found in landcover)
landscape.raw <- landscape.raw[landscape.raw$site_year %in% climate.raw$site_year, ]
climate.raw <- climate.raw[climate.raw$site_year %in% unique(landscape.raw$site_year), ]

# Aggregate the continuous topographic variables and bind to the climate data
name.store <- colnames(topo.continuous)
topo.continuous <- aggregate(topo.continuous[, 2:6], by = list(topo.continuous[, "UID"]), FUN = mean)
colnames(topo.continuous) <- c("site_year", name.store[-1])
topo.continuous$site_year <- as.factor(topo.continuous$site_year)
climate.raw <- merge.data.frame(climate.raw, topo.continuous, by = "site_year")

save(climate.raw, file = "data/processed/landcover/spatial-climate_2021-11-08.Rdata")

rm(climate.raw, topo.continuous, name.store)

###################
# Exploring layer #
###################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Land facet quadrants
facet.summary <- data.frame(Facet = names(table(landscape.raw$land_facet)),
                            Count = as.numeric(table(landscape.raw$land_facet)))
facet.summary$Facet <- c("Water", "Flat_wet", "Flat_mesic", "Flat_dry", "Sloped_wet",
                         "Sloped_mesic", "Cool_dry", "Warm_dry", "Dunes", "Saline")

# Each quadrant is 2500m
facet.area <- aggregate(landscape.raw$shape_area, 
                        by = list(Facet = landscape.raw$land_facet),
                        FUN = sum)
facet.area$Facet <- c("Water", "Flat_wet", "Flat_mesic", "Flat_dry", "Sloped_wet",
                       "Sloped_mesic", "Cool_dry", "Warm_dry", "Dunes", "Saline")
colnames(facet.area)[2] <- "Area"
facet.area$Area <- facet.area$Area / 2500
facet.summary <- merge.data.frame(facet.summary, facet.area, by = "Facet")

png(filename = "results/figures/land-facet-site-distribution_2021-10-15.png",
    height = 2400,
    width = 2400,
    res = 300)

ggplot(data = facet.summary) +
        geom_bar(stat="identity", aes(x = Facet, y = Area, color = "Area", fill = "Area"), show.legend = FALSE) +
        scale_color_manual(values = c("#2D415B")) +
        scale_fill_manual(values = c("#2D415B")) +
        xlab("Land Facet") +
        ylab("Equivalent Quadrants") + 
        theme_light()

dev.off()

rm(facet.summary, facet.area)

# Solar radiation by soil quadrants
facet.summary <- aggregate(landscape.raw$shape_area, by = list(Soil = landscape.raw$soil_class,
                                                               Solar = landscape.raw$solar_123),
                           FUN = sum)

colnames(facet.summary)[3] <- "Area"
facet.summary$Area <- facet.summary$Area / 2500

soil.facet.hf.lookup <- read.csv("data/lookup/soil-facet-hf-lookup-v61_2021.csv")

soil.category <- c("Loamy", "SandyLoam", 
                   "ClaySub", "RapidDrain", "Blowout", 
                   "ThinBreak", "Other")

simplified.soil <- NULL

for (id in soil.category) {
        
        soil.subset <- soil.facet.hf.lookup[soil.facet.hf.lookup$Class_Use %in% id, "Class_Fine"]
        facet.area <- facet.summary[facet.summary$Soil %in% soil.subset, ]
        facet.area$Soil <- id
        facet.area <- aggregate(facet.area$Area, by = list(Soil = facet.area$Soil,
                                                                       Solar = facet.area$Solar),
                                   FUN = sum)
        colnames(facet.area)[3] <- "Area"
        
        simplified.soil <- rbind(simplified.soil, facet.area)
        
}

simplified.soil$Solar <- rep(c("Cool", "Normal", "Warm"), 7)

simplified.soil$Solar <- factor(simplified.soil$Solar)

png(filename = "results/figures/land-facet-soil-solar-site-distribution_2021-10-15.png",
    height = 2400,
    width = 2400,
    res = 300)

ggplot(data = simplified.soil) +
        geom_bar(position="dodge", stat="identity", aes(x = Soil, y = Area, fill = Solar)) +
        scale_color_manual(values = c("#62929a", "#ffbb44", "#ee8577")) +
        scale_fill_manual(values = c("#62929a", "#ffbb44", "#ee8577")) +
        xlab("Soil") +
        ylab("Equivalent Quadrants") + 
        theme_light()

dev.off()

rm(facet.summary, facet.area)

# Wetness by soil quadrants
facet.summary <- aggregate(landscape.raw$shape_area, by = list(Soil = landscape.raw$soil_class,
                                                               Wetness = landscape.raw$cti_123),
                           FUN = sum)

colnames(facet.summary)[3] <- "Area"
facet.summary$Area <- facet.summary$Area / 2500

soil.facet.hf.lookup <- read.csv("data/lookup/soil-facet-hf-lookup-v61_2021.csv")

soil.category <- c("Loamy", "SandyLoam", 
                   "ClaySub", "RapidDrain", "Blowout", 
                   "ThinBreak", "Other")

simplified.soil <- NULL

for (id in soil.category) {
        
        soil.subset <- soil.facet.hf.lookup[soil.facet.hf.lookup$Class_Use %in% id, "Class_Fine"]
        facet.area <- facet.summary[facet.summary$Soil %in% soil.subset, ]
        facet.area$Soil <- id
        facet.area <- aggregate(facet.area$Area, by = list(Soil = facet.area$Soil,
                                                           Wetness = facet.area$Wetness),
                                FUN = sum)
        colnames(facet.area)[3] <- "Area"
        
        simplified.soil <- rbind(simplified.soil, facet.area)
        
}

simplified.soil$Wetness <- rep(c("Dry", "Mesic", "Wet"), 7)

simplified.soil$Wetness <- factor(simplified.soil$Wetness)

png(filename = "results/figures/land-facet-soil-wetness-site-distribution_2021-10-15.png",
    height = 2400,
    width = 2400,
    res = 300)

ggplot(data = simplified.soil) +
        geom_bar(position="dodge", stat="identity", aes(x = Soil, y = Area, fill = Wetness)) +
        scale_color_manual(values = c("#ee8577", "#ffbb44", "#62929a")) +
        scale_fill_manual(values = c("#ee8577", "#ffbb44", "#62929a")) +
        xlab("Soil") +
        ylab("Equivalent Quadrants") + 
        theme_light()

dev.off()

# Create Terrain layer
landscape.raw["terrain"] <- paste(landscape.raw$cti_123, landscape.raw$solar_123, sep = "_")
landscape.raw$terrain[landscape.raw$terrain == "1_1"] <- "Cool_slopes"
landscape.raw$terrain[landscape.raw$terrain == "1_2"] <- "Dry"
landscape.raw$terrain[landscape.raw$terrain == "1_3"] <- "Warm_slopes"
landscape.raw$terrain[landscape.raw$terrain == "2_1"] <- "Sloped_mesic"
landscape.raw$terrain[landscape.raw$terrain == "2_2"] <- "Mesic"
landscape.raw$terrain[landscape.raw$terrain == "2_3"] <- "Sloped_mesic"
landscape.raw$terrain[landscape.raw$terrain == "3_1"] <- "Sloped_wet"
landscape.raw$terrain[landscape.raw$terrain == "3_2"] <- "Wet"
landscape.raw$terrain[landscape.raw$terrain == "3_3"] <- "Sloped_wet"

# Rule set for creating the terrain x Water land facet layer 
# terrain x water_01 x springs_01 x snowice_01 = terrain
# If any of the three water groups = 1, replace value from the terrain model

# Create terrain + water layer
landscape.raw["terrain_water"] <- landscape.raw$terrain
landscape.raw$terrain_water[landscape.raw$water_01 == 1] <- "Water"
landscape.raw$terrain_water[landscape.raw$snowice_01 == 1] <- "Snowice"
landscape.raw$terrain_water[landscape.raw$springs_01 == 1] <- "Springs"

# Wetness by soil quadrants
facet.summary <- aggregate(landscape.raw$shape_area, by = list(Terrain = landscape.raw$terrain_water),
                           FUN = sum)

colnames(facet.summary)[2] <- "Area"
facet.summary$Area <- facet.summary$Area / 2500

png(filename = "results/figures/terrain-water-site-distribution_2021-10-15.png",
    height = 2400,
    width = 2400,
    res = 300)

ggplot(data = facet.summary) +
    geom_bar(stat="identity", aes(x = Terrain, y = Area, color = "Area", fill = "Area"), show.legend = FALSE) +
    scale_color_manual(values = c("#2D415B")) +
    scale_fill_manual(values = c("#2D415B")) +
    xlab("Terrain by Water") +
    ylab("Equivalent Quadrants") + 
    theme_light()

dev.off()

# # Flatness by soil quadrants
# landscape.raw["Flatness"] <- paste(landscape.raw$cti_123, landscape.raw$solar_123, sep = "_")
# landscape.raw$Flatness[landscape.raw$Flatness == "1_1"] <- "Sloped"
# landscape.raw$Flatness[landscape.raw$Flatness == "1_2"] <- "Flat"
# landscape.raw$Flatness[landscape.raw$Flatness == "1_3"] <- "Sloped"
# landscape.raw$Flatness[landscape.raw$Flatness == "2_1"] <- "Sloped"
# landscape.raw$Flatness[landscape.raw$Flatness == "2_2"] <- "Flat"
# landscape.raw$Flatness[landscape.raw$Flatness == "2_3"] <- "Sloped"
# landscape.raw$Flatness[landscape.raw$Flatness == "3_1"] <- "Sloped"
# landscape.raw$Flatness[landscape.raw$Flatness == "3_2"] <- "Flat"
# landscape.raw$Flatness[landscape.raw$Flatness == "3_3"] <- "Sloped"
# 
# facet.summary <- aggregate(landscape.raw$shape_area, by = list(Soil = landscape.raw$soil_class,
#                                                                Flatness = landscape.raw$Flatness),
#                            FUN = sum)
# 
# colnames(facet.summary)[3] <- "Area"
# facet.summary$Area <- facet.summary$Area / 2500
# 
# soil.facet.hf.lookup <- read.csv("data/lookup/soil-facet-hf-lookup-v61_2021.csv")
# 
# soil.category <- c("Loamy", "SandyLoam", 
#                    "ClaySub", "RapidDrain", "Blowout", 
#                    "ThinBreak", "Other")
# 
# simplified.soil <- NULL
# 
# for (id in soil.category) {
#     
#     soil.subset <- soil.facet.hf.lookup[soil.facet.hf.lookup$Class_Use %in% id, "Class_Fine"]
#     facet.area <- facet.summary[facet.summary$Soil %in% soil.subset, ]
#     facet.area$Soil <- id
#     facet.area <- aggregate(facet.area$Area, by = list(Soil = facet.area$Soil,
#                                                        Flatness = facet.area$Flatness),
#                             FUN = sum)
#     colnames(facet.area)[3] <- "Area"
#     
#     simplified.soil <- rbind(simplified.soil, facet.area)
#     
# }
# 
# simplified.soil$Flatness <- factor(simplified.soil$Flatness)
# 
# png(filename = "results/figures/land-facet-soil-flatness-site-distribution_2021-10-15.png",
#     height = 2400,
#     width = 2400,
#     res = 300)
# 
# ggplot(data = simplified.soil) +
#     geom_bar(position="dodge", stat="identity", aes(x = Soil, y = Area, fill = Flatness)) +
#     xlab("Soil") +
#     ylab("Equivalent Quadrants") + 
#     theme_light()
# 
# dev.off()

# # Kgrid distribution
# facet.grid <- read.csv("data/base/habitat/land-facet-kgrid.csv")
# colnames(facet.grid) <- c("LinkID", "Water_only", "Springs", "SnowIce", 
#                           "Wet", "Mesic", "Dry", "Sloped_wet", 
#                           "Sloped_mesic", "Cool_slopes", "Warm_slopes")
# facet.grid["Water"] <- rowSums(facet.grid[, 2:4]) # Combine water categories
# facet.grid <- facet.grid[, c(1, 5:12)]
# facet.grid[, 2:9] <- facet.grid[, 2:9] / rowSums(facet.grid[, 2:9]) # Convert to proportions (~km2)
# 
# load("data/base/habitat/kgrid_table_km.Rdata")
# kgrid <- kgrid[, c("Row_Col", "NSRNAME", "NRNAME")]
# colnames(kgrid)[1] <- "LinkID"
# facet.grid <- merge.data.frame(kgrid, facet.grid, by = "LinkID")
# rm(kgrid)
# 
# facet.grid <- facet.grid[facet.grid$NRNAME %in% c("Grassland", "Parkland") | facet.grid$NSRNAME == "Dry Mixedwood", ]
# 
# grid.summary <- colSums(facet.grid[, 4:11])
# grid.summary <- data.frame(Facet = names(grid.summary),
#                            Area = as.numeric(grid.summary))
# 
# png(filename = "results/figures/land-facet-kgrid-distribution_2021-10-15.png",
#     height = 4800,
#     width = 2400,
#     res = 300)
# 
# ggplot(data = grid.summary) +
#         geom_bar(stat="identity", aes(x = Facet, y = Area, color = "Area", fill = "Area"), show.legend = FALSE) +
#         scale_color_manual(values = c("#2D415B")) +
#         scale_fill_manual(values = c("#2D415B")) +
#         xlab("Land Facet") +
#         ylab("Area (km2)") + 
#         theme_light()
# 
# dev.off()

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
landscape.raw["terrain"] <- paste(landscape.raw$cti_123, landscape.raw$solar_123, sep = "_")
landscape.raw$terrain[landscape.raw$terrain == "1_1"] <- "Cool_slopes"
landscape.raw$terrain[landscape.raw$terrain == "1_2"] <- "Dry"
landscape.raw$terrain[landscape.raw$terrain == "1_3"] <- "Warm_slopes"
landscape.raw$terrain[landscape.raw$terrain == "2_1"] <- "Sloped_mesic"
landscape.raw$terrain[landscape.raw$terrain == "2_2"] <- "Mesic"
landscape.raw$terrain[landscape.raw$terrain == "2_3"] <- "Sloped_mesic"
landscape.raw$terrain[landscape.raw$terrain == "3_1"] <- "Sloped_wet"
landscape.raw$terrain[landscape.raw$terrain == "3_2"] <- "Wet"
landscape.raw$terrain[landscape.raw$terrain == "3_3"] <- "Sloped_wet"

# Rule set for creating the terrain x Water land facet layer 
# terrain x water_01 x springs_01 x snowice_01 = terrain
# If any of the three water groups = 1, replace value from the terrain model

# Create terrain + water layer
landscape.raw["terrain_water"] <- landscape.raw$terrain
landscape.raw$terrain_water[landscape.raw$water_01 == 1] <- "Water"
landscape.raw$terrain_water[landscape.raw$snowice_01 == 1] <- "Snowice"
landscape.raw$terrain_water[landscape.raw$springs_01 == 1] <- "Springs"

# Reclass land facet into names
landscape.raw$land_facet[landscape.raw$land_facet == 0] <- "Water"
landscape.raw$land_facet[landscape.raw$land_facet == 3] <- "Wet"
landscape.raw$land_facet[landscape.raw$land_facet == 4] <- "Mesic"
landscape.raw$land_facet[landscape.raw$land_facet == 5] <- "Dry"
landscape.raw$land_facet[landscape.raw$land_facet == 6] <- "Sloped_wet"
landscape.raw$land_facet[landscape.raw$land_facet == 7] <- "Sloped_mesic"
landscape.raw$land_facet[landscape.raw$land_facet == 8] <- "Cool_slopes"
landscape.raw$land_facet[landscape.raw$land_facet == 9] <- "Warm_slopes"
landscape.raw$land_facet[landscape.raw$land_facet == 14] <- "Dunes"
landscape.raw$land_facet[landscape.raw$land_facet == 15] <- "Saline"

# Factor the land facet and soil categories
landscape.raw$land_facet <- factor(landscape.raw$land_facet)
landscape.raw$terrain_water <- factor(landscape.raw$terrain_water)
landscape.raw$soil_class <- factor(landscape.raw$soil_class)
landscape.raw$soil_wet <- factor(paste0(landscape.raw$soil_class, "_", landscape.raw$cti_123))


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
                        as.factor(sort(unique(landscape.raw$terrain_water))), 
                        droplevels(sort(unique(landscape.raw$soil_class))),
                        droplevels(sort(unique(landscape.raw$feature_ty))),
                        factor(sort(unique(landscape.raw$soil_wet))))
names(features.lookup) <- c("land_facet", "terrain_water", "soil_class", "feature_ty", "soil_wet")
save(features.lookup, file = "data/lookup/features-lookup_2021-10-07.Rdata")

# Add the site and year information for continuity
load("data/processed/landcover/spatial-climate_2021-11-08.Rdata")
climate.raw <- climate.raw[, c("abmi_site", "year", "site_year")]
landscape.raw <- merge.data.frame(landscape.raw, climate.raw, by = "site_year")

# Create and save long forms for the four combinations (land facet quadrant, land facet site, soil quadrant, soil site)

landfacet.long.form <- list(aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "abmi_site"], 
                                                landscape.raw[, "year"], 
                                                landscape.raw[, "site_year"], 
                                                landscape.raw[, "quadrant"],
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"], 
                                                landscape.raw[, "land_facet"]), 
                                      FUN = sum),
                            aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "abmi_site"], 
                                                landscape.raw[, "year"], 
                                                landscape.raw[, "site_year"],
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"], 
                                                landscape.raw[, "land_facet"]), 
                                      FUN = sum),
                            aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "abmi_site"], 
                                                landscape.raw[, "year"], 
                                                landscape.raw[, "site_year"], 
                                                landscape.raw[, "quadrant"], 
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"], 
                                                landscape.raw[, "terrain_water"]), 
                                      FUN = sum),
                            aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "abmi_site"], 
                                                landscape.raw[, "year"], 
                                                landscape.raw[, "site_year"], 
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"],
                                                landscape.raw[, "terrain_water"]), 
                                      FUN = sum),
                            
                            aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "abmi_site"], 
                                                landscape.raw[, "year"], 
                                                landscape.raw[, "site_year"], 
                                                landscape.raw[, "quadrant"], 
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"], 
                                                landscape.raw[, "soil_class"]), 
                                      FUN = sum),
                            aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "abmi_site"], 
                                                landscape.raw[, "year"], 
                                                landscape.raw[, "site_year"],  
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"],
                                                landscape.raw[, "soil_class"]), 
                                      FUN = sum),
                            
                            aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "abmi_site"], 
                                                landscape.raw[, "year"], 
                                                landscape.raw[, "site_year"], 
                                                landscape.raw[, "quadrant"],
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"], 
                                                landscape.raw[, "soil_wet"]), 
                                      FUN = sum),
                            aggregate(landscape.raw[, "shape_area"], 
                                      by = list(landscape.raw[, "abmi_site"], 
                                                landscape.raw[, "year"], 
                                                landscape.raw[, "site_year"],  
                                                landscape.raw[, "nr"], 
                                                landscape.raw[, "nsr"],
                                                landscape.raw[, "feature_ty"],
                                                landscape.raw[, "soil_wet"]), 
                                      FUN = sum))

names(landfacet.long.form) <- c("landfacet.quadrant", "landfacet.site", "terrain.quadrant", "terrain.site", 
                                "soil.quadrant", "soil.site", "soil.wet.quadrant", "soil.wet.site")
names(landfacet.long.form$landfacet.quadrant) <- c("site", "year", "site_year", "quadrant", "nr", "nsr", "feature_ty", "land_facet", "area")
names(landfacet.long.form$landfacet.site) <- c("site", "year", "site_year", "nr", "nsr", "feature_ty", "land_facet", "area")
names(landfacet.long.form$terrain.quadrant) <- c("site", "year", "site_year", "quadrant", "nr", "nsr", "feature_ty", "terrain_water", "area")
names(landfacet.long.form$terrain.site) <- c("site", "year", "site_year", "nr", "nsr", "feature_ty", "terrain_water", "area")
names(landfacet.long.form$soil.quadrant) <- c("site", "year", "site_year", "quadrant", "nr", "nsr", "feature_ty", "soil_class", "area")
names(landfacet.long.form$soil.site) <- c("site", "year", "site_year", "nr", "nsr", "feature_ty", "soil_class", "area")
names(landfacet.long.form$soil.wet.quadrant) <- c("site", "year", "site_year", "quadrant", "nr", "nsr", "feature_ty", "soil_wet", "area")
names(landfacet.long.form$soil.wet.site) <- c("site", "year", "site_year", "nr", "nsr", "feature_ty", "soil_wet", "area")

save(landfacet.long.form, file = "data/processed/landcover/facet-soil-longform_2021-11-08.Rdata")
rm(landfacet.long.form, landscape.raw)
gc()

####################
# Site proportions # 
####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load the landcover layer 
load("data/processed/landcover/facet-soil-longform_2021-11-08.Rdata")

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

# Terrain Quadrant
terrain.quad <- matrix_creation(site.list = unique(paste(landfacet.long.form$terrain.quadrant$site_year, 
                                                       landfacet.long.form$terrain.quadrant$quadrant, sep = "_")), 
                              feature.lookup = features.lookup$feature_ty[!(features.lookup$feature_ty %in% "NATIVE")], 
                              landscape.lookup = features.lookup$terrain_water)
names(terrain.quad) <- c("curr", "ref")

# Terrain Site
terrain.site <- matrix_creation(site.list = unique(paste(landfacet.long.form$terrain.site$site_year)), 
                              feature.lookup = features.lookup$feature_ty[!(features.lookup$feature_ty %in% "NATIVE")], 
                              landscape.lookup = features.lookup$terrain_water)
names(terrain.site) <- c("curr", "ref")

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

# Soil Wetness Quadrant
soil.wet.quad <- matrix_creation(site.list = unique(paste(landfacet.long.form$soil.wet.quadrant$site_year, 
                                                      landfacet.long.form$soil.wet.quadrant$quadrant, sep = "_")), 
                             feature.lookup = features.lookup$feature_ty[!(features.lookup$feature_ty %in% "NATIVE")], 
                             landscape.lookup = features.lookup$soil_wet)
names(soil.wet.quad) <- c("curr", "ref")

# Soil Wetness Site
soil.wet.site <- matrix_creation(site.list = unique(paste(landfacet.long.form$soil.wet.site$site_year)), 
                             feature.lookup = features.lookup$feature_ty[!(features.lookup$feature_ty %in% "NATIVE")], 
                             landscape.lookup = features.lookup$soil_wet)
names(soil.wet.site) <- c("curr", "ref")

# Load relevant lookup tables
soil.facet.hf.lookup <- read.csv("data/lookup/soil-facet-hf-lookup-v61_2021.csv")
soil.moisture.hf.lookup <- read.csv("data/lookup/soil-moisture-hf-lookup-v61_2021.csv")

# Create proportion summaries
facet.quad <- proportion_summary(facet.quad, landfacet.long.form$landfacet.quadrant, soil.facet.hf.lookup, landscape.ty = "land_facet", TRUE)
facet.site <- proportion_summary(facet.site, landfacet.long.form$landfacet.site, soil.facet.hf.lookup, landscape.ty = "land_facet", TRUE)
terrain.quad <- proportion_summary(terrain.quad, landfacet.long.form$terrain.quadrant, soil.facet.hf.lookup, landscape.ty = "terrain_water", TRUE)
terrain.site <- proportion_summary(terrain.site, landfacet.long.form$terrain.site, soil.facet.hf.lookup, landscape.ty = "terrain_water", TRUE)
soil.quad <- proportion_summary(soil.quad, landfacet.long.form$soil.quadrant, soil.facet.hf.lookup, landscape.ty = "soil_class", TRUE)
soil.site <- proportion_summary(soil.site, landfacet.long.form$soil.site, soil.facet.hf.lookup, landscape.ty = "soil_class", TRUE)
soil.wet.quad <- proportion_summary(soil.wet.quad, landfacet.long.form$soil.wet.quadrant, soil.moisture.hf.lookup, landscape.ty = "soil_wet", TRUE)
soil.wet.site <- proportion_summary(soil.wet.site, landfacet.long.form$soil.wet.site, soil.moisture.hf.lookup, landscape.ty = "soil_wet", TRUE)

landscape.summaries <- list(facet.quad, facet.site, terrain.quad, terrain.site, 
                            soil.quad, soil.site, soil.wet.quad, soil.wet.site)
names(landscape.summaries) <- c("landfacet.quadrant", "landfacet.site", "terrain.quadrant", "terrain.site",
                                "soil.quadrant", "soil.site", "soil.wet.quadrant", "soil.wet.site")

save(landscape.summaries, file = "data/processed/landcover/facet-soil-proportions_2021-11-08.Rdata")

rm(list=ls())
gc()

#####################
# Kgrid proportions # 
#####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# facet.grid <- read.csv("data/base/habitat/land-facet-kgrid.csv")
# colnames(facet.grid) <- c("LinkID", "Water_only", "Springs", "SnowIce", 
#                           "Wet", "Mesic", "Dry", "Sloped_wet", 
#                           "Sloped_mesic", "Cool_slopes", "Warm_slopes")
# facet.grid["Water"] <- rowSums(facet.grid[, 2:4]) # Combine water categories
# facet.grid <- facet.grid[, c(1, 5:12)]
# facet.grid[, 2:9] <- facet.grid[, 2:9] / rowSums(facet.grid[, 2:9])
# 
# save(facet.grid, file = "data/processed/facet-kgrid_2021-10-07.Rdata")
# 
# load("data/base/habitat/veg-hf_1kmgrid_v6-fixage0.Rdata")
# soil.ref <- as.data.frame(as.matrix(dd1km_pred$soil_reference))
# 
# rm(dd1km_pred)
# 
# # Soil grid
# soil.grid <- as.data.frame(landscape_hf_reclass(soil.ref, soil.facet.hf.lookup))
# soil.grid <- soil.grid / rowSums(soil.grid)
# soil.grid["LinkID"] <- rownames(soil.grid)
# soil.grid <- soil.grid[, c(7, 1:6)]
# 
# save(soil.grid, file = "data/processed/soil-kgrid_2021-10-07.Rdata")
# 
# rm(list=ls())
# gc()

###############
# Occurrences # 
###############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#
# Site
#

# Load the landcover data
load("data/processed/landcover/facet-soil-proportions_2021-11-08.Rdata")

# Define the site list
site.list <- rownames(landscape.summaries$landfacet.site$curr)

# Define the species files
species.list <- list.files("data/base/species/", full.names = TRUE)
species.list <- species.list[grep(".RData", species.list)] # Grab only the R data ones
names(species.list) <- c("lichen", "mite", "moss", "vplant") # Name

# For each species data, standardize

for(taxon in names(species.list)) {
        
        # Load taxon
        load(species.list[taxon])
        
        # Create quandrant column
        d$nQuadrant <- 1
        
        # Aggregate by site_year
        occurrence.in <- aggregate(d[, c(FirstSpCol:(LastSpCol), ncol(d))], by = list(d$SiteYear), FUN = sum)
        colnames(occurrence.in)[1] <- "site_year"
        occurrence.in <- occurrence.in[, c(1, ncol(occurrence.in), 2:(ncol(occurrence.in) - 1))]
        
        # Filter to match sites in the landscape summary
        occurrence.in <- occurrence.in[occurrence.in$site_year %in% site.list, ]
        
        # Save results
        save(occurrence.in, file = paste0("data/processed/occurrence/", taxon, "-site-occurrence.RData"))
        
        # Remove old information
        rm(d, pm, FirstSpCol, LastSpCol, SpTable, SpTable.ua, occurrence.in)
        
}

rm(list=ls())
gc()

#
# Quadrant
#

# Load the landcover data
load("data/processed/landcover/facet-soil-proportions_2021-11-08.Rdata")

# Define the site list
site.list <- rownames(landscape.summaries$landfacet.quadrant$curr)

# Define the species files
species.list <- list.files("data/base/species/", full.names = TRUE)
species.list <- species.list[grep(".RData", species.list)] # Grab only the R data ones
names(species.list) <- c("lichen", "mite", "moss", "vplant") # Name

# For each species data, standardize

for(taxon in names(species.list)) {
    
    # Load taxon
    load(species.list[taxon])
    
    # Create quandrant column
    d$nQuadrant <- 1
    occurrence.in <- d[, c(1:5, ncol(d), FirstSpCol:(LastSpCol))]
    
    # Filter to match sites in the landscape summary
    occurrence.in <- occurrence.in[occurrence.in$SiteYearQu %in% site.list, ]
    
    # Save results
    save(occurrence.in, file = paste0("data/processed/occurrence/", taxon, "-quadrant-occurrence.RData"))
    
    # Remove old information
    rm(d, pm, FirstSpCol, LastSpCol, SpTable, SpTable.ua, occurrence.in)
    
}

rm(list=ls())
gc()
