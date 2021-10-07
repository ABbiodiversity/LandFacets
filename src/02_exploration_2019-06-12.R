#
# Title: Visual exploration of land facet + soil data
# Created: June 12th, 2019
# Last Updated: October 7th, 2021
# Author: Brandon Allen
# Objectives: Visually explore the similarities between the land facet and soils layers
# Keywords: Initialization, Exploration
#

##################
# Initialization # 
##################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
gc()

# Load relevant libraries, data sets, format data for exploration
library(ggplot2)
library(reshape2)

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

load("data/processed/facet-soil-proportions_2021-10-07.Rdata")

# Remove the duplicated HF information from one of the current conditions 
facet.temp <- as.data.frame(landscape.summaries$landfacet.quadrant$curr)
facet.temp["ID"] <- rownames(facet.temp)
facet.temp <- facet.temp[, -c(9:16)]
soil.temp <- as.data.frame(landscape.summaries$soil.quadrant$curr)
soil.temp["ID"] <- rownames(soil.temp)

site.curr <- merge.data.frame(soil.temp, facet.temp, "ID")

facet.temp <- as.data.frame(landscape.summaries$landfacet.quadrant$ref)
facet.temp["ID"] <- rownames(facet.temp)
soil.temp <- as.data.frame(landscape.summaries$soil.quadrant$ref)
soil.temp["ID"] <- rownames(soil.temp)

site.ref <- merge.data.frame(soil.temp, facet.temp, "ID")

rm(soil.temp, facet.temp, landscape.summaries)

###############
# Exploration # 
###############~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Correlation for current landscape

cormat <- cor(site.curr[, -1])

# Get upper triangle of the correlation matrix

get_lower_tri<-function(cormat){
        cormat[upper.tri(cormat)] <- NA
        return(cormat)
}

upper_tri <- get_lower_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Plot heatmap
png(filename = "results/figures/facet-soil-hf-curr_corr_2021-10-07.png", width = 800, height = 800)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + # Visualization
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, limit = c(-1,1), space = "Lab", 
                             name="Pearson\nCorrelation") +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 12, hjust = 1))+
        coord_fixed()
dev.off()

# Correlation for reference landscape

cormat <- cor(site.ref[, -1])

# Get upper triangle of the correlation matrix

get_lower_tri<-function(cormat){
        cormat[upper.tri(cormat)] <- NA
        return(cormat)
}

upper_tri <- get_lower_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)
rownames(melted_cormat) <- seq(from = 1, to = nrow(melted_cormat), by = 1)

# Plot heatmap
png(filename = "results/figures/facet-soil-hf-ref_corr_2021-10-07.png", width = 800, height = 800)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + # Visualization
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, limit = c(-1,1), space = "Lab", 
                             name="Pearson\nCorrelation") +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 12, hjust = 1))+
        coord_fixed()
dev.off()

# Number of sites for each landscape class
variable.occurrence <- apply(site.curr[, -1], MARGIN = 2, FUN = function(x) sum(ifelse(x > 0, 1, 0)))

png(filename = "results/figures/facet-soil-hf_distributions_2021-10-07.png",
    height = 3600,
    width = 3600,
    res = 300)

# Loop through each plot\
plot.list <- list()
for (x in 1:length(variable.occurrence)) {
    
    temp.data <- data.frame(Site = 1:length(site.curr[, (x + 1)][site.curr[, (x + 1)] > 0]),
                            Area = sort(site.curr[, (x + 1)][site.curr[, (x + 1)] > 0]))
    
    plot.list[[names(variable.occurrence)[x]]] <- ggplot(data = temp.data) +
        geom_point(aes(x = Site, y = Area, color = "Area"), show.legend = FALSE) +
        scale_color_manual(values = c("#2D415B")) +
        xlab("Index") +
        ylab("% Area") + 
        ylim(c(0,1)) +
        ggtitle(paste(names(variable.occurrence)[x], variable.occurrence[x], "sites", sep = " ")) +
        theme_light()

    
    rm(temp.data)
    
}

multiplot(plot.list$Blowout, plot.list$RapidDrain, plot.list$ClayLoamSand, plot.list$Other, plot.list$ThinBreak, plot.list$Water, plot.list$Hwater,
      plot.list$Crop, plot.list$RoughP, plot.list$TameP, plot.list$UrbInd, plot.list$SoftLin, plot.list$HardLin, plot.list$HFor, 
      plot.list$Cool_slopes, plot.list$Dry, plot.list$Mesic, plot.list$Sloped_mesic, plot.list$Sloped_wet, plot.list$Warm_slopes, plot.list$Wet,
      cols = 3)

dev.off()

#
# Climate Correlations
# Prediction grid
pm <- read.csv("data/lookup/prediction-matrix/predictionmatrix_terrain-water-dunes.csv")
rownames(pm) <- colnames(pm)

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
model.data <- merge.data.frame(model.data, climate.raw, by = "site_year")

# Correlation for current landscape

cormat <- cor(model.data[, -c(1, 16:18, 33)])

# Get upper triangle of the correlation matrix

get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

upper_tri <- get_lower_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Plot heatmap
png(filename = "results/figures/facet-climate_corr_2019-09-19.png", width = 800, height = 800)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + # Visualization
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
dev.off()

