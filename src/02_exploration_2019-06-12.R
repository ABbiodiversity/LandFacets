#
# Title: Visual exploration of land facet + soil data
# Created: June 12th, 2019
# Last Updated: August 19th, 2019
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

load("data/processed/facet-soil-proportions_2019-08-16.Rdata")
load("data/processed/soil-raw_2019-08-16.Rdata")

facet.temp <- as.data.frame(landscape.summaries$landfacet.quadrant$curr)
facet.temp["ID"] <- rownames(facet.temp)
soil.temp <- as.data.frame(soil.quad.raw$curr)
soil.temp["ID"] <- rownames(soil.temp)

site.curr <- merge.data.frame(soil.temp, facet.temp, "ID")

facet.temp <- as.data.frame(landscape.summaries$landfacet.quadrant$ref)
facet.temp["ID"] <- rownames(facet.temp)
soil.temp <- as.data.frame(soil.quad.raw$ref)
soil.temp["ID"] <- rownames(soil.temp)

site.ref <- merge.data.frame(soil.temp, facet.temp, "ID")

rm(soil.temp, facet.temp, landscape.summaries, soil.quad.raw)

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
png(filename = "results/figures/facet-soil-hf-curr_corr_2019-08-19.png", width = 800, height = 800)
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
png(filename = "results/figures/facet-soil-hf-ref_corr_2019-08-19.png", width = 800, height = 800)
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

# Number of sites for each landscape class
variable.occurrence <- apply(site.curr[, -1], MARGIN = 2, FUN = function(x) sum(ifelse(x > 0, 1, 0)))

png(filename = "results/figures/facet-soil-hf_distributions_2019-08-19.png",
    height = 2400,
    width = 2400)

par(mfrow = c(5, 4), cex = 2)

for (x in 1:length(variable.occurrence)) {
        
        plot(sort(site.curr[, (x + 1)][site.curr[, (x + 1)] > 0]), 
             ylim = c(0,1), 
             xlab = "Index",
             ylab = "% Area",
             main = paste(names(variable.occurrence)[x], variable.occurrence[x], "sites", sep = " "))
        
}

dev.off()




