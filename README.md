# Land facets versus soil composition: how are our models impacted?
Brandon Allen & Ermias Azeria

October 7th, 2021

## Introduction
The southern habitat models created by the ABMI use soil information aggregated from the AGRASID database. This database is composed of complex polygons which can contain multiple soil classes. As we don’t know how these soil classes are distributed within the polygon, we have chosen to reclassify the polygon by the dominate soil class. However, this results in lost information, and does not guarantee the soil summaries are accurrate at our sites. This affects both the quality of our models and our provincial predictions onto the backfilled soil layer. An alternative to this soil layer could be the recently created land facet layer developed by Scott Nielsen. This layer categorizes 15 m2 pixels into 16 categories by combining three sub-models. 

1) The Water sub-model is comprised of surface water (ABMI 2010 landcover), springs (Alberta Geological Survey) and snow and ice (ABMI 2010 landcover). 
2) The Terrain sub-model is based on four topographic indices (topographic wetness, solar radiation, aspect, and slope).
3) The Geology-soil features sub-model is comprised of exposed rock (ABMI 2010 landcover and Alberta Geological Survey), dunes (Alberta Geological Survey), and saline soils (AGRASID). 

With a higher spatial resolution, and larger number of classes compared to our soil layer (16 vs 4), there may be benefits to incorporating this GIS layer into our habitat suitability modeling framework.

As this layer combines several data sources, there are a few issues we should be aware of. First, the saline soils layer is derived from the AGRASID database, where polygons are defined as saline if they are greater than 30% saline soil. This means the soil composition at our sites may be incorrect as we don’t know the distribution of saline soil in the polygon. Second, the dunes and saline soils layers have pixels masked if the topographic layers indicate the pixel is wetter or drier than expected. As we don’t have ground truthing, it makes it difficult to know which classification is correct. As Scott created the layer in a modular framework, we could choose to ignore the saline soils and/or dunes layers or create a different ruleset for stacking layers that would better meet our needs. For this exercise, I chose to remove the saline soils layer and maintain the dunes layer based on Scott's stacking process. This processes involved reclassifying Mesic, Dry, and Sloped pixels from the Terrain sub-model that overlapped with the Alberta Geological Surveys dunes layer.

## Landcover Distribution
Landcover data was summarized for all ABMI sites surveyed between 2003-2018. At each site, both land facet and soil information were summarized at the 1ha and quadrant scales. Visual inspection of the data shows that under the current landscape condition (including human footprint), there are positive correlations between the land facet and soil classes (Figure 1; e.g., ClayLoamSand and Dry, Blowout and Mesic, ThinBreak and Warm slope). However, many of those relationships were not observed under the reference condition (Figure 2). 


![Heatmap of correlation between soil and land facets under current conditions](results/figures/facet-soil-hf-curr_corr_2021-10-07.png)
Figure 1 - Heatmap of correlation between soil and land facets under current conditions.

![Heatmap of correlation between soil and land facets under reference conditions](results/figures/facet-soil-hf-ref_corr_2021-10-07.png)
Figure 2 - eatmap of correlation between soil and land facets under reference conditions.

The proportion of sites occupied by landcover classes varied between soil and land facet categories. Based on the soil data, many sites were observed to be 100% one class. However, we observed few sites with 100% coverage of a single land facet (Figure 3). The increased complexity of each site could affect our models, but may represent the landscape more accurately.

![Distribution of landcover variables across ABMI sites (2003-2018)](results/figures/facet-soil-hf_distributions_2021-10-07.png)
Figure 3 - Distribution of landcover variables across ABMI sites (2003-2018).

## Model Results

THIS IS UNDERGOING A REVISION
I explored two different approaches for using the land facet layer in our vascular plant habitat suitability models. The first option uses the ABMI soils layer but adds the continuous topographic indices from the Terrain sub-model as another spatial and climate covariate. The second option replaces the ABMI soils layer with the categories defined in the modified land facet layer. Model fit of these approaches were compared to our current models built on only the ABMI soils layer. AUC scores were calculated using all available vascular plant data and cross validation was performed using a "leave-one-out" approach.  
