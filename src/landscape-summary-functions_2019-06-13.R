#
# Title: Landscape summary functions
# Created: June 13th, 2019
# Last Updated: June 13th, 2019
# Author: Brandon Allen
# Objectives: Functions required to summarize the landscape long-form data
# Keywords: Matrix Creation, Proportion Summaries, Soil + HF Reclass
#

###################
# Matrix Creation # 
###################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to create current and reference matrix based on site list and landscape lookups (soil, facet, HF)
matrix_creation <- function(site.list, feature.lookup, landscape.lookup) {
        
        matrix.curr <- matrix(nrow = length(unique(site.list)),
                              ncol = (length(landscape.lookup) + length(feature.lookup)))
        rownames(matrix.curr) <- unique(site.list)
        colnames(matrix.curr) <- c(as.character(landscape.lookup), as.character(feature.lookup))
        matrix.curr[is.na(matrix.curr)] <- 0 
        
        matrix.ref <- matrix(nrow = length(unique(site.list)),
                             ncol = (length(landscape.lookup)))
        rownames(matrix.ref) <- unique(site.list)
        colnames(matrix.ref) <- c(as.character(landscape.lookup))
        matrix.ref[is.na(matrix.ref)] <- 0 
        
        return(list(matrix.curr, matrix.ref))
        
}


########################
# Proportion Summaries # 
########################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to calculate the proportion matrix of each landscape class for the scale of interest.
# This is the dataset required for the down stream species modeling. 
proportion_summary <- function(blank.matrix, longform.matrix, feature.reclass, facet.soil, type.proportion) {
        
        landscape.ty <- ifelse(facet.soil == TRUE, "land_facet", "soil_class")
        
        temp.cur.native <- longform.matrix[longform.matrix$feature_ty %in% "NATIVE", ]
        temp.cur.HF  <- longform.matrix[!(longform.matrix$feature_ty %in% "NATIVE"), ]
        matrix.curr <- blank.matrix$curr
        matrix.ref <- blank.matrix$ref
        
        # Current HF
        if(nrow(temp.cur.HF) != 0) {
                
                if(table(colnames(temp.cur.HF) %in% c("quadrant", "site_year"))["TRUE"] == 2) {
                        
                        temp.cur.HF <- aggregate(temp.cur.HF[, "area"],
                                                 by = list(temp.cur.HF[, "quadrant"], 
                                                           temp.cur.HF[, "site_year"],
                                                           temp.cur.HF[, "feature_ty"]),
                                                 FUN = sum)
                        
                        colnames(temp.cur.HF) <- c("quadrant", "site_year", "feature_ty", "area")
                        temp.cur.HF["site_id"] <- paste(temp.cur.HF$site_year, temp.cur.HF$quadrant, sep = "_")
                        
                        for (y in 1:nrow(temp.cur.HF)) {
                                
                                matrix.curr[as.character(temp.cur.HF[y, "site_id"]), as.character(temp.cur.HF[y, "feature_ty"])] <- temp.cur.HF[y, "area"]
                                
                        }
                        
                } else {
                        
                        temp.cur.HF <- aggregate(temp.cur.HF[, "area"],
                                                 by = list(temp.cur.HF[, "site_year"],
                                                           temp.cur.HF[, "feature_ty"]),
                                                 FUN = sum)
                        
                        colnames(temp.cur.HF) <- c("site_year", "feature_ty", "area")
                        temp.cur.HF["site_id"] <- temp.cur.HF$site_year
                        
                        for (y in 1:nrow(temp.cur.HF)) {
                                
                                matrix.curr[as.character(temp.cur.HF[y, "site_id"]), as.character(temp.cur.HF[y, "feature_ty"])] <- temp.cur.HF[y, "area"]
                                
                        }
                        
                }
                
                
                
                
        }
        
        
        # Current Native
        if(nrow(temp.cur.native) != 0) {
                
                if(as.numeric(table(colnames(temp.cur.native) %in% c("quadrant", "site_year"))["TRUE"]) == 2) {
                        
                        temp.cur.native <- aggregate(temp.cur.native[, "area"],
                                                     by = list(temp.cur.native[, "quadrant"], 
                                                               temp.cur.native[, "site_year"],
                                                               temp.cur.native[, landscape.ty]),
                                                     FUN = sum)
                        
                        colnames(temp.cur.native) <- c("quadrant", "site_year", landscape.ty, "area")
                        temp.cur.native["site_id"] <- paste(temp.cur.native$site_year, temp.cur.native$quadrant, sep = "_")
                        
                        for (y in 1:nrow(temp.cur.native)) {
                                
                                matrix.curr[as.character(temp.cur.native[y, "site_id"]), as.character(temp.cur.native[y, landscape.ty])] <- temp.cur.native[y, "area"]
                                
                        }
                        
                        
                } else {
                        
                        temp.cur.native <- aggregate(temp.cur.native[, "area"],
                                                     by = list(temp.cur.native[, "site_year"],
                                                               temp.cur.native[, landscape.ty]),
                                                     FUN = sum)
                        
                        colnames(temp.cur.native) <- c("site_year", landscape.ty, "area")
                        temp.cur.native["site_id"] <- temp.cur.native$site_year
                        
                        for (y in 1:nrow(temp.cur.native)) {
                                
                                matrix.curr[as.character(temp.cur.native[y, "site_id"]), as.character(temp.cur.native[y, landscape.ty])] <- temp.cur.native[y, "area"]
                                
                        }
                        
                }
                
        }
        
        # Reference Native
        
        if(table(colnames(longform.matrix) %in% c("quadrant", "site_year"))["TRUE"] == 2) {
                
                longform.matrix <- aggregate(longform.matrix[, "area"],
                                             by = list(longform.matrix[, "quadrant"], 
                                                       longform.matrix[, "site_year"],
                                                       longform.matrix[, landscape.ty]),
                                             FUN = sum)
                
                colnames(longform.matrix) <- c("quadrant", "site_year", landscape.ty, "area")
                longform.matrix["site_id"] <- paste(longform.matrix$site_year, longform.matrix$quadrant, sep = "_")
                
                for (y in 1:nrow(longform.matrix)) {
                        
                        matrix.ref[as.character(longform.matrix[y, "site_id"]), as.character(longform.matrix[y, landscape.ty])] <- longform.matrix[y, "area"]
                        
                }
                
                
        } else {
                
                longform.matrix <- aggregate(longform.matrix[, "area"],
                                             by = list(longform.matrix[, "site_year"],
                                                       longform.matrix[, landscape.ty]),
                                             FUN = sum)
                
                colnames(longform.matrix) <- c("site_year", landscape.ty, "area")
                longform.matrix["site_id"] <- longform.matrix$site_year
                
                for (y in 1:nrow(longform.matrix)) {
                        
                        matrix.ref[as.character(longform.matrix[y, "site_id"]), as.character(longform.matrix[y, landscape.ty])] <- longform.matrix[y, "area"]
                        
                }
                
        }
        
        if (type.proportion == TRUE) {
                
                # Convert to proportions
                matrix.curr <- matrix.curr / rowSums(matrix.curr)
                matrix.ref <- matrix.ref / rowSums(matrix.ref)
                
        }
        
        # # Reclassify feature types to those used in the analysis
        matrix.curr <- landscape_hf_reclass(matrix.curr, feature.reclass)
        matrix.ref <- landscape_hf_reclass(matrix.ref, feature.reclass)
        
        matrix.results <- list(matrix.curr, matrix.ref)
        names(matrix.results) <- c("curr", "ref")
        
        return(matrix.results)
        
        rm(temp_cur_HF, temp_cur_native, temp_cur, temp_ref, matrix.results, matrix.curr, matrix.ref)
        
}

#####################
# Soil + HF Reclass # 
#####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

landscape_hf_reclass <- function(data.in, landscape.lookup) {
        
        # Matching of lookup tables and merging native features
        landscape.lookup <- landscape.lookup[landscape.lookup$Class_Fine %in% colnames(data.in), ]
        landscape.clean <- matrix(nrow = nrow(data.in), ncol = length(unique(landscape.lookup[, "Class_Use"])))
        
        for (abmi.coef in 1:length(unique(landscape.lookup[, "Class_Use"]))) {
                
                coef.temp <- as.character(landscape.lookup[landscape.lookup[, "Class_Use"] %in% as.character(unique(landscape.lookup[, "Class_Use"]))[abmi.coef], "Class_Fine"])
                
                if(length(coef.temp) == 1) {
                        
                        landscape.clean[, abmi.coef] <- data.in[, coef.temp]
                        
                } else {
                        
                        landscape.clean[, abmi.coef] <- rowSums(data.in[, coef.temp])
                        
                }
                
        }
        
        colnames(landscape.clean) <- as.character(unique(landscape.lookup[, "Class_Use"]))
        rownames(landscape.clean) <- rownames(data.in)
        
        return(landscape.clean)
        
}


